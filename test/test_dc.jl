@testset "DC types" begin
    g = PWAFunc{Convex,1}([Plane((1,), 0), Plane((-1,), 0)])
    h = PWAFunc{Concave,1}([Plane((0.5,), -1.0)])
    f = DCFunc(g, h)
    @test f isa DCFunc{1,Float64}
    @test (@inferred evaluate(f, (2.0,))) == 2.0 + 0.0
    @test evaluate(f, (-2,)) == 2.0 + 2.0
    @test PWA._planes(f) == 3
    @test summary(f) == "DCFunc{1, Float64} with 2 convex and 1 concave planes"
    @test repr(MIME("text/plain"), f) ==
          "DCFunc{1, Float64} with 2 convex and 1 concave planes:\n" *
          " z = z₁ + z₂\n z₁ ≥ 1.0 x₁ + 0.0\n z₁ ≥ -1.0 x₁ + 0.0\n" *
          " z₂ ≤ -0.5 x₁ + 1.0\n"

    # An empty part is zero
    e = DCFunc(
        PWAFunc{Convex,2}(),
        PWAFunc{Concave,2}([Plane((1.0, 2.0), 1.0)]),
    )
    @test (@inferred evaluate(e, (1.0, 1.0))) == -4.0
    @test evaluate(DCFunc(PWAFunc{Convex,2}(), PWAFunc{Concave,2}()), (1, 1)) ==
          0

    @test_throws ArgumentError approx(
        FunctionEvaluations([(0.0,), (1.0,)], [0.0, 1.0]),
        Convex(),
        Alternating(convex = Cluster()),
    )

    # Alternating is printed in keyword form
    a = Alternating(
        convex = Interpol(planes = 3),
        concave = Interpol(planes = 2),
    )
    @test repr(a) ==
          "Alternating(convex = $(repr(Interpol(planes = 3))), " *
          "concave = $(repr(Interpol(planes = 2))), iterations = 5)"
    @test occursin(
        "optimizer = ",
        repr(Alternating(convex = Cluster(); optimizer)),
    )
end

@testset "DC serialization" begin
    f = DCFunc(
        PWAFunc{Convex,2}([Plane((1.0, 2.0), 1.0), Plane((-3.0, 1.0), -2.0)]),
        PWAFunc{Concave,2}([Plane((0.5, 0.5), 0.0)]),
    )
    tmp = mktempdir()
    fn = joinpath(tmp, "dc3.json")
    JSON3.write(fn, f)
    @test JSON3.read(fn, DCFunc) == f
    fn = joinpath(tmp, "dc.json")
    JSON.json(fn, f)
    @test JSON.parsefile(fn, DCFunc) == f
end

@testset "DC decomposition" begin
    x = collect(range(-1, 1; length = 25))
    z = sin.(3x)
    f = FunctionEvaluations(tuple.(x), z)
    slopes(v) = diff(v) ./ diff(x)
    # Closed form in one dimension, also for unsorted points
    perm = shuffle(StableRNG(2), collect(eachindex(x)))
    for fe ∈ (f, FunctionEvaluations(tuple.(x[perm]), z[perm]))
        parts = dc_decomposition(fe)
        @test parts.convex.points === fe.points
        @test parts.convex.values .+ parts.concave.values ≈ fe.values
        o = sortperm(first.(fe.points))
        @test minimum(diff(slopes(parts.convex.values[o]))) ≥ -1e-12
        @test maximum(diff(slopes(parts.concave.values[o]))) ≤ 1e-12
        # The concave part has no affine component
        @test abs(sum(parts.concave.values)) ≤ 1e-12
        @test abs(sum(first.(fe.points) .* parts.concave.values)) ≤ 1e-12
    end
    # The linear program gives the same decomposition up to its affine part
    lp = invoke(
        dc_decomposition,
        Tuple{FunctionEvaluations{1},Any},
        f,
        optimizer,
    )
    parts = dc_decomposition(f)
    @test lp.convex.values .+ lp.concave.values ≈ z
    @test diff(slopes(lp.concave.values)) ≈ diff(slopes(parts.concave.values)) atol =
        1e-6

    # Convex data has no concave part
    parts = dc_decomposition(FunctionEvaluations(tuple.(x), x .^ 2))
    @test all(iszero, parts.concave.values)
    parts = dc_decomposition(FunctionEvaluations(tuple.(x), x .^ 2), optimizer)
    @test all(iszero, parts.concave.values)

    # Two dimensions: the parts can be interpolated by convex and concave
    # functions, i.e. enforce_curvature does not change them
    s(x) = x[1]^2 - x[2]^2
    f2 = PWA._sample_uniform(s, [(-1, 1), (-1, 1)], 8)
    for metric ∈ (:l1, :max)
        parts = dc_decomposition(f2, optimizer; metric)
        @test parts.convex.values .+ parts.concave.values ≈ f2.values
        @test enforce_curvature(parts.convex, Convex(), optimizer).values ≈
              parts.convex.values
        @test enforce_curvature(parts.concave, Concave(), optimizer).values ≈
              parts.concave.values
    end
    @test_throws ArgumentError dc_decomposition(f2)
    @test_throws ArgumentError dc_decomposition(f2, optimizer; metric = :none)
    # Repeated points with different values have no decomposition
    @test_throws ErrorException dc_decomposition(
        FunctionEvaluations([(0.0,), (0.0,), (1.0,)], [0.0, 1.0, 0.0]),
        optimizer,
    )
end

@testset "DC approximation" begin
    x = collect(range(-1, 1; length = 25))
    f = FunctionEvaluations(tuple.(x), sin.(3x))
    err(fe, pwa) = approximation_error(pwa, fe, :l1)
    for alg ∈ (
        Cluster(planes = 4, rng = StableRNG(1)),
        Interpol(planes = 4, optimizer = optimizer),
        Interpol(planes = 4),
        MILP(planes = 3, optimizer = optimizer),
    )
        dc = approx(f, DC(), alg)
        @test dc isa DCFunc{1,Float64}
        @test PWA._planes(dc.convex) ≤ alg.planes
        @test PWA._planes(dc.concave) ≤ alg.planes
        # Better than a convex and a concave approximation
        @test err(f, dc) <
              0.5 * min(
            err(f, approx(f, Convex(), alg)),
            err(f, approx(f, Concave(), alg)),
        )
    end
    # Different algorithms and numbers of planes for the two parts
    alg = Alternating(
        convex = Cluster(planes = 5, rng = StableRNG(1)),
        concave = Interpol(planes = 2),
    )
    dc = approx(f, DC(), alg)
    @test PWA._planes(dc.convex) == 5
    @test PWA._planes(dc.concave) == 2
    @test err(f, dc) < 0.2
    @test alg.iterations == 5
    @test Alternating(convex = Progressive(; optimizer)).iterations == 0

    # Algorithms requiring convex data are used on the decomposed parts
    dc = approx(f, DC(), Progressive(; optimizer, tolerance = 0.05))
    @test err(f, dc) ≤ 0.1
    dc = approx(f, DC(), FullOrder(; optimizer))
    @test err(f, dc) ≤ 1e-6

    # Sampling a function
    dc = approx(x -> sin(3x[1]), [(-1, 1)], DC(), Cluster(planes = 4))
    @test err(f, dc) < 0.1
    dc = approx(x -> sin(3x[1]), [(-1, 1)], DC(), alg)
    @test err(f, dc) < 0.2

    # Two dimensions
    s(x) = x[1]^2 - x[2]^2
    f2 = PWA._sample_uniform(s, [(-1, 1), (-1, 1)], 8)
    for alg ∈ (
        Cluster(planes = 3, rng = StableRNG(1)),
        Cluster(planes = 3, optimizer = optimizer, rng = StableRNG(1)),
    )
        dc = approx(f2, DC(), alg)
        @test err(f2, dc) < 0.25 * err(f2, approx(f2, Convex(), alg))
    end
end

@testset "DC in JuMP" begin
    x = collect(range(0, 2π; length = 40))
    f = FunctionEvaluations(tuple.(x), sin.(x))
    dc = approx(f, DC(), Interpol(planes = 3; optimizer))
    @test PWA._planes(dc.convex) == 3 && PWA._planes(dc.concave) == 3

    # The value of z at fixed points, in both directions of optimization
    function solve(dc, x̂, sense; kwargs...)
        m = Model(optimizer)
        @variable(m, 0 ≤ xv ≤ 2π)
        z = pwaffine(m, (xv,), dc; kwargs...)
        @constraint(m, xv == x̂)
        @objective(m, sense, z)
        optimize!(m)
        @test is_solved_and_feasible(m)
        return value(z)
    end
    for x̂ ∈ (0.3, 1.5, π, 4.4, 6.0)
        v = evaluate(dc, (x̂,))
        @test solve(dc, x̂, MOI.MIN_SENSE) ≈ v atol = 1e-6
        @test solve(dc, x̂, MOI.MAX_SENSE) ≈ v atol = 1e-6
        @test solve(dc, x̂, MOI.MIN_SENSE; direction = :epigraph) ≈ v atol = 1e-6
        @test solve(dc, x̂, MOI.MAX_SENSE; direction = :hypograph) ≈ v atol =
            1e-6
    end

    # Binary variables only for the parts that need them
    nbin(; kwargs...) = begin
        m = Model()
        @variable(m, 0 ≤ xv ≤ 2π)
        pwaffine(m, (xv,), dc; kwargs...)
        return count(is_binary, all_variables(m))
    end
    @test nbin() == 6
    @test nbin(direction = :epigraph) == 3
    @test nbin(direction = :hypograph) == 3
    @test nbin(selection = :indicator) == 6
    m = Model()
    @variable(m, xv)
    pwaffine(m, (xv,), dc; selection = :indicator)
    @test num_constraints(m, VariableRef, MOI.ZeroOne) == 6
    @test num_constraints(
        m,
        Vector{AffExpr},
        MOI.Indicator{MOI.ACTIVATE_ON_ONE,MOI.LessThan{Float64}},
    ) == 3

    # Bounds are required for big-M, from the variables or the keyword
    m = Model()
    @variable(m, xv)
    @test_throws ArgumentError pwaffine(m, (xv,), dc)
    @test_throws ArgumentError pwaffine(m, (xv,), dc; bounds = [(0, Inf)])
    @test_throws ArgumentError pwaffine(m, (xv,), dc; bounds = [(0, 1), (0, 1)])
    @test_throws ArgumentError pwaffine(m, (xv,), dc; direction = :none)
    @test_throws ArgumentError pwaffine(
        m,
        (xv,),
        dc;
        formulation = LambdaFormulation(),
        bounds = [(0, 1)],
    )
    @test_throws ArgumentError pwaffine(m, (xv,), dc; selection = :sos)
    @test_throws ArgumentError pwaffine(m, (xv, xv), dc)
    # Expressions as inputs with explicit bounds, and a given z
    m = Model(optimizer)
    @variable(m, 0 ≤ t ≤ π)
    @variable(m, zv)
    z = pwaffine(m, (2t,), dc; z = zv, bounds = [(0, 2π)])
    @test z === zv
    @constraint(m, t == 0.75)
    @objective(m, Min, zv)
    optimize!(m)
    @test value(zv) ≈ evaluate(dc, (1.5,)) atol = 1e-6

    # Parts with a single plane or no planes need no binary variables
    g = PWAFunc{Convex,1}([Plane((1.0,), 0.0)])
    h = PWAFunc{Concave,1}([Plane((0.5,), -1.0), Plane((-1.0,), 0.0)])
    for dc₁ ∈ (
        DCFunc(g, h),
        DCFunc(PWAFunc{Convex,1}(), h),
        DCFunc(g, PWAFunc{Concave,1}()),
    )
        m = Model(optimizer)
        @variable(m, -2 ≤ xv ≤ 2)
        z = pwaffine(m, (xv,), dc₁)
        @test count(is_binary, all_variables(m)) ==
              (isempty(dc₁.concave.planes) ? 0 : 2)
        @constraint(m, xv == 0.5)
        @objective(m, Max, z)
        optimize!(m)
        @test value(z) ≈ evaluate(dc₁, (0.5,)) atol = 1e-6
    end

    # Two dimensions, with the fit from data through pwaffine
    s(x) = x[1]^2 - x[2]^2
    f2 = PWA._sample_uniform(s, [(-1, 1), (-1, 1)], 8)
    m = Model(optimizer)
    @variable(m, -1 ≤ xv[1:2] ≤ 1)
    z = pwaffine(m, xv, f2, DC(), Cluster(planes = 3, rng = StableRNG(1)))
    dc2 = approx(f2, DC(), Cluster(planes = 3, rng = StableRNG(1)))
    for x̂ ∈ ((0.3, -0.7), (-0.9, 0.9), (0.0, 0.0))
        c = @constraint(m, xv .== collect(x̂))
        @objective(m, Min, z)
        optimize!(m)
        @test value(z) ≈ evaluate(dc2, x̂) atol = 1e-6
        delete.(m, c)
    end
end

@testset "DC strict" begin
    x = collect(range(-1, 1; length = 25))
    f = FunctionEvaluations(tuple.(x), sin.(3x))
    s(x) = sin(3x[1]) * cos(2x[2])
    f2 = PWA._sample_uniform(s, [(-1, 1), (-1, 1)], 10)
    gap(fe, pwa) = extrema(evaluate(pwa, p) - v for (p, v) ∈ fe)
    err(fe, pwa) = approximation_error(pwa, fe, :l1)

    # strict applies to the sum of the parts, below (:outer) or above (:inner)
    # the data, and is interchanged for the concave part
    @test PWA._dc_swap_strict(Cluster(strict = :outer)).strict == :inner
    @test PWA._dc_swap_strict(Cluster(strict = :inner)).strict == :outer
    @test PWA._dc_swap_strict(Cluster(strict = :none)).strict == :none
    @test PWA._dc_swap_strict(Interpol()) === Interpol()
    # The concave part has concave_planes planes
    @test PWA._dc_concave_algorithm(Cluster(planes = 5)).planes == 5
    @test PWA._dc_concave_algorithm(Cluster(planes = 5, concave_planes = 2)).planes ==
          2
    a = PWA._dc_concave_algorithm(
        MILP(planes = 5, concave_planes = 2, strict = :outer; optimizer),
    )
    @test a == MILP(planes = 2, concave_planes = 2, strict = :inner; optimizer)
    @test PWA._dc_concave_algorithm(FullOrder(; optimizer)) ===
          FullOrder(; optimizer)
    for alg ∈ (
        Cluster(planes = 5, concave_planes = 2, rng = StableRNG(1)),
        Interpol(planes = 5, concave_planes = 2),
        MILP(planes = 3, concave_planes = 2; optimizer),
    )
        dc = approx(f, DC(), alg)
        @test PWA._planes(dc.convex) == alg.planes
        @test PWA._planes(dc.concave) == 2
        @test err(f, dc) < 0.15
    end
    a = MILP(strict = :outer; optimizer, planes = 3)
    @test PWA._dc_swap_strict(a) == MILP(strict = :inner; optimizer, planes = 3)

    for (fe, alg) ∈ (
        (
            f,
            Cluster(
                planes = 4,
                metric = :l1,
                strict = :outer;
                optimizer,
                rng = StableRNG(1),
            ),
        ),
        (
            f,
            Cluster(
                planes = 4,
                metric = :l1,
                strict = :inner;
                optimizer,
                rng = StableRNG(1),
            ),
        ),
        (
            f,
            Cluster(
                planes = 4,
                metric = :max,
                strict = :outer;
                optimizer,
                rng = StableRNG(1),
            ),
        ),
        (
            f2,
            Cluster(
                planes = 4,
                metric = :l1,
                strict = :outer;
                optimizer,
                rng = StableRNG(1),
            ),
        ),
        (
            f2,
            Cluster(
                planes = 4,
                metric = :l1,
                strict = :inner;
                optimizer,
                rng = StableRNG(1),
            ),
        ),
        (f, MILP(planes = 2, strict = :outer; optimizer)),
        (f, MILP(planes = 2, strict = :inner; optimizer)),
        (f, Progressive(; optimizer, tolerance = 0.1, strict = :outer)),
        (f, Progressive(; optimizer, tolerance = 0.1, strict = :inner)),
    )
        dc = approx(fe, DC(), alg)
        lo, hi = gap(fe, dc)
        if alg.strict == :outer
            @test hi ≤ 1e-6
        else
            @test lo ≥ -1e-6
        end
        @test err(fe, dc) < 0.3
    end
    # Also with different algorithms for the parts
    alg = Alternating(
        convex = Cluster(
            planes = 4,
            metric = :l1,
            strict = :outer;
            optimizer,
            rng = StableRNG(1),
        ),
        concave = MILP(planes = 2, strict = :outer; optimizer),
    )
    dc = approx(f, DC(), alg)
    @test gap(f, dc)[2] ≤ 1e-6
    @test PWA._planes(dc.convex) == 4 && PWA._planes(dc.concave) == 2
end

@testset "DC 1D formulations" begin
    x = collect(range(0, 2π; length = 40))
    f = FunctionEvaluations(x, sin.(x))
    dc = approx(f, DC(), Interpol(planes = 3; optimizer))
    @test PWA._planes(dc.convex) == 3 && PWA._planes(dc.concave) == 3

    # The breakpoints over the bounds are the kinks of both parts
    xp, zp = PWA._dc_breakpoints(dc, [(0, 2π)])
    @test xp[1] == 0 && xp[end] ≈ 2π
    @test length(xp) == 6
    @test zp ≈ [evaluate(dc, (xᵢ,)) for xᵢ ∈ xp]
    slope = diff(zp) ./ diff(xp)
    @test length(PWA._dc_binary_kinks(slope, :graph)) == 4
    @test length(PWA._dc_binary_kinks(slope, :epigraph)) == 2
    @test length(PWA._dc_binary_kinks(slope, :hypograph)) == 2
    # Narrower bounds keep only the kinks inside
    xp₂, _ = PWA._dc_breakpoints(dc, [(1, 2)])
    @test xp₂[1] == 1 && xp₂[end] == 2 && length(xp₂) < length(xp)

    function solve(x̂, sense; kwargs...)
        m = Model(optimizer)
        @variable(m, 0 ≤ xv ≤ 2π)
        z = pwaffine(m, xv, dc; kwargs...)
        @constraint(m, xv == x̂)
        @objective(m, sense, z)
        optimize!(m)
        @test is_solved_and_feasible(m)
        return value(z)
    end
    nbin(; kwargs...) = begin
        m = Model()
        @variable(m, 0 ≤ xv ≤ 2π)
        pwaffine(m, xv, dc; kwargs...)
        return count(is_binary, all_variables(m))
    end
    for formulation ∈ (DeltaFormulation(), PhiFormulation())
        # Exact values in both directions of optimization
        for x̂ ∈ (0.0, 0.3, 1.5, π, 4.4, 6.0, 2π)
            v = evaluate(dc, (x̂,))
            @test solve(x̂, MOI.MIN_SENSE; formulation) ≈ v atol = 1e-6
            @test solve(x̂, MOI.MAX_SENSE; formulation) ≈ v atol = 1e-6
            @test solve(x̂, MOI.MIN_SENSE; formulation, direction = :epigraph) ≈
                  v atol = 1e-6
            @test solve(x̂, MOI.MAX_SENSE; formulation, direction = :hypograph) ≈
                  v atol = 1e-6
        end
        # Binary variables only at the kinks against the direction
        @test nbin(; formulation) == 4
        @test nbin(; formulation, direction = :epigraph) == 2
        @test nbin(; formulation, direction = :hypograph) == 2
        # Free optimization finds the minimum and maximum of the function
        m = Model(optimizer)
        @variable(m, 0 ≤ xv ≤ 2π)
        z = pwaffine(m, xv, dc; formulation, direction = :epigraph)
        @objective(m, Min, z)
        optimize!(m)
        @test value(z) ≈ minimum(zp) atol = 1e-6
        @objective(m, Max, z)
        @test_throws ArgumentError pwaffine(
            m,
            xv,
            dc;
            formulation = LambdaFormulation(),
        )
    end
    # The relaxation of the Δ-formulation never undercuts the function: a
    # cheap segment after a concave kink cannot be reached without filling the
    # earlier segments. With z only bounded below, every feasible z is at
    # least f(x).
    g = DCFunc(
        PWAFunc{Convex,1}([Plane((5.0,), 0.0), Plane((6.0,), -1.0)]),
        PWAFunc{Concave,1}([Plane((0.0,), 0.0), Plane((4.0,), -8.0)]),
    )
    xp, zp = PWA._dc_breakpoints(g, [(0, 4)])
    for formulation ∈ (DeltaFormulation(), PhiFormulation())
        m = Model(optimizer)
        @variable(m, 0 ≤ xv ≤ 4)
        z = pwaffine(m, xv, g; formulation, direction = :epigraph)
        for x̂ ∈ range(0, 4; length = 9)
            c = @constraint(m, xv == x̂)
            @objective(m, Min, z)
            optimize!(m)
            @test value(z) ≈ evaluate(g, (x̂,)) atol = 1e-6
            delete(m, c)
        end
    end
    # Two dimensions are not supported by the one-dimensional formulations
    s(x) = x[1]^2 - x[2]^2
    f2 = PWA._sample_uniform(s, [(-1, 1), (-1, 1)], 6)
    dc2 = approx(f2, DC(), Cluster(planes = 2, rng = StableRNG(1)))
    m = Model()
    @variable(m, -1 ≤ xv[1:2] ≤ 1)
    @test_throws ArgumentError pwaffine(
        m,
        xv,
        dc2;
        formulation = DeltaFormulation(),
    )
    @test_throws ArgumentError pwaffine(
        m,
        xv,
        dc2;
        formulation = PhiFormulation(),
    )
end
