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
    err(fe, pwa) = PWA._tolerance_error(fe, pwa, :l1)
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
