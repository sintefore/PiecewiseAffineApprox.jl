@testset "add_pwaffine and records" begin
    x = collect(range(0, 2π; length = 40))
    cvx = approx(
        FunctionEvaluations(x, (x .- 3) .^ 2),
        Convex(),
        Interpol(planes = 3),
    )
    dc = approx(
        FunctionEvaluations(x, sin.(x)),
        DC(),
        Interpol(planes = 3; optimizer),
    )
    ncons(m) = num_constraints(m; count_variable_in_set_constraints = false)

    # The record holds what was added, and deleting it restores the model
    for (f, kwargs) ∈ (
        (cvx, (;)),
        (cvx, (formulation = LambdaFormulation(),)),
        (cvx, (direction = Graph(),)),
        (dc, (;)),
        (dc, (formulation = DeltaFormulation(), direction = Epigraph())),
        (dc, (formulation = PhiFormulation(),)),
    )
        m = Model()
        @variable(m, 0 ≤ xv ≤ 2π)
        nv, nc = num_variables(m), ncons(m)
        z, r = add_pwaffine(m, xv, f; name = "cost", kwargs...)
        @test r isa PWARecord
        @test r.output[] === z
        @test z ∈ r.variables
        @test name(z) == "z_cost"
        @test num_variables(m) == nv + length(r.variables)
        @test ncons(m) == nc + length(r.constraints)
        @test occursin("cost", repr(r))
        delete(m, r)
        @test num_variables(m) == nv
        @test ncons(m) == nc
        @test isempty(r.variables) && isempty(r.constraints)
    end
    # A given output is kept when the record is deleted, but its bound
    # constraint is removed
    m = Model()
    @variable(m, 0 ≤ xv ≤ 2π)
    @variable(m, zv)
    z, r = add_pwaffine(m, xv, cvx; z = zv)
    @test z === zv && zv ∉ r.variables
    @test !isnothing(constraint_by_name(m, "pwa_1_zbounds"))
    delete(m, r)
    @test is_valid(m, zv)
    @test ncons(m) == 0
    @test_throws ArgumentError delete(Model(), r)

    # An error leaves the model unchanged: the output is created before the
    # missing bounds are found
    m = Model()
    @variable(m, xv)
    g = PWAFunc{Convex,1}([Plane((1.0,), 0.0), Plane((-1.0,), 0.0)])
    @test_throws ArgumentError add_pwaffine(
        m,
        xv,
        g;
        formulation = LambdaFormulation(),
    )
    @test num_variables(m) == 1
    @test num_constraints(m; count_variable_in_set_constraints = true) == 0
    # ... and the next function can be added
    z = pwaffine(m, xv, g; bounds = [(-1, 1)])
    @test lower_bound(z) == 0.0
end

@testset "Callable functions" begin
    x = collect(range(0, 2π; length = 40))
    cvx = approx(
        FunctionEvaluations(x, (x .- 3) .^ 2),
        Convex(),
        Interpol(planes = 3),
    )
    dc = approx(
        FunctionEvaluations(x, sin.(x)),
        DC(),
        Interpol(planes = 3; optimizer),
    )

    # Numbers are evaluated
    @test cvx(1.0) == evaluate(cvx, (1.0,))
    @test dc(2.0) == evaluate(dc, (2.0,))
    @test cvx((1.0,)) == cvx(1.0)
    @test_throws ArgumentError cvx(1.0; formulation = DeltaFormulation())
    s(x) = x[1]^2 + x[2]^2
    f2 = PWA._sample_uniform(s, [(-1, 1), (-1, 1)], 6)
    p2 = approx(f2, Convex(), Cluster(planes = 4, rng = StableRNG(1)))
    @test p2(0.5, -0.5) == evaluate(p2, (0.5, -0.5)) == p2([0.5, -0.5])

    # JuMP inputs add the function to their model
    m = Model(optimizer)
    @variable(m, 0 ≤ xv[1:3] ≤ 2π)
    @objective(m, Min, sum(dc(xv[i]) for i ∈ 1:3))
    @constraint(m, sum(xv) ≥ 9)
    optimize!(m)
    @test is_solved_and_feasible(m)
    @test objective_value(m) ≈ sum(dc(value(xv[i])) for i ∈ 1:3) atol = 1e-6
    # Keyword arguments, and vector inputs
    m = Model()
    @variable(m, -1 ≤ y[1:2] ≤ 1)
    z = p2(y; name = "p2")
    @test name(z) == "z_p2"
    p2(y[1], y[2]; formulation = LambdaFormulation())
    @test any(v -> startswith(name(v), "_λ_"), all_variables(m))
    @test_throws ArgumentError cvx(y[1], @variable(Model()))

    # Formulated stores the options
    cost =
        Formulated(dc; formulation = DeltaFormulation(), direction = Epigraph())
    @test cost(1.0) == dc(1.0)
    @test evaluate(cost, (1.0,)) == dc(1.0)
    @test domain(cost) == domain(dc)
    @test occursin("formulation = ", repr(cost))
    @test_throws ArgumentError Formulated(dc; foo = 1)
    m = Model(optimizer)
    @variable(m, 0 ≤ xv ≤ 2π)
    z = cost(xv)
    @test any(v -> startswith(name(v), "_Δ_"), all_variables(m))
    @test count(is_binary, all_variables(m)) == 2
    @objective(m, Min, z)
    optimize!(m)
    @test objective_value(m) ≈
          minimum(dc(t) for t ∈ range(0, 2π; length = 20001)) atol = 1e-4
    # Options in the call override the stored ones
    m = Model()
    @variable(m, 0 ≤ xv ≤ 2π)
    cost(xv; direction = Graph())
    @test count(is_binary, all_variables(m)) == 4
    z, r = add_pwaffine(m, xv, cost)
    @test r isa PWARecord
    @test pwaffine(m, xv, cost) isa VariableRef
end

@testset "Output bounds by direction" begin
    x = collect(range(0, 2π; length = 40))
    cvx = approx(
        FunctionEvaluations(x, (x .- 3) .^ 2),
        Convex(),
        Interpol(planes = 3),
    )
    ccv = approx(
        FunctionEvaluations(x, -(x .- 3) .^ 2),
        Concave(),
        Interpol(planes = 3),
    )
    dc = approx(
        FunctionEvaluations(x, sin.(x)),
        DC(),
        Interpol(planes = 3; optimizer),
    )
    sides(z) = (has_lower_bound(z), has_upper_bound(z))
    grid = range(0, 2π; length = 20001)

    # Convex: from below in the epigraph direction, both sides for the graph
    m = Model()
    @variable(m, 0 ≤ xv ≤ 2π)
    z = pwaffine(m, xv, cvx)
    @test sides(z) == (true, false)
    @test minimum(cvx(t) for t ∈ grid) - 1e-3 ≤
          lower_bound(z) ≤
          minimum(cvx(t) for t ∈ grid)
    @test sides(pwaffine(m, xv, cvx; direction = Graph())) == (true, true)
    @test sides(pwaffine(m, xv, cvx; direction = Hypograph())) == (false, true)
    # Concave: from above in the hypograph direction
    @test sides(pwaffine(m, xv, ccv)) == (false, true)
    # Difference of convex functions, for all formulations
    for formulation ∈ (PsiFormulation(), DeltaFormulation(), PhiFormulation())
        z = pwaffine(m, xv, dc; formulation)
        @test sides(z) == (true, true)
        @test lower_bound(z) ≤ minimum(dc(t) for t ∈ grid) + 1e-9
        @test upper_bound(z) ≥ maximum(dc(t) for t ∈ grid) - 1e-9
        z = pwaffine(m, xv, dc; formulation, direction = Epigraph())
        @test sides(z) == (true, false)
        z = pwaffine(m, xv, dc; formulation, direction = Hypograph())
        @test sides(z) == (false, true)
    end
    # The 1D formulations bound z exactly
    z = pwaffine(m, xv, dc; formulation = DeltaFormulation())
    @test minimum(dc(t) for t ∈ grid) - 1e-3 ≤
          lower_bound(z) ≤
          minimum(dc(t) for t ∈ grid)

    # Without bounds on x, the Ψ-formulation of a convex function adds no
    # bounds, as x may leave the domain
    m = Model()
    @variable(m, xv)
    z = pwaffine(m, xv, cvx)
    @test sides(z) == (false, false)
    # ... while a difference of convex function restricts x to the domain used
    # for the big-M values, and bounds z
    z = pwaffine(m, xv, dc)
    @test sides(z) == (true, true)
    @test !isnothing(constraint_by_name(m, "pwa_2_bounds[1]"))
    m = Model(optimizer)
    @variable(m, xv)
    z = pwaffine(m, xv, dc)
    @constraint(m, xv ≥ 7)
    optimize!(m)
    @test termination_status(m) ∈ (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)

    # In two dimensions, the bounds are valid
    s(x) = sin(2x[1]) * x[2]
    f2 = PWA._sample_uniform(s, [(-1, 1), (-1, 1)], 7)
    dc2 = approx(f2, DC(), Cluster(planes = 3, rng = StableRNG(1)))
    m = Model()
    @variable(m, -1 ≤ y[1:2] ≤ 1)
    z = pwaffine(m, y, dc2)
    pts = [
        (a, b) for a ∈ range(-1, 1; length = 41) for
        b ∈ range(-1, 1; length = 41)
    ]
    @test lower_bound(z) ≤ minimum(dc2(p) for p ∈ pts) + 1e-9
    @test upper_bound(z) ≥ maximum(dc2(p) for p ∈ pts) - 1e-9
end
