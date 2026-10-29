@testset "add_pwaformulation and records" begin
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
        z, r = add_pwaformulation(m, xv, f; name = "cost", kwargs...)
        @test r isa PWAFormulation
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
    z, r = add_pwaformulation(m, xv, cvx; z = zv)
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
    @test_throws ArgumentError add_pwaformulation(
        m,
        xv,
        g;
        formulation = LambdaFormulation(),
    )
    @test num_variables(m) == 1
    @test num_constraints(m; count_variable_in_set_constraints = true) == 0
    # ... and the next function can be added
    z = g(xv; bounds = [(-1, 1)])
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
    z, r = add_pwaformulation(m, xv, cost)
    @test r isa PWAFormulation
    @test cost(xv) isa VariableRef
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
    z = cvx(xv)
    @test sides(z) == (true, false)
    @test minimum(cvx(t) for t ∈ grid) - 1e-3 ≤
          lower_bound(z) ≤
          minimum(cvx(t) for t ∈ grid)
    @test sides(cvx(xv; direction = Graph())) == (true, true)
    @test sides(cvx(xv; direction = Hypograph())) == (false, true)
    # Concave: from above in the hypograph direction
    @test sides(ccv(xv)) == (false, true)
    # Difference of convex functions, for all formulations
    for formulation ∈ (PsiFormulation(), DeltaFormulation(), PhiFormulation())
        z = dc(xv; formulation)
        @test sides(z) == (true, true)
        @test lower_bound(z) ≤ minimum(dc(t) for t ∈ grid) + 1e-9
        @test upper_bound(z) ≥ maximum(dc(t) for t ∈ grid) - 1e-9
        z = dc(xv; formulation, direction = Epigraph())
        @test sides(z) == (true, false)
        z = dc(xv; formulation, direction = Hypograph())
        @test sides(z) == (false, true)
    end
    # The 1D formulations bound z exactly
    z = dc(xv; formulation = DeltaFormulation())
    @test minimum(dc(t) for t ∈ grid) - 1e-3 ≤
          lower_bound(z) ≤
          minimum(dc(t) for t ∈ grid)

    # Without bounds on x, the Ψ-formulation of a convex function adds no
    # bounds, as x may leave the domain
    m = Model()
    @variable(m, xv)
    z = cvx(xv)
    @test sides(z) == (false, false)
    # ... while a difference of convex function restricts x to the domain used
    # for the big-M values, and bounds z
    z = dc(xv; formulation = PsiFormulation())
    @test sides(z) == (true, true)
    @test !isnothing(constraint_by_name(m, "pwa_2_bounds[1]"))
    m = Model(optimizer)
    @variable(m, xv)
    z = dc(xv)
    # (the default Δ-formulation restricts x to the breakpoints)
    @test any(v -> startswith(name(v), "_Δ_"), all_variables(m))
    @constraint(m, xv ≥ 7)
    optimize!(m)
    @test termination_status(m) ∈ (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)

    # In two dimensions, the bounds are valid
    s(x) = sin(2x[1]) * x[2]
    f2 = PWA._sample_uniform(s, [(-1, 1), (-1, 1)], 7)
    dc2 = approx(f2, DC(), Cluster(planes = 3, rng = StableRNG(1)))
    m = Model()
    @variable(m, -1 ≤ y[1:2] ≤ 1)
    z = dc2(y)
    pts = [
        (a, b) for a ∈ range(-1, 1; length = 41) for
        b ∈ range(-1, 1; length = 41)
    ]
    @test lower_bound(z) ≤ minimum(dc2(p) for p ∈ pts) + 1e-9
    @test upper_bound(z) ≥ maximum(dc2(p) for p ∈ pts) - 1e-9
end

@testset "Default formulations" begin
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
    s(x) = sin(2x[1]) * x[2]
    dc2 = approx(
        PWA._sample_uniform(s, [(-1, 1), (-1, 1)], 6),
        DC(),
        Cluster(planes = 2, rng = StableRNG(1)),
    )
    uses(m, prefix) = any(v -> startswith(name(v), prefix), all_variables(m))
    # Convex in the natural direction: Ψ, a linear program
    m = Model()
    @variable(m, 0 ≤ xv ≤ 2π)
    cvx(xv)
    @test !uses(m, "_Δ_") && count(is_binary, all_variables(m)) == 0
    # One variable with binaries: Δ, with the domain, explicit bounds or the
    # bounds of the variable
    for (f, kw) ∈ ((dc, (;)), (cvx, (direction = Graph(),)))
        m = Model()
        @variable(m, xv)
        f(xv; kw...)
        @test uses(m, "_Δ_")
    end
    g = PWA._with_domain(dc, nothing)
    m = Model()
    @variable(m, 0 ≤ xv ≤ 2π)
    g(xv)
    @test uses(m, "_Δ_")
    m = Model()
    @variable(m, xv)
    g(xv; bounds = [(0, 2π)])
    @test uses(m, "_Δ_")
    # ... and Ψ without any bounds, or in two dimensions
    m = Model()
    @variable(m, xv)
    g(xv; selection = :indicator)
    @test !uses(m, "_Δ_") && count(is_binary, all_variables(m)) == 6
    m = Model()
    @variable(m, -1 ≤ y[1:2] ≤ 1)
    dc2(y)
    @test !uses(m, "_Δ_")
    # The direction types are not exported
    @test :Graph ∉ names(PiecewiseAffineApprox)
end

@testset "pwaffine is deprecated" begin
    x = collect(range(0, 2π; length = 40))
    data = FunctionEvaluations(x, (x .- 3) .^ 2)
    cvx = approx(data, Convex(), Interpol(planes = 3))
    m = Model()
    @variable(m, 0 ≤ xv ≤ 2π)
    z = @test_deprecated pwaffine(m, xv, cvx)
    @test z isa VariableRef
    z = @test_deprecated pwaffine(m, xv, data, Convex(), Interpol(planes = 3))
    @test z isa VariableRef
    # Constant inputs with a given output add the function to its model
    m = Model(optimizer)
    @variable(m, zv)
    z = cvx(1.5; z = zv)
    @test z === zv
    @objective(m, Min, zv)
    optimize!(m)
    @test value(zv) ≈ cvx(1.5) atol = 1e-6
    @test_throws ArgumentError cvx(@variable(Model()); z = zv)
end

# A formulation method added outside the package: the Ψ-formulation of a
# convex function as one constraint per plane, in a single container
struct VectorPsi <: PWAMethod end
function PWA.formulate!(
    m,
    x,
    z,
    f::PWAFunc{Convex,D},
    ::VectorPsi,
    ::Epigraph;
    name,
    kwargs...,
) where {D}
    P = f.planes
    PWA.record!(
        m,
        @constraint(
            m,
            [k ∈ eachindex(P)],
            z ≥ sum(P[k].α[i] * x[i] for i ∈ 1:D) + P[k].β,
            base_name = "vpsi_$(name)"
        )
    )
    return nothing
end

# ... and one that forgets to record a constraint
struct Unrecorded <: PWAMethod end
function PWA.formulate!(m, x, z, f::PWAFunc, ::Unrecorded, ::Direction; kwargs...)
    @constraint(m, z ≥ 0)
    return nothing
end

@testset "Adding a formulation method" begin
    x = collect(range(0, 2π; length = 40))
    cvx = approx(FunctionEvaluations(x, (x .- 3) .^ 2), Convex(), Interpol(planes = 3))
    dc = approx(FunctionEvaluations(x, sin.(x)), DC(), Interpol(planes = 3; optimizer))
    ncons(m) = num_constraints(m; count_variable_in_set_constraints = false)

    # The new method is used like the built-in ones, and is recorded
    m = Model(optimizer)
    @variable(m, 0 ≤ xv ≤ 2π)
    z, r = @test_logs add_pwaformulation(m, xv, cvx; formulation = VectorPsi())
    @test length(r.constraints) == length(cvx.planes)
    # the output is bounded over the bounds of the variable
    @test has_lower_bound(z)
    @objective(m, Min, z)
    @constraint(m, c, xv == 1.5)
    optimize!(m)
    @test value(z) ≈ cvx(1.5) atol = 1e-6
    delete(m, c)
    delete(m, r)
    @test num_variables(m) == 1 && ncons(m) == 0

    # Without a method for the function, formulation and direction, the model
    # is left unchanged
    for (f, kw) ∈ ((dc, (;)), (cvx, (direction = Graph(),)))
        m = Model()
        @variable(m, 0 ≤ xv ≤ 2π)
        @test_throws ArgumentError add_pwaformulation(m, xv, f; formulation = VectorPsi(), kw...)
        @test num_variables(m) == 1 && ncons(m) == 0
    end
    m = Model()
    @variable(m, -1 ≤ y[1:2] ≤ 1)
    @test_throws ArgumentError cvx(y)

    # Unrecorded variables or constraints give a warning
    m = Model()
    @variable(m, 0 ≤ xv ≤ 2π)
    @test_logs (:warn, r"not recorded") add_pwaformulation(m, xv, cvx; formulation = Unrecorded())

    # The built-in formulations record everything they add
    s(x) = sin(2x[1]) * x[2]
    dc2 = approx(PWA._sample_uniform(s, [(-1, 1), (-1, 1)], 6), DC(), Cluster(planes = 2, rng = StableRNG(1)))
    for (f, kw) ∈ (
        (cvx, (;)),
        (cvx, (bounds = [(0, 2π)],)),
        (cvx, (formulation = LambdaFormulation(),)),
        (cvx, (formulation = DeltaFormulation(),)),
        (cvx, (formulation = PhiFormulation(),)),
        (cvx, (direction = Graph(),)),
        (cvx, (direction = Hypograph(), formulation = PsiFormulation())),
        (dc, (;)),
        (dc, (formulation = PsiFormulation(),)),
        (dc, (formulation = PsiFormulation(), selection = :indicator)),
        (dc, (formulation = PhiFormulation(), direction = Epigraph())),
        (dc, (formulation = DeltaFormulation(), direction = Hypograph())),
    )
        m = Model()
        @variable(m, 0 ≤ xv ≤ 2π)
        @test_logs add_pwaformulation(m, xv, f; kw...)
        @variable(m, zv)
        @test_logs add_pwaformulation(m, xv, f; z = zv, kw...)
    end
    m = Model()
    @variable(m, -1 ≤ y[1:2] ≤ 1)
    @test_logs add_pwaformulation(m, y, dc2)
end
