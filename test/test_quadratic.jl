using Clarabel

@testset "Quadratic pieces" begin
    q = Quadratic([2.0 1.0; 0.0 2.0], (1.0, -1.0), 0.5)
    @test q.Q == [2.0 0.5; 0.5 2.0]                 # symmetrized
    @test evaluate(q, (1.0, 2.0)) ≈ 2 + 2 * 0.5 * 2 + 8 + 1 - 2 + 0.5
    @test_throws ArgumentError Quadratic([1.0 0.0; 0.0 -1.0], (0.0, 0.0), 0.0)
    @test_throws ArgumentError Quadratic([1.0;;], (0.0, 0.0), 0.0)
    @test Quadratic(zeros(1, 1), (2,), 1) == Quadratic([0.0;;], (2.0,), 1.0)

    # The maximum of convex pieces, and the minimum of concave ones
    p = Quadratic([1.0;;], (0.0,), 0.0)              # x²
    r = Quadratic([0.0;;], (1.0,), 1.0)              # x + 1
    f = PWQFunc{Convex,1}([p, r]; domain = [(-2, 3)])
    @test evaluate(f, (0.5,)) ≈ 1.5
    @test evaluate(f, (3.0,)) ≈ 9.0
    @test f(0.5) ≈ 1.5
    g = PWQFunc{Concave,1}([p, r])
    @test evaluate(g, (0.5,)) ≈ -1.5
    @test domain(f) == [(-2.0, 3.0)]
    @test f == PWQFunc{Convex,1}([p, r]; domain = [(-2, 3)])
    @test occursin("2 quadratic pieces", sprint(show, MIME("text/plain"), f))
    @test_throws ArgumentError PWQFunc{DC,1}([p])
end

@testset "Cluster with quadratic pieces" begin
    pts = [(a, b) for a ∈ range(-1, 2, length = 9) for b ∈ range(0, 1, length = 7)]
    fe = FunctionEvaluations(pts, [1 + a^2 + 2b^2 + a for (a, b) ∈ pts])

    # A quadratic is recovered by one piece, with its curvature
    f = approx(fe, Convex(), Cluster(; planes = 1, pieces = :quadratic, trials = 2, rng = StableRNG(1)))
    @test f isa PWQFunc{Convex,2}
    @test approximation_error(f, fe, :max) < 1e-6
    @test only(f.pieces).Q ≈ [1.0 0.0; 0.0 2.0] atol = 1e-6
    @test domain(f) == [(-1.0, 2.0), (0.0, 1.0)]

    # Concave by negation
    g = approx(FunctionEvaluations(pts, -fe.values), Concave(), Cluster(; planes = 2, pieces = :quadratic, trials = 2, rng = StableRNG(1)))
    @test g isa PWQFunc{Concave,2}
    @test approximation_error(g, FunctionEvaluations(pts, -fe.values), :max) < 1e-6

    # One-sided at the data, with an optimizer
    for strict ∈ (:inner, :outer)
        alg = Cluster(; planes = 3, pieces = :quadratic, strict, metric = :l1, optimizer, trials = 2, rng = StableRNG(1))
        f = approx(fe, Convex(), alg)
        dev = [evaluate(f, x) - v for (x, v) ∈ fe]
        @test strict == :inner ? all(dev .≥ -1e-6) : all(dev .≤ 1e-6)
    end

    # A smooth function needs fewer quadratic pieces than planes
    pts = [(a, b) for a ∈ range(-1, 1.5, length = 20) for b ∈ range(-1, 1, length = 20)]
    h = FunctionEvaluations(pts, [exp(a) + b^2 + 0.5a * b + 0.2a^4 for (a, b) ∈ pts])
    rms(pieces, K) = approximation_error(approx(h, Convex(), Cluster(; planes = K, pieces, trials = 4, rng = StableRNG(1))), h, :rms)
    @test rms(:quadratic, 3) < rms(:affine, 8)
    @test rms(:quadratic, 3) < 0.25 * rms(:affine, 3)

    # Constrained least squares (Lawson & Hanson), against a direct check:
    # min (θ₁ - 1)² + (θ₂ - 2)² with θ₁ + θ₂ ≤ 1 gives θ = (0, 1)
    θ = PWA._lsi([1.0 0.0; 0.0 1.0], [1.0, 2.0], [-1.0 -1.0], [-1.0])
    @test θ ≈ [0.0, 1.0] atol = 1e-6

    # Not available
    q = Cluster(; planes = 2, pieces = :quadratic)
    @test_throws ArgumentError approx(fe, DC(), q)
    @test_throws ArgumentError approx(fe, DC(), Alternating(convex = q, concave = Cluster(planes = 2)))
    @test_throws ArgumentError fit((x, y) -> x^2 + y^2, [(0, 1), (0, 1)], Convex(), q; everywhere = true)
    @test_throws ArgumentError approx(fe, Convex(), Cluster(; pieces = :cubic))
end

@testset "Quadratic pieces in a model" begin
    # max(x², x + 1) and its concave mirror image
    p = Quadratic([1.0;;], (0.0,), 0.0)
    r = Quadratic([0.0;;], (1.0,), 1.0)
    f = PWQFunc{Convex,1}([p, r])
    g = PWQFunc{Concave,1}([p, r])
    solver = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)
    for (x0, expected) ∈ ((0.5, 1.5), (2.5, 6.25), (-1.5, 2.25))
        m = Model(solver)
        @variable(m, x == x0)
        z, record = add_pwaformulation(m, x, f; name = "f")
        @test length(record.constraints) == 2
        @objective(m, Min, z)
        optimize!(m)
        @test value(z) ≈ expected atol = 1e-5

        m = Model(solver)
        @variable(m, x == x0)
        @objective(m, Max, g(x))
        optimize!(m)
        @test objective_value(m) ≈ -expected atol = 1e-5
    end

    # Rolled back, with the bounds on the inputs
    m = Model()
    @variable(m, x)
    z, record = add_pwaformulation(m, x, f; bounds = [(-1, 2)])
    @test num_constraints(m; count_variable_in_set_constraints = false) == 3
    delete(m, record)
    @test num_constraints(m; count_variable_in_set_constraints = false) == 0

    # Only the Ψ-formulation in the natural direction
    @test_throws ArgumentError f(x; formulation = LambdaFormulation())
    @test_throws ArgumentError f(x; direction = PWA.Graph())
    @test_throws ArgumentError g(x; direction = PWA.Epigraph())
    @test_throws ArgumentError add_pwaformulation(m, [x, x], f)
    @test num_constraints(m; count_variable_in_set_constraints = false) == 0
end

# The Weymouth equation for the flow in a pipe, with a hydrogen fraction x₃,
# approximated by planes and by quadratic pieces in a small gas network
const M₁, M₂, M₃, Cw = 16.042, 2.016, 28.96, 0.1375
weymouth(p₁, p₂, x₃) = sqrt(Cw * (p₁^2 - p₂^2) * M₃ / (M₁ * (1 - x₃) + M₂ * x₃))

@testset "Gas network with a Weymouth approximation" begin
    X = [(p₁, p₂, x₃) for p₁ ∈ 145:2.5:175 for p₂ ∈ 145:2.5:175 for x₃ ∈ 0:0.05:0.2 if p₁ > p₂]
    data = FunctionEvaluations(X, [weymouth(x...) for x ∈ X])
    # Source 1, junction 2 and sinks 3 and 4, with a capacity factor for each
    # pipe, at a hydrogen fraction of 0.1. With a fixed fraction, the exact
    # capacity q ≤ K κ √(pᵢ² - pⱼ²) is a second order cone.
    pipes = [(1, 2, 2.0), (2, 3, 1.0), (2, 4, 0.8)]
    H₂ = 0.1
    κ = weymouth(1.0, 0.0, H₂)
    solver = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => true)
    function delivered(f, opt)
        m = Model(opt)
        @variable(m, 145 ≤ p[1:4] ≤ 175)
        @variable(m, q[1:3] ≥ 0)
        @variable(m, x₃ == H₂)
        for (k, (i, j, K)) ∈ enumerate(pipes)
            @constraint(m, p[i] ≥ p[j] + 2.5)
            if isnothing(f)
                @constraint(m, [p[i], q[k] / (K * κ), p[j]] ∈ SecondOrderCone())
            else
                @constraint(m, q[k] ≤ K * f([p[i], p[j], x₃]))
            end
        end
        @constraint(m, q[1] == q[2] + q[3])
        @constraint(m, q[3] ≤ 25)
        @objective(m, Max, q[2] + q[3])
        optimize!(m)
        @test is_solved_and_feasible(m)
        return objective_value(m)
    end
    exact = delivered(nothing, solver)
    @test exact ≈ 64.6 atol = 0.1

    flow_error(pieces, K; kwargs...) = begin
        alg = Cluster(; planes = K, pieces, trials = 4, rng = StableRNG(1), optimizer, kwargs...)
        f = approx(data, Concave(), alg)
        (delivered(f, pieces == :quadratic ? solver : optimizer) - exact) / exact, f
    end
    # With the least squares fit, the delivered flow is close to the exact
    e, f = flow_error(:quadratic, 3)
    @test f isa PWQFunc{Concave,3}
    @test abs(e) < 0.005
    @test abs(first(flow_error(:affine, 3))) < 0.01

    # Conservative fits, below the capacity at the data: one quadratic piece
    # is close, one plane is far off
    e, f = flow_error(:quadratic, 1; strict = :inner, metric = :l1)
    @test all(evaluate(f, x) ≤ v + 1e-6 for (x, v) ∈ data)
    @test -0.03 < e < 0
    @test first(flow_error(:affine, 1; strict = :inner, metric = :l1)) < -0.05
end
