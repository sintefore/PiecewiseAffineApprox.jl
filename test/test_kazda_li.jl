
@testset "Progressive fitting" begin

    # 1D
    x = collect(range(-1, 1; length = 30))
    z = x .^ 2
    f = FunctionEvaluations(tuple.(x), z)
    for (tol, pen) ∈ [(0.05, :max), (0.02, :l1), (0.03, :l2)]
        pwa_red = approx(
            f,
            Convex(),
            Progressive(optimizer = optimizer, tolerance = tol, metric = pen),
        )
        for (x, z) ∈ f
            pv = evaluate(pwa_red, x)
            @test z ≥ pv || z ≈ pv
        end
        # The tolerance is met for the error per data point
        @test PWA._tolerance_error(f, pwa_red, pen) ≤ tol
        @test length(pwa_red.planes) < length(f)
    end

    # The tolerance is per data point, so the number of planes does not grow
    # with the number of data points
    planes = map((30, 120)) do N
        xN = collect(range(-1, 1; length = N))
        fN = FunctionEvaluations(tuple.(xN), xN .^ 2)
        alg = Progressive(optimizer = optimizer, tolerance = 0.01, metric = :l1)
        return length(approx(fN, Convex(), alg).planes)
    end
    @test planes[2] ≤ planes[1] + 1

    # Metrics without a tolerance interpretation are rejected
    @test_throws ArgumentError approx(
        f,
        Convex(),
        Progressive(optimizer = optimizer, tolerance = 0.5, metric = :none),
    )

    # 2D
    g(x) = x[1]^2 + x[2]^2
    vals = PiecewiseAffineApprox._sample_uniform(g, [(-1, 1), (-1, 1)], 10)
    pwa_red = approx(
        vals,
        Convex(),
        Progressive(optimizer = optimizer, tolerance = 0.2, metric = :max),
    )
    # The number of planes and the value at the origin depend on which of
    # several equally good solutions is found, so test the guaranteed properties
    @test PWA._tolerance_error(vals, pwa_red, :max) ≤ 0.2
    @test all(evaluate(pwa_red, x) ≤ z + 1e-6 for (x, z) ∈ vals)
    @test abs(evaluate(pwa_red, (0, 0))) ≤ 0.025

    @testset "Nonconvex" begin
        h(x) = sin(5 * x[1]) * (x[1]^2 + x[2]^2)
        vals = PiecewiseAffineApprox._sample_uniform(h, [(-1, 1), (-1, 1)], 10)

        @test_throws ErrorException approx(
            vals,
            Convex(),
            Progressive(optimizer = optimizer, tolerance = 2.0, metric = :l2),
        )

        vals_c = enforce_curvature(vals, Convex(), optimizer, :l1)
        @test length(vals_c) == length(vals)

        pwa_con = approx(
            vals_c,
            Convex(),
            Progressive(optimizer = optimizer, tolerance = 0.1, metric = :l2),
        )

        x = collect(range(-1, 1; length = 30))
        z = x .^ 2
        z[5] += 0.1
        f = FunctionEvaluations(tuple.(x), z)

        fc = enforce_curvature(f, Convex(), optimizer)
        @test length(fc) == length(f)
        @test fc.values[5] ≈ 0.5291 atol = 0.001

        pwa_con = approx(
            fc,
            Convex(),
            Progressive(optimizer = optimizer, tolerance = 0.02, metric = :l1),
        )
        @test PWA._tolerance_error(fc, pwa_con, :l1) ≤ 0.02
        @test 1 < PWA._planes(pwa_con) < length(fc)

        # Convex data is returned unchanged, and :none is not a valid metric
        @test enforce_curvature(fc, Convex(), optimizer) === fc
        @test_throws ArgumentError enforce_curvature(f, Convex(), optimizer, :none)

        # The error for non-convex data suggests enforcing curvature
        err = try
            approx(f, Convex(), FullOrder(optimizer = optimizer))
        catch e
            sprint(showerror, e)
        end
        @test occursin("enforce_curvature", err)

        x = collect(range(-1, 1; length = 30))
        z = -x .^ 2
        z[8] -= 0.2
        z[12] += 0.1
        f = FunctionEvaluations(tuple.(x), z)
        fc = enforce_curvature(f, Concave(), optimizer)
        @test fc.values[8] ≈ -0.2722 atol = 0.001
    end
end

@testset "Full Order fitting" begin
    # 1D
    x = collect(range(-1, 1; length = 30))
    z = x .^ 2
    f = FunctionEvaluations(tuple.(x), z)
    for (tol, pen) ∈ [(0.05, :max), (2.5, :l1), (0.5, :l2), (0.5, :none)]
        pwa_red =
            approx(f, Convex(), FullOrder(optimizer = optimizer, metric = pen))
        for (x, z) ∈ f
            pv = evaluate(pwa_red, x)
            @test z ≥ pv || z ≈ pv
        end
    end

    # 2D
    g(x) = x[1]^2 + x[2]^2
    vals = PiecewiseAffineApprox._sample_uniform(g, [(-1, 1), (-1, 1)], 10)
    pwa_red =
        approx(vals, Convex(), FullOrder(optimizer = optimizer, metric = :max))
    @test length(pwa_red.planes) == 100
    @test evaluate(pwa_red, (0, 0)) ≈ 0.024 atol = 0.001

    @testset "Nonconvex" begin
        h(x) = sin(5 * x[1]) * (x[1]^2 + x[2]^2)
        vals = PiecewiseAffineApprox._sample_uniform(h, [(-1, 1), (-1, 1)], 10)

        @test_throws ErrorException approx(
            vals,
            Convex(),
            FullOrder(optimizer = optimizer, metric = :l2),
        )

        vals_c = PWA.convexify(vals, optimizer, :l1)
        @test length(vals_c) == length(vals)

        pwa_con = approx(
            vals_c,
            Convex(),
            FullOrder(optimizer = optimizer, metric = :l2),
        )

        x = collect(range(-1, 1; length = 30))
        z = x .^ 2
        z[5] += 0.1
        f = FunctionEvaluations(tuple.(x), z)

        fc = PWA.convexify(f, optimizer)
        @test length(fc) == length(f)
        @test fc.values[5] ≈ 0.5291 atol = 0.001

        pwa_con =
            approx(fc, Convex(), FullOrder(optimizer = optimizer, metric = :l1))
        @test PWA._planes(pwa_con) == 30
    end
end

@testset "Separable and one-dimensional formulations" begin
    # Lower convex hull: collinear points and duplicates are removed
    hx, hz = PWA._lower_hull([0.0, 1.0, 2.0, 3.0, 1.0, 4.0], [0.0, 1.0, 2.0, 0.0, 3.0, 4.0])
    @test hx == [0.0, 3.0, 4.0]
    @test hz == [0.0, 0.0, 4.0]

    x = sort(2 .* rand(StableRNG(1), 40) .- 1)
    z = exp.(x) .+ 0.02 .* sin.(25 .* x)
    f = FunctionEvaluations(tuple.(x), z)

    # The supporting plane of the lower convex envelope agrees with the
    # optimization problem it replaces: maximize the value at c below all points,
    # for c within the range of the data (as for weighted centroids)
    support = PWA._envelope_support(f, optimizer)
    for c ∈ (x[1], -0.3, 0.0, 0.45, x[end], (x[11] + x[12]) / 2)
        m = Model(optimizer)
        @variable(m, a)
        @variable(m, b)
        @constraint(m, [l ∈ eachindex(x)], a * x[l] + b ≤ z[l])
        @objective(m, Max, a * c + b)
        optimize!(m)
        @test evaluate(support([c]), (c,)) ≈ objective_value(m) atol = 1e-8
    end

    # One-dimensional full order planes are optimal, compared with solving the
    # optimization problem for each plane
    fc = enforce_curvature(f, Convex(), optimizer)
    X = reshape(first.(fc.points), 1, :)
    for metric ∈ (:l1, :l2, :max), k ∈ (1, 7, 20, 40)
        p₁ = PWA._full_order_plane_1d(X, fc.values, k, metric)
        p₂ = PWA._full_order_plane(X, fc.values, k, metric, optimizer)
        gaps(p) = [z - evaluate(p, x) for (x, z) ∈ fc]
        obj(g) = metric == :l1 ? sum(g) : metric == :l2 ? sum(abs2, g) : maximum(g)
        @test minimum(gaps(p₁)) ≥ -1e-6
        @test obj(gaps(p₁)) ≈ obj(gaps(p₂)) rtol = 1e-6 atol = 1e-8
    end

    # Enforcing curvature in one dimension gives the same total change as the
    # general formulation (here applied to the same points embedded in 2D)
    f2 = FunctionEvaluations([(xᵢ, 0.0) for xᵢ ∈ x], z)
    for metric ∈ (:l1, :max)
        v₁ = enforce_curvature(f, Convex(), optimizer, metric).values
        v₂ = enforce_curvature(f2, Convex(), optimizer, metric).values
        change(v) = metric == :l1 ? sum(abs, v .- z) : maximum(abs, v .- z)
        @test change(v₁) ≈ change(v₂) rtol = 1e-5
        s = diff(v₁) ./ diff(x)
        @test all(diff(s) .≥ -1e-6)
    end

    # Larger one-dimensional problems are solved quickly without an LP solver
    # for the full order planes and the reduced order planes
    xN = collect(range(-1, 1; length = 400))
    fN = FunctionEvaluations(tuple.(xN), xN .^ 2 .+ 0.3 .* xN)
    pwa = approx(fN, Convex(), FullOrder(optimizer = optimizer))
    @test length(pwa.planes) == 400
    alg = Progressive(optimizer = optimizer, tolerance = 0.005, metric = :max)
    pwa = approx(fN, Convex(), alg)
    @test PWA._tolerance_error(fN, pwa, :max) ≤ 0.005
end

@testset "Progressive fitting within a band" begin
    x = collect(range(-1, 1; length = 50))
    f = FunctionEvaluations(tuple.(x), x .^ 2 .+ 0.3 .* x)
    alg(δ, strict) = Progressive(; optimizer, tolerance = δ, metric = :max, strict)
    dev(f, pwa) = [evaluate(pwa, xᵢ) - zᵢ for (xᵢ, zᵢ) ∈ f]

    # Two-sided band: within the tolerance above and below, with fewer planes
    # than when all planes must lie below the data
    none = approx(f, Convex(), alg(0.01, :none))
    outer = approx(f, Convex(), alg(0.01, :outer))
    @test maximum(abs, dev(f, none)) ≤ 0.01 + 1e-6
    @test length(none.planes) ≤ length(outer.planes)
    @test all(dev(f, outer) .≤ 1e-6)

    # One-sided band above the data
    inner = approx(f, Convex(), alg(0.01, :inner))
    @test all(-1e-6 .≤ dev(f, inner) .≤ 0.01 + 1e-6)

    # Noisy data is handled without enforcing curvature, and the result is
    # within the tolerance of the original data
    xn = collect(range(0, 2; length = 100))
    fn = FunctionEvaluations(tuple.(xn), exp.(xn) .+ 0.02 .* randn(StableRNG(1), 100))
    pwa = approx(fn, Convex(), alg(0.05, :none))
    @test maximum(abs, dev(fn, pwa)) ≤ 0.05 + 1e-6
    @test length(pwa.planes) < 20

    # No convex function above the noisy data within the tolerance
    err = try
        approx(fn, Convex(), alg(0.05, :inner))
    catch e
        sprint(showerror, e)
    end
    @test occursin("tolerance", err)

    # Concave approximations: the band is mirrored
    fc = FunctionEvaluations(tuple.(x), -x .^ 2)
    pwa = approx(fc, Concave(), alg(0.02, :none))
    @test maximum(abs, dev(fc, pwa)) ≤ 0.02 + 1e-6
    pwa = approx(fc, Concave(), alg(0.02, :inner))
    @test all(-0.02 - 1e-6 .≤ dev(fc, pwa) .≤ 1e-6)

    # Two-dimensional data
    vals = PWA._sample_uniform(
        x -> x[1]^2 + x[2]^2 + 0.5x[1] * x[2],
        [(-1, 1), (-1, 1)],
        6,
    )
    pwa = approx(vals, Convex(), alg(0.05, :none))
    @test maximum(abs, dev(vals, pwa)) ≤ 0.05 + 1e-6

    # All planes are active at some data point
    for p ∈ (none, outer, inner)
        U = PWA._update_allocation(f, p)
        @test sort(unique(last.(U))) == 1:length(p.planes)
    end

    @test_throws ArgumentError approx(f, Convex(), alg(0.01, :foo))
end

@testset "Progressive robustness" begin
    # Tolerances that require many planes, where planes are also removed
    # along the way (previously an index error)
    x = collect(range(-1, 1; length = 30))
    f = FunctionEvaluations(tuple.(x), x .^ 2)
    for metric ∈ (:l1, :l2, :max), tol ∈ (0.0, 1e-12)
        alg = Progressive(optimizer = optimizer, tolerance = tol, metric = metric)
        pwa = approx(f, Convex(), alg)
        @test PWA._tolerance_error(f, pwa, metric) ≤ 1e-9
        @test length(pwa.planes) ≤ length(f)
    end

    # Two-dimensional data with a moderate tolerance
    vals = PWA._sample_uniform(x -> exp(x[1]) + x[2]^2, [(-1, 1), (-1, 1)], 10)
    alg = Progressive(optimizer = optimizer, tolerance = 0.05, metric = :max)
    pwa = approx(vals, Convex(), alg)
    @test PWA._tolerance_error(vals, pwa, :max) ≤ 0.05

    # A single data point
    one = FunctionEvaluations([(0.5,)], [0.25])
    pwa = approx(one, Convex(), FullOrder(optimizer = optimizer))
    @test all(isfinite, pwa.planes[1].α) && isfinite(pwa.planes[1].β)
    @test evaluate(pwa, (0.5,)) ≈ 0.25
    alg = Progressive(optimizer = optimizer, tolerance = 0.01, metric = :max)
    @test evaluate(approx(one, Convex(), alg), (0.5,)) ≈ 0.25
    alg = Progressive(optimizer = optimizer, tolerance = 0.01, strict = :none)
    @test abs(evaluate(approx(one, Convex(), alg), (0.5,)) - 0.25) ≤ 0.01
end

@testset "Progressive band independent of data scale" begin
    # The band is met regardless of offsets and scales of points and values
    x = collect(range(-1, 1; length = 100))
    for (xs, zs, zo) ∈ ((1.0, 1.0, 1e5), (1e4, 1.0, 1e5), (1.0, 1e-3, 0.0), (1.0, 1e4, -1e6))
        f = FunctionEvaluations(tuple.(xs .* x .+ 3xs), zs .* (x .^ 2 .+ 0.3 .* x) .+ zo)
        for strict ∈ (:none, :inner)
            alg = Progressive(optimizer = optimizer, tolerance = 0.01zs, strict = strict)
            pwa = approx(f, Convex(), alg)
            @test PWA._tolerance_error(f, pwa, :max) ≤ 0.01zs * (1 + 1e-5)
            @test length(pwa.planes) ≤ 25
        end
    end

    # Normalization of planes is undone exactly
    pts = FunctionEvaluations([(1.0, 10.0), (3.0, -2.0), (2.0, 4.0)], [5.0, 7.0, 6.0])
    X, μx, σx = PWA._normalize_points(pts)
    zn, z̄, s = PWA._normalize_values(pts.values, 0.1)
    q = Plane((0.3, -0.7), 0.2)
    qo = PWA._denormalize(q, μx, σx, z̄, s)
    for (k, pt) ∈ enumerate(pts.points)
        @test evaluate(qo, pt) ≈ s * evaluate(q, Tuple(X[:, k])) + z̄
    end
end
