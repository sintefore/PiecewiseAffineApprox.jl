
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
    # Disabled as it gives inconsistent results locally and in CI
    # @test length(pwa_red.planes) == 10
    @test evaluate(pwa_red, (0, 0)) ≈ 0.024 atol = 0.001

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
