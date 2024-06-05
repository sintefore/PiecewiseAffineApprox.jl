
@testset "Progressive fitting" begin

    # 1D
    x = collect(range(-1, 1; length = 30))
    z = x .^ 2
    f = FunctionEvaluations(tuple.(x), z)
    for (tol, pen) in [(0.05, :max), (2.5, :l1), (0.5, :l2)]
        pwa_red = approx(
            f,
            Convex(),
            ProgressiveFitting(
                optimizer = optimizer,
                tolerance = tol,
                pen = pen,
            ),
        )
        for (x, z) in f
            pv = evaluate(pwa_red, x)
            @test z ≥ pv || z ≈ pv
        end
    end

    # 2D
    g(x) = x[1]^2 + x[2]^2
    vals = PiecewiseAffineApprox._sample_uniform(g, [(-1, 1), (-1, 1)], 10)
    pwa_red = approx(
        vals,
        Convex(),
        ProgressiveFitting(optimizer = optimizer, tolerance = 0.2, pen = :max),
    )
    @test length(pwa_red.planes) == 9
    @test evaluate(pwa_red, (0, 0)) ≈ -0.024 atol = 0.001

    @testset "Nonconvex" begin
        h(x) = sin(5 * x[1]) * (x[1]^2 + x[2]^2)
        vals = PiecewiseAffineApprox._sample_uniform(h, [(-1, 1), (-1, 1)], 10)

        @test_throws ErrorException approx(
            vals,
            Convex(),
            ProgressiveFitting(
                optimizer = optimizer,
                tolerance = 2.0,
                pen = :l2,
            ),
        )

        vals_c = PWA.convexify(vals, optimizer, :l1)
        @test length(vals_c) == length(vals)

        pwa_con = approx(
            vals_c,
            Convex(),
            ProgressiveFitting(
                optimizer = optimizer,
                tolerance = 0.1,
                pen = :l2,
            ),
        )

        x = collect(range(-1, 1; length = 30))
        z = x .^ 2
        z[5] += 0.1
        f = FunctionEvaluations(tuple.(x), z)

        fc = PWA.convexify(f, optimizer)
        @test length(fc) == length(f)
        @test fc.values[5] ≈ 0.5291 atol = 0.001

        pwa_con = approx(
            fc,
            Convex(),
            ProgressiveFitting(
                optimizer = optimizer,
                tolerance = 0.2,
                pen = :l1,
            ),
        )
        @test PWA._planes(pwa_con) == 8
    end
end
