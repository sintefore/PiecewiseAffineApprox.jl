
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
        for (x, z) in PWA.point_vals(f)
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
    @test evaluate(pwa_red, (0, 0)) ≈ 0.0082 atol = 0.0001

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
    end
end
