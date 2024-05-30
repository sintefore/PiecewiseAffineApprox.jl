@testset "Big-M calculations" begin
    function mat2tuples(x)
        return collect(Tuple(x'[:, i]) for i ∈ 1:size(x', 2))
    end

    # big-M for a 1D function
    x = collect(range(-1, 1; length = 10))
    z = x .^ 2

    M = PWA._linear_big_M(FunctionEvaluations(tuple.(x), z))
    @test M ≈ 3.91 atol = 0.01

    # big-M for a 2D function
    xg = collect(range(-1, 1; length = 4))
    yg = collect(range(-1, 1; length = 4))
    X = [repeat(xg, inner = [size(yg, 1)]) repeat(yg, outer = [size(xg, 1)])]
    z = X[:, 1] .^ 2 + X[:, 2] .^ 2
    np = 5

    M = PWA._linear_big_M(FunctionEvaluations(mat2tuples(X), z))
    @test M ≈ 14.67 atol = 0.01

    @testset "Duplicate inputs" begin
        # 1D
        x = [302.3633312955, 300.6680553405, 274.281686739, 274.281686739]
        z = [-184.9368247, -183.9371662, -167.94263, -167.94263]

        M¹ = PWA._linear_big_M(FunctionEvaluations(tuple.(x), z))
        @test M¹ ≈ 0.479 atol = 0.001

        M² =
            PWA._linear_big_M(FunctionEvaluations(tuple.(unique(x)), unique(z)))
        @test M¹ == M²

        # Duplicate with different function value
        z[4] = -166.3
        M³ = PWA._linear_big_M(FunctionEvaluations(tuple.(x), z))
        @test M³ ≈ 2.286 atol = 0.001

        # 2D
        xg = collect(range(-1, 1; length = 4))
        yg = collect(range(-1, 1; length = 4))

        X =
            [repeat(xg, inner = [size(yg, 1)]) repeat(yg, outer = [size(xg, 1)])]
        z = X[:, 1] .^ 2 + X[:, 2] .^ 2

        # Duplicate point with different function value
        X = vcat(X, [-1 1 / 3])
        z = vcat(z, 1.5)

        M = PWA._linear_big_M(FunctionEvaluations(mat2tuples(X), z))
        @test M ≈ 15.95 atol = 0.01
    end

    @testset "Close points" begin
        # 1D
        x = collect(range(-1, 1; length = 10))
        z = x .^ 2

        # Not close enough
        x1 = vcat(x, 0.55)
        z1 = vcat(z, 0.41)
        M¹ = PWA._linear_big_M(FunctionEvaluations(tuple.(x1), z1))
        @test M¹ ≈ 9.68 atol = 0.001

        # Close points
        x2 = vcat(x, 0.5555)
        z2 = vcat(z, 0.41)
        M² = PWA._linear_big_M(FunctionEvaluations(tuple.(x2), z2))
        @test M² ≈ 3.911 atol = 0.001
    end
end
