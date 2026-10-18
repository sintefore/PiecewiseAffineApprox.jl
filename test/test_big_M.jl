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

    @testset "Scale invariance 2D" begin
        # big-M should scale with the data, not drop to zero for small domains
        g = range(-1, 1; length = 4)
        pts = vec([(a, b) for a ∈ g, b ∈ g])
        z = [a^2 + b^2 for (a, b) ∈ pts]
        M = PWA._linear_big_M(FunctionEvaluations(pts, z))
        for s ∈ (10.0, 0.1, 0.01)
            Mˢ = PWA._linear_big_M(
                FunctionEvaluations([s .* p for p ∈ pts], s^2 .* z),
            )
            @test Mˢ ≈ s^2 * M
        end
    end

    @testset "Same as brute force" begin
        # Direct implementation of the definitions: all segments (1D) or
        # triangles (2D) without data points in the interior
        function brute_force_1d(x, z)
            δ = (maximum(x) - minimum(x)) / (100 * length(x))
            M = zeros(length(x))
            for i ∈ eachindex(x), j ∈ eachindex(x)
                x[j] - x[i] > δ || continue
                any(p -> x[i] + δ < p < x[j] - δ, x) && continue
                a = (z[j] - z[i]) / (x[j] - x[i])
                M .= max.(M, z .- a .* x .- (z[i] - a * x[i]))
            end
            return 1.1 .* M
        end
        function brute_force_2d(pts, z)
            lo, hi = minimum.(zip(pts...)), maximum.(zip(pts...))
            x = [(p .- lo) ./ (hi .- lo) for p ∈ pts]
            M = zeros(length(x))
            for i ∈ eachindex(x), j ∈ i+1:length(x), k ∈ j+1:length(x)
                v = (x[i], x[j], x[k])
                any(p -> p ∉ v && PWA._point_in_triangle(p, v...), x) && continue
                PWA._sliver_parameter(v...) > 0.025 || continue
                t = [(x[l]..., z[l]) for l ∈ (i, j, k)]
                n = cross(collect(t[1] .- t[2]), collect(t[1] .- t[3]))
                M .= max.(M, z .- [PWA._plane_f(p, n, dot(n, t[1])) for p ∈ x])
            end
            return 1.1 .* M
        end

        rng = StableRNG(1)
        for n ∈ (2, 5, 20)
            x = rand(rng, n)
            x = vcat(x, x[1], x[2] + 1e-6)
            z = x .^ 2 .+ 0.1 .* randn(rng, length(x))
            M = PWA._linear_big_M_points(FunctionEvaluations(tuple.(x), z))
            @test M ≈ brute_force_1d(x, z)

            pts = [(rand(rng), rand(rng)) for _ ∈ 1:n+1]
            pts = vcat(pts, pts[1:1], [(0.5, 0.5), (0.5, 0.25), (0.5, 0.75)])
            z = [a^2 + b^2 + 0.1 * randn(rng) for (a, b) ∈ pts]
            M = PWA._linear_big_M_points(FunctionEvaluations(pts, z))
            @test M ≈ brute_force_2d(pts, z)
        end
    end
end
