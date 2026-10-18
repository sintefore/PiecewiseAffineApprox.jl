@testset "1D Convex Approx" begin
    x = collect(range(-1, 1; length = 10))
    z = x .^ 2

    pwa1 = approx(
        FunctionEvaluations(Tuple.(x), z),
        Convex(),
        MILP(optimizer = optimizer, metric = :l1, planes = 5),
    )
    @test length(pwa1.planes) == 5
    @test issorted((p.α for p ∈ pwa1.planes))
    @test isapprox(PWA.evaluate(pwa1, 0.4), 0.16, atol = 0.035)

    pwa2 = approx(
        x -> x[1]^2,
        [(-1, 1)],
        Convex(),
        MILP(optimizer = optimizer, metric = :l1, planes = 5),
    )
    @test length(pwa2.planes) == 5
    @test isapprox(PWA.evaluate(pwa2, 0.4), 0.1733, atol = 0.015)

    pwa3 = approx(
        x -> x[1]^2,
        [(-1, 1)],
        Convex(),
        MILP(optimizer = optimizer, metric = :l1, planes = 5, strict = :strict),
    )
    @test length(pwa2.planes) == 5
    @test isapprox(PWA.evaluate(pwa3, 0.4), 0.16, atol = 0.03)

    # pwa4 = convex_linearization(x -> x^2, -1, 1, optimizer; planes=5, method=:ipol)
    pwa4 = approx(
        x -> x[1]^2,
        [(-1, 1)],
        Convex(),
        Interpol(optimizer = optimizer, planes = 10),
    )

    @test isapprox(PWA.evaluate(pwa4, 0.4), 0.16, atol = 0.015)

    # Test correctness of strictness for :outer and :inner
    pwa_outer = approx(
        FunctionEvaluations(Tuple.(x), z),
        Convex(),
        MILP(optimizer = optimizer, metric = :l1, planes = 5, strict = :outer),
    )

    pwa_inner = approx(
        FunctionEvaluations(Tuple.(x), z),
        Convex(),
        MILP(optimizer = optimizer, metric = :l1, planes = 5, strict = :inner),
    )

    for i ∈ eachindex(x)
        # Convex, :outer should be below original points
        @test PWA.evaluate(pwa_outer, x[i]) ≤ z[i] + 1e-6
        # Convex, :inner should be above original points
        @test PWA.evaluate(pwa_inner, x[i]) ≥ z[i] - 1e-6
    end

    #Check approximation for flat function with releatively high values:
    I = 10
    xmat = 2 * rand(rng, 2, I) .- 1
    x = [Tuple(xmat[:, i]) for i ∈ 1:size(xmat, 2)]
    z = [10_000 + 0.001 * p[1]^2 + 0.001 * p[2]^2 for p ∈ x]
    vals = FunctionEvaluations(x, z)
    pwa = approx(
        vals,
        Convex(),
        MILP(optimizer = optimizer, planes = 4, metric = :l1),
    )
    @test PWA.evaluate(pwa, (0, 0)) ≈ 10_000.0 atol = 0.1
end

@testset "MILP robustness" begin
    x = collect(range(-1, 1; length = 10))

    @testset "Invariance to shift in function values" begin
        for strict ∈ (:none, :outer, :inner)
            err = map((0.0, 10.0)) do c
                pwa = approx(
                    FunctionEvaluations(Tuple.(x), c .- x .^ 2),
                    Concave(),
                    MILP(; optimizer, planes = 3, strict, metric = :l1),
                )
                return sum(abs(PWA.evaluate(pwa, xᵢ) - (c - xᵢ^2)) for xᵢ ∈ x)
            end
            @test err[1] ≈ err[2] atol = 1e-6
        end
    end

    @testset "Duplicate points" begin
        pwa = approx(
            FunctionEvaluations(Tuple.([0.0, 0.5, 0.5, 1.0]), [0.0, 0.25, 0.3, 1.0]),
            Convex(),
            MILP(; optimizer, planes = 2, metric = :l1),
        )
        @test length(pwa.planes) == 2
    end

    @testset "Small domain 2D" begin
        g = range(0, 0.1; length = 5)
        pts = vec([(a, b) for a ∈ g, b ∈ g])
        z = [a^2 + b^2 for (a, b) ∈ pts]
        pwa = approx(
            FunctionEvaluations(pts, z),
            Convex(),
            MILP(; optimizer, planes = 4, metric = :max),
        )
        err = maximum(abs(PWA.evaluate(pwa, p) - zᵢ) for (p, zᵢ) ∈ zip(pts, z))
        @test err ≈ 0.03125 * (maximum(z) - minimum(z)) rtol = 1e-4
    end
end
