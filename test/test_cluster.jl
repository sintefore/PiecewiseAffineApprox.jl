@testset "Cluster reproducibility" begin
    pts = [(a, b) for a ∈ range(-1, 1, length = 6) for b ∈ range(-1, 1, length = 6)]
    vals = [a^2 + b^2 + a * b for (a, b) ∈ pts]
    fe = FunctionEvaluations(pts, vals)

    alg(seed) = Cluster(; optimizer, planes = 4, trials = 5, rng = StableRNG(seed))
    pwa1 = approx(fe, Convex(), alg(42))
    pwa2 = approx(fe, Convex(), alg(42))
    @test pwa1.planes == pwa2.planes

    # Concave approximations go through the same code path
    cv1 = approx(FunctionEvaluations(pts, -vals), Concave(), alg(7))
    cv2 = approx(FunctionEvaluations(pts, -vals), Concave(), alg(7))
    @test cv1.planes == cv2.planes

    # Reusing the algorithm object advances its generator, but stays reproducible
    a1 = alg(1)
    a2 = alg(1)
    first1, second1 = approx(fe, Convex(), a1), approx(fe, Convex(), a1)
    first2, second2 = approx(fe, Convex(), a2), approx(fe, Convex(), a2)
    @test first1.planes == first2.planes
    @test second1.planes == second2.planes

    # Default constructor still works without an explicit generator
    @test Cluster(; optimizer).rng === Random.default_rng()
end

@testset "Cluster random sampling" begin
    # Samples match the mean and covariance of the data
    X = [1.0 2.0 4.0 7.0 3.0; 0.5 -1.0 2.0 1.0 0.0]
    P = PWA._sample_normal(StableRNG(1), X, 100_000)
    @test vec(PWA.mean(P, dims = 2)) ≈ vec(PWA.mean(X, dims = 2)) rtol = 0.02
    @test PWA.cov(P, dims = 2) ≈ PWA.cov(X, dims = 2) rtol = 0.02

    # Collinear points give a singular covariance
    pts = [(t, 2t) for t ∈ range(0, 1, length = 10)]
    fe = FunctionEvaluations(pts, [a^2 + b^2 for (a, b) ∈ pts])
    pwa = approx(fe, Convex(), Cluster(; optimizer, planes = 3, rng = StableRNG(3)))
    @test maximum(abs(evaluate(pwa, p) - v) for (p, v) ∈ fe) < 0.25
end

@testset "Cluster normalization" begin
    # Normalized data has zero mean and identity covariance
    X = [1.0 2.0 4.0 7.0 3.0; 0.5 -1.0 2.0 1.0 0.0]
    μ, T = PWA._normalization(X)
    X̃ = T * (X .- μ)
    @test vec(PWA.mean(X̃, dims = 2)) ≈ zeros(2) atol = 1e-12
    @test PWA.cov(X̃, dims = 2) ≈ I

    # Back transformation gives the same function values
    pwa = PWAFunc{Convex,2}([Plane([1.0, -2.0], 0.5), Plane([-0.5, 3.0], 1.0)])
    pwaₒ = PWA._denormalize(pwa, μ, T)
    @test all(evaluate(pwaₒ, x) ≈ evaluate(pwa, x̃) for (x, x̃) ∈ zip(eachcol(X), eachcol(X̃)))

    # Result is (close to) invariant to scaling and shifting the input
    pts = [(a, b) for a ∈ range(-2, 2, length = 8) for b ∈ range(-2, 2, length = 8)]
    vals = [log(exp(a) + exp(b) + 1) for (a, b) ∈ pts]
    alg = Cluster(; optimizer, planes = 5, trials = 5, metric = :l2, rng = StableRNG(1))
    rms(pts, pwa) = sqrt(PWA.mean((evaluate(pwa, p) - v)^2 for (p, v) ∈ zip(pts, vals)))

    pwa = approx(FunctionEvaluations(pts, vals), Convex(), alg)
    scaled = [(1000a + 50, b) for (a, b) ∈ pts]
    alg = Cluster(; optimizer, planes = 5, trials = 5, metric = :l2, rng = StableRNG(1))
    pwaₛ = approx(FunctionEvaluations(scaled, vals), Convex(), alg)
    @test length(pwaₛ.planes) == 5
    @test rms(scaled, pwaₛ) ≤ 1.25 * rms(pts, pwa)
end

@testset "Cluster least squares" begin
    X = [1.0 2.0 4.0 7.0 3.0 -1.0; 0.5 -1.0 2.0 1.0 0.0 0.3]
    z = [2.0, 1.0, 5.0, 9.0, 2.5, 0.7]

    # Closed form agrees with solving the least squares problem with a QP
    m = Model(optimizer)
    @variable(m, a[1:2])
    @variable(m, b)
    @objective(m, Min, sum((z[i] - a' * X[:, i] - b)^2 for i ∈ 1:6))
    optimize!(m)
    â, b̂ = PWA._least_squares_fit(X, z)
    @test â ≈ value.(a) atol = 1e-6
    @test b̂ ≈ value(b) atol = 1e-6

    # Underdetermined subsets give zero slope where the data has no variation
    â, b̂ = PWA._least_squares_fit(X[:, 1:1], z[1:1])
    @test â == zeros(2)
    @test b̂ == z[1]
    â, b̂ = PWA._least_squares_fit([1.0 3.0; 2.0 2.0], [1.0, 5.0])
    @test â ≈ [2.0, 0.0] atol = 1e-6
    @test b̂ ≈ -1.0 atol = 1e-6

    # No optimizer needed for :l2, but required for the other metrics
    pts = [(a, b) for a ∈ range(-1, 1, length = 6) for b ∈ range(-1, 1, length = 6)]
    fe = FunctionEvaluations(pts, [a^2 + b^2 for (a, b) ∈ pts])
    pwa = approx(fe, Convex(), Cluster(; metric = :l2, planes = 4, rng = StableRNG(1)))
    @test length(pwa.planes) == 4
    @test_throws ErrorException approx(fe, Convex(), Cluster(; metric = :l1))
    @test_throws ErrorException approx(
        fe,
        Convex(),
        Cluster(; metric = :l2, strict = :outer),
    )
end
