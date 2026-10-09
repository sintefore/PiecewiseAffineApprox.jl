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
