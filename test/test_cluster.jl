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
