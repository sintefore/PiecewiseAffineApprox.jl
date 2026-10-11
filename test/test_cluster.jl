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

@testset "Cluster refill of subsets" begin
    # Subsets must spread out in all directions spanned by the data
    X = randn(StableRNG(1), 3, 200)
    μ, T = PWA._normalization(X)
    X̃ = T * (X .- μ)
    determined = PWA._determined_test(X̃, 1e-4)
    @test determined(collect(1:50))
    @test !determined([1, 2, 3])                    # too few points
    Y = X̃[:, 1:10]
    Y[3, :] .= 0.5                                  # coplanar points
    @test !PWA._determined_test(hcat(X̃, Y), 1e-4)(collect(201:210))

    # For data on a line, subsets only need to spread along the line
    L = [1.0, 2.0, -1.0] * range(0, 1, length = 20)'
    μ, T = PWA._normalization(L)
    @test PWA._determined_test(T * (L .- μ), 1e-4)([1, 5, 9])

    # Refilling splits the subsets with the largest weight
    𝒫 = [collect(1:150), collect(151:200)]
    w = [i ≤ 150 ? 0.0 : 1.0 for i ∈ 1:200]
    PWA._refill!(𝒫, X̃, w, 4, determined)
    @test length(𝒫) == 4
    @test all(determined, 𝒫)
    @test sort(reduce(vcat, 𝒫)) == 1:200
    @test 𝒫[1] == collect(1:150)                    # zero weight, not split

    # The requested number of planes is returned
    pts = [(a, b) for a ∈ range(-2, 2, length = 25) for b ∈ range(-2, 2, length = 25)]
    fe = FunctionEvaluations(pts, [a^2 + exp(b) + a * b for (a, b) ∈ pts])
    for metric ∈ (:l2, :l1)
        alg = Cluster(; optimizer, metric, planes = 8, trials = 3, rng = StableRNG(2))
        @test length(approx(fe, Convex(), alg).planes) == 8
    end

    # No trial ends in a poor local optimum caused by a dominating plane. Before,
    # with the example from Magnani & Boyd, some trials gave rms ≈ 3.
    V = -5:5
    X = Float64[c[i] for i ∈ 1:3, c ∈ vec(collect(Iterators.product(V, V, V)))]
    z = [log(sum(exp, c)) for c ∈ eachcol(X)]
    μ, T = PWA._normalization(X)
    X̃ = T * (X .- μ)
    errors = map(1:20) do i
        𝒫 = PWA._random_partition(StableRNG(i), X̃, 12)
        pwa = PWA._refine_partition(X̃, z, 𝒫, 50, :l2, nothing, :none)
        @test length(pwa.planes) == 12
        return PWA._approx_error(X̃, z, pwa, :rms)
    end
    @test maximum(errors) < 0.1
end

@testset "Cluster best fit and defaults" begin
    @test Cluster().metric == :l2

    # Each trial returns the best fit over its iterations, so allowing more
    # iterations never gives a worse result
    V = -5:5
    X = Float64[c[i] for i ∈ 1:3, c ∈ vec(collect(Iterators.product(V, V, V)))]
    z = [log(sum(exp, c)) for c ∈ eachcol(X)]
    μ, T = PWA._normalization(X)
    X̃ = T * (X .- μ)
    for seed ∈ 1:5
        𝒫 = PWA._random_partition(StableRNG(seed), X̃, 12)
        errors = map(0:15) do itlim
            pwa = PWA._refine_partition(X̃, z, 𝒫, itlim, :l2, nothing, :none)
            return PWA._approx_error(X̃, z, pwa, :l2)
        end
        @test issorted(errors, rev = true)
    end
end

@testset "Cluster validation and edge cases" begin
    pts = [(a, b) for a ∈ range(-1, 1, length = 6) for b ∈ range(-1, 1, length = 6)]
    fe = FunctionEvaluations(pts, [a^2 + b^2 for (a, b) ∈ pts])
    @test_throws ArgumentError approx(fe, Convex(), Cluster(; metric = :foo))
    @test_throws ArgumentError approx(fe, Convex(), Cluster(; optimizer, metric = :foo))
    @test_throws ArgumentError approx(fe, Convex(), Cluster(; optimizer, strict = :bar))
    @test_throws ArgumentError approx(fe, Convex(), Cluster(; planes = 0))
    @test_throws ArgumentError PWA._approx_error(zeros(2, 1), [0.0], PWAFunc{Convex,2}(), :foo)

    # Errors in the trials are reported directly, not as failed tasks
    failing = () -> error("no solver available")
    err = try
        approx(fe, Convex(), Cluster(; optimizer = failing, metric = :l1, trials = 2))
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("no solver available", sprint(showerror, err))

    # A single point, and identical points
    one = FunctionEvaluations([(0.5, 0.5)], [1.0])
    @test evaluate(approx(one, Convex(), Cluster(; planes = 2, trials = 2)), (0.5, 0.5)) ≈ 1.0
    same = FunctionEvaluations(fill((1.0, 2.0), 5), ones(5))
    @test evaluate(approx(same, Convex(), Cluster(; planes = 2, trials = 2)), (1.0, 2.0)) ≈ 1.0

    # A coordinate with a small but genuine range is not treated as degenerate:
    # the fit is as good as with the coordinate rescaled
    pts3 = [(a, b) for a ∈ range(-1, 1, length = 10) for b ∈ range(0, 1e-6, length = 10)]
    z3 = [a^2 + (b / 1e-6)^2 for (a, b) ∈ pts3]
    rms(f, pwa) = sqrt(PWA.mean((evaluate(pwa, x) - v)^2 for (x, v) ∈ f))
    f3 = FunctionEvaluations(pts3, z3)
    f3r = FunctionEvaluations([(a, 1e6 * b) for (a, b) ∈ pts3], z3)
    alg() = Cluster(; planes = 6, trials = 5, rng = StableRNG(2))
    e3 = rms(f3, approx(f3, Convex(), alg()))
    e3r = rms(f3r, approx(f3r, Convex(), alg()))
    @test e3 ≤ 1.1 * e3r
    @test e3 < 0.1
end
