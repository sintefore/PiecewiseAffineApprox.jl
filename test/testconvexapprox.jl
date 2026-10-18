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
        MILP(optimizer = optimizer, metric = :l1, planes = 5, strict = :outer),
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

    @testset "Strict and warm start" begin
        z = x .^ 2
        fe = FunctionEvaluations(Tuple.(x), z)
        for metric ∈ (:l1, :max), strict ∈ (:none, :outer, :inner)
            obj = map((false, true)) do warmstart
                pwa = approx(
                    fe,
                    Convex(),
                    MILP(; optimizer, planes = 3, metric, strict, warmstart),
                )
                e = [PWA.evaluate(pwa, xᵢ) - zᵢ for (xᵢ, zᵢ) ∈ zip(x, z)]
                strict == :outer && @test maximum(e) ≤ 1e-6
                strict == :inner && @test minimum(e) ≥ -1e-6
                return metric == :l1 ? sum(abs, e) : maximum(abs, e)
            end
            @test obj[1] ≈ obj[2] atol = 1e-6
        end
    end

    @testset "Inner as good as outer for symmetric data" begin
        # The approximation, not each plane, should lie above the data
        e = map((:outer, :inner)) do strict
            pwa = approx(
                FunctionEvaluations(Tuple.(x), x .^ 2),
                Convex(),
                MILP(; optimizer, planes = 3, metric = :max, strict),
            )
            return maximum(abs(PWA.evaluate(pwa, xᵢ) - xᵢ^2) for xᵢ ∈ x)
        end
        @test e[1] ≈ e[2] atol = 1e-6
    end

    @testset "Three dimensions" begin
        g = range(-1, 1; length = 3)
        pts = vec([(a, b, c) for a ∈ g, b ∈ g, c ∈ g])
        z = [sum(abs2, p) for p ∈ pts]
        fe = FunctionEvaluations(pts, z)
        pwa = approx(fe, Convex(), MILP(; optimizer, planes = 2))
        @test length(pwa.planes) == 2
        @test_throws ArgumentError approx(
            fe,
            Convex(),
            MILP(; optimizer, planes = 2, warmstart = false),
        )
    end

    @testset "Validation" begin
        fe = FunctionEvaluations(Tuple.(x), x .^ 2)
        @test_throws ArgumentError approx(
            fe,
            Convex(),
            MILP(; optimizer, strict = :strict),
        )
        @test_throws ArgumentError approx(
            fe,
            Convex(),
            MILP(; optimizer, metric = :rms),
        )
        @test_throws ArgumentError approx(fe, Convex(), MILP(; optimizer, planes = 0))
    end

    @testset "Invariance to scaling of the data" begin
        maxerr(pwa, fe) = maximum(
            abs(PWA.evaluate(pwa, p) - zᵢ) for (p, zᵢ) ∈ zip(fe.points, fe.values)
        )
        g = range(-1, 1; length = 4)
        pts = vec([(a, b) for a ∈ g, b ∈ g])
        for (points, f) ∈ ((Tuple.(x), p -> p[1]^2), (pts, p -> p[1]^2 + p[2]^2))
            err = map(((1.0, 0.0, 0.0), (1e-3, 5.0, 0.0), (1e3, 0.0, 1e6))) do (s, c, d)
                fe = FunctionEvaluations(
                    [c .+ s .* p for p ∈ points],
                    [d + s^2 * f(p) for p ∈ points],
                )
                pwa = approx(fe, Convex(), MILP(; optimizer, planes = 3, metric = :max))
                return maxerr(pwa, fe) / s^2
            end
            @test err[1] ≈ err[2] rtol = 1e-4
            @test err[1] ≈ err[3] rtol = 1e-4
        end
    end

    @testset "Recovery from too small big-M" begin
        fe = FunctionEvaluations(Tuple.(x), (x .^ 2 .+ 1) ./ 2)
        options = MILP(; optimizer, planes = 4, metric = :max, warmstart = false)
        pwaᵒᵖᵗ = approx(fe, Convex(), options)
        errᵒᵖᵗ = maximum(abs(PWA.evaluate(pwaᵒᵖᵗ, xᵢ) - zᵢ) for (xᵢ, zᵢ) ∈ fe)
        a, b = PWA._milp(fe, options; Mᵖ = fill(1e-3, length(x)))
        pwa = PWAFunc{Convex,1}([Plane((a[1, k],), b[k]) for k ∈ eachindex(b)])
        err = maximum(abs(PWA.evaluate(pwa, xᵢ) - zᵢ) for (xᵢ, zᵢ) ∈ fe)
        @test err ≈ errᵒᵖᵗ atol = 1e-6
        # Warns if the big-M may still be too small after the last doubling
        @test_logs (:warn, r"big-M") match_mode = :any PWA._milp(
            fe,
            options;
            Mᵖ = fill(1e-3, length(x)),
            rounds = 0,
        )
    end

    @testset "Symmetry breaking and warm start" begin
        function objective(a, b, fe, metric)
            e = [
                maximum(sum(a[:, k] .* p) + b[k] for k ∈ eachindex(b)) - zᵢ for
                (p, zᵢ) ∈ fe
            ]
            return metric == :l1 ? sum(abs, e) : maximum(abs, e)
        end
        rng = StableRNG(3)
        g = range(0, 1; length = 4)
        pts = vec([(a, b) for a ∈ g, b ∈ g])
        noisy = [abs(xᵢ) + 0.05 * randn(rng) for xᵢ ∈ x]
        cases = (
            FunctionEvaluations(Tuple.(x), noisy),
            FunctionEvaluations(pts, [(a - 0.4)^2 + (b - 0.5)^2 for (a, b) ∈ pts]),
        )
        for fe ∈ cases, metric ∈ (:l1, :max)
            options = MILP(; optimizer, planes = 3, metric)
            start = PWA._milp_start(fe, options)
            if length(fe.points[1]) == 1
                # The symmetry breaking constraints do not cut off optimal
                # solutions, and are satisfied by the start
                obj = map((false, true)) do symmetry
                    a, b = PWA._milp(fe, options; symmetry)
                    return objective(a, b, fe, metric)
                end
                @test obj[1] ≈ obj[2] rtol = 1e-3
                order = sortperm(eachindex(fe.points); by = p -> fe.points[p][1])
                @test issorted(start.active[order])
            end

            # The start is the best of the heuristic fits
            objᴴ = [
                PWA._milp_start(fe, options, planes).obj for
                planes ∈ PWA._heuristic_fits(fe, options)
            ]
            e = start.ẑ .- fe.values
            @test (metric == :l1 ? sum(abs, e) : maximum(abs, e)) ≈ minimum(objᴴ)
        end
    end
end

@testset "Sampling of functions" begin
    milp(K) = MILP(; optimizer, planes = K)
    # MILP: about 3 D K points in total, unchanged in 1D
    @test PWA._default_samples(milp(5), 1) == 15
    @test PWA._default_samples(milp(6), 2) == 6
    @test PWA._default_samples(milp(4), 3) == 4
    @test PWA._default_samples(Cluster(; planes = 4), 2) == 12
    # Algorithms without a number of planes: about 100 points in total
    progressive = Progressive(; optimizer, tolerance = 0.01)
    @test PWA._default_samples(progressive, 1) == 100
    @test PWA._default_samples(progressive, 2) == 10

    f(x) = x[1]^2 + x[2]^2
    bbox = [(-1, 1), (-1, 1)]
    pwa = approx(f, bbox, Convex(), MILP(; optimizer, planes = 3, metric = :max))
    @test length(pwa.planes) == 3
    pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), progressive)
    @test PWA.evaluate(pwa, 0.5) ≈ 0.25 atol = 0.02

    # The number of samples can be given; with 2 samples in each dimension the
    # corners are interpolated
    pwa = approx(f, bbox, Convex(), Cluster(; planes = 1); nsamples = 2)
    @test PWA.evaluate(pwa, (1.0, 1.0)) ≈ 2.0 atol = 1e-8
    @test_throws ArgumentError approx(f, bbox, Convex(), milp(3); nsamples = 1)
end
