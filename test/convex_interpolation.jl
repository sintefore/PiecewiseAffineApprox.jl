
@testset "Convex Interpolation" begin
    # Points to approximate
    x = [i for i ∈ -2:0.1:2]
    z = x .^ 2

    # Test interpolation routine
    pwa = PWA._interpolatepw(x, z, Interpol(optimizer = optimizer, planes = 5))

    @test length(pwa.x) == 6
    @test length(pwa.z) == 6

    @test pwa.z[2] ≈ 1.44
    @test pwa.z[3] ≈ 0.16

    # Test convexification
    cpwa = PWA._convexify1D(pwa, Interpol(optimizer = optimizer))

    @test cpwa.x[2] ≈ -1.2
    @test cpwa.x[3] ≈ -0.4
end

@testset "Interpolation input handling" begin
    x = collect(-2:0.25:2)
    z = x .^ 2 .+ 0.3 .* x
    alg = Interpol(; optimizer, planes = 4)
    ref = approx(FunctionEvaluations(Tuple.(x), z), Convex(), alg)
    same(pwa) = all(evaluate(pwa, (t,)) ≈ evaluate(ref, (t,)) for t ∈ x)

    # Unsorted input gives the same result as sorted input
    p = [9, 3, 16, 1, 12, 5, 14, 7, 2, 11, 17, 4, 8, 15, 6, 10, 13]
    @test same(approx(FunctionEvaluations(Tuple.(x[p]), z[p]), Convex(), alg))

    # Duplicate points are merged
    xd = vcat(x, x[[3, 8]])
    zd = vcat(z, z[[3, 8]])
    @test same(approx(FunctionEvaluations(Tuple.(xd), zd), Convex(), alg))

    # Duplicates with different values are replaced by their average
    xs, zs = PWA._sort_unique([2.0, 1.0, 2.0, 0.0, 2.0], [1.0, 5.0, 2.0, 7.0, 6.0])
    @test xs == [0.0, 1.0, 2.0]
    @test zs ≈ [7.0, 5.0, 3.0]

    # Concave approximations use the same path
    cv = approx(FunctionEvaluations(Tuple.(x[p]), -z[p]), Concave(), alg)
    @test all(evaluate(cv, (t,)) ≈ -evaluate(ref, (t,)) for t ∈ x)

    # Number of planes must be supported by the number of distinct points
    few = FunctionEvaluations([(0.0,), (1.0,), (1.0,), (2.0,)], [0.0, 1.0, 1.0, 4.0])
    @test length(approx(few, Convex(), Interpol(; optimizer, planes = 2)).planes) == 2
    @test_throws ArgumentError approx(few, Convex(), Interpol(; optimizer, planes = 3))
    @test_throws ArgumentError approx(few, Convex(), Interpol(; optimizer, planes = 0))

    # Only 1D is supported
    pts = [(a, b) for a ∈ 0:2 for b ∈ 0:2]
    fe2 = FunctionEvaluations(pts, [a^2 + b^2 for (a, b) ∈ pts])
    @test_throws ArgumentError approx(fe2, Convex(), alg)
end

@testset "Interpolation dynamic programming" begin
    rng = StableRNG(7)

    # Squared deviations by prefix sums equal direct computation, also with
    # large offsets in the data
    x = sort(rand(rng, 30)) .+ 1e3
    z = (x .- 1e3) .^ 2 .+ 0.1 .* randn(rng, 30) .+ 1e4
    p = PWA._segment_penalties(x, z, :l2)
    for i ∈ 1:29, j ∈ i+1:30
        c = (z[j] - z[i]) / (x[j] - x[i])
        direct = sum((z[i] + c * (x[k] - x[i]) - z[k])^2 for k ∈ i:j)
        @test p[i, j] ≈ direct atol = 1e-8
    end

    # The dynamic program finds the optimal knots for any number of segments,
    # compared with enumerating all choices of interior knots. For :max the
    # largest segment penalty is minimized, with ties broken by the sum.
    N = 9
    x = sort(rand(rng, N))
    z = sin.(4 .* x) .+ 0.2 .* randn(rng, N)
    for metric ∈ (:l1, :l2, :max)
        p = PWA._segment_penalties(x, z, metric)
        K = 5
        minimax = metric == :max
        path, costs = PWA._shortest_path_segments(p, K; minimax)
        edges(q) = [p[q[k], q[k+1]] for k ∈ 1:length(q)-1]
        pathcost(q) = minimax ? maximum(edges(q)) : sum(edges(q))
        paths(k) = ([1; collect(c); N] for c ∈ PWA.combinations(2:N-1, k - 1))
        for k ∈ 1:K
            @test costs[k] ≈ minimum(pathcost, paths(k))
        end
        @test length(path) == K + 1
        @test path[1] == 1 && path[end] == N && issorted(path)
        @test pathcost(path) ≈ costs[K]
        if minimax
            optimal = [q for q ∈ paths(K) if pathcost(q) == costs[K]]
            @test sum(edges(path)) ≈ minimum(q -> sum(edges(q)), optimal)
        end
    end

    # Minimax gives a smaller maximum deviation than minimizing the sum of
    # the segment maxima
    x = collect(range(0, 3, length = 40))
    z = exp.(x)
    p = PWA._segment_penalties(x, z, :max)
    mm, _ = PWA._shortest_path_segments(p, 4; minimax = true)
    sm, _ = PWA._shortest_path_segments(p, 4; minimax = false)
    maxdev(q) = maximum(p[q[k], q[k+1]] for k ∈ 1:4)
    @test maxdev(mm) < maxdev(sm)

    # No optimization solver is used for selecting the knots
    pwa = PWA._interpolatepw(x, z, Interpol(; optimizer = nothing, planes = 3))
    @test length(pwa.x) == 4
end
