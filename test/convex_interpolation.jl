
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

@testset "Interpolation convexification" begin
    rng = StableRNG(11)

    # The convex dynamic program finds the best convex interpolant through a
    # subset of the points, compared with enumerating all subsets
    N = 9
    x = sort(rand(rng, N))
    z = sin.(5 .* x) .+ 0.3 .* randn(rng, N)
    σ = [i < j ? (z[j] - z[i]) / (x[j] - x[i]) : 0.0 for i ∈ 1:N, j ∈ 1:N]
    K = 4
    for metric ∈ (:l1, :l2, :max)
        p = PWA._segment_penalties(x, z, metric)
        minimax = metric == :max
        edges(q) = [p[q[k], q[k+1]] for k ∈ 1:length(q)-1]
        cost(q) = minimax ? maximum(edges(q)) : sum(edges(q))
        convex(q) = all(
            PWA._leq(σ[q[k], q[k+1]], σ[q[k+1], q[k+2]]) for k ∈ 1:length(q)-2
        )
        feasible = [
            [1; collect(c); N] for k ∈ 1:K for
            c ∈ PWA.combinations(2:N-1, k - 1) if convex([1; collect(c); N])
        ]
        path = PWA._convex_shortest_path(σ, p, K; minimax)
        @test convex(path)
        @test length(path) ≤ K + 1
        @test cost(path) ≈ minimum(cost, feasible)
        if minimax
            optimal = [q for q ∈ feasible if cost(q) == cost(path)]
            @test sum(edges(path)) ≈ minimum(q -> sum(edges(q)), optimal)
        end
    end

    # Data made convex has nondecreasing slopes, and convex data is unchanged
    x = collect(range(0, 3, length = 30))
    z = sin.(3 .* x)
    for metric ∈ (:l1, :l2, :max)
        v = PWA._convexify_values(x, z, metric, optimizer)
        @test PWA._isconvex(x, v)
        @test !PWA._isconvex(x, z)
        w = PWA._convexify_values(x, x .^ 2, metric, optimizer)
        @test maximum(abs.(w .- x .^ 2)) < 1e-5
    end

    fe = FunctionEvaluations(Tuple.(x), z)
    error(pwa, metric) = PWA._ipol_error(pwa, x, z, metric)
    for metric ∈ (:l1, :l2, :max)
        alg(c; opt = optimizer) =
            Interpol(; optimizer = opt, planes = 5, metric, convexify = c)
        results = Dict(c => approx(fe, Convex(), alg(c)) for c ∈
                       (:auto, :dp, :data, :knots))
        e = Dict(c => error(pwa, metric) for (c, pwa) ∈ results)
        @test all(length(pwa.planes) ≤ 5 for pwa ∈ values(results))

        # The dynamic program interpolates the end points of the data
        @test evaluate(results[:dp], (x[1],)) ≈ z[1]
        @test evaluate(results[:dp], (x[end],)) ≈ z[end]

        # Auto selects the best of :dp and :data, which for these non-convex
        # data is better than adjusting the knots after interpolation
        @test e[:auto] ≈ min(e[:dp], e[:data])
        @test e[:auto] < e[:knots]

        # Without an optimizer, auto uses the dynamic program
        pwa = approx(fe, Convex(), alg(:auto; opt = nothing))
        @test error(pwa, metric) ≈ e[:dp]
    end

    # For convex data all methods give the same result
    zc = exp.(x)
    fec = FunctionEvaluations(Tuple.(x), zc)
    ec = [
        PWA._ipol_error(
            approx(fec, Convex(), Interpol(; optimizer, planes = 4, convexify = c)),
            x,
            zc,
            :l1,
        ) for c ∈ (:auto, :dp, :data, :knots)
    ]
    @test all(e ≈ ec[1] for e ∈ ec)

    # Options requiring an optimizer, and unknown options
    @test Interpol().optimizer === nothing
    @test Interpol().convexify == :auto
    @test_throws ErrorException approx(fe, Convex(), Interpol(; convexify = :data))
    @test_throws ErrorException approx(fe, Convex(), Interpol(; convexify = :knots))
    @test_throws ArgumentError approx(fe, Convex(), Interpol(; convexify = :foo))
end

@testset "Interpolation selection of planes" begin
    rng = StableRNG(5)

    # Paths for all numbers of edges agree with single dynamic programs
    x = sort(rand(rng, 12))
    z = sin.(4 .* x) .+ 0.2 .* randn(rng, 12)
    for metric ∈ (:l1, :l2, :max)
        p = PWA._segment_penalties(x, z, metric)
        op = metric == :max ? max : +
        paths = PWA._dp_all_paths(p, 6, op)
        @test all(paths[s] == PWA._dp_segments(p, s, op)[1] for s ∈ 1:6)

        # Convex paths with at most k edges improve with k, and agree with
        # the single convex dynamic program
        σ = [i < j ? (z[j] - z[i]) / (x[j] - x[i]) : 0.0 for i ∈ 1:12, j ∈ 1:12]
        cpaths = PWA._convex_dp_all_paths(σ, p, 6, op)
        cost(q) = reduce(op, p[q[k], q[k+1]] for k ∈ 1:length(q)-1)
        @test issorted(cost.(cpaths), rev = true)
        @test all(
            cost(cpaths[k]) ≈ cost(PWA._convex_dp_segments(σ, p, k, op)) for
            k ∈ 1:6
        )
        @test all(length(cpaths[k]) ≤ k + 1 for k ∈ 1:6)
    end

    # Noise-free data from 3 planes with kinks at data points (x = 0 and x = 1)
    x = collect(range(-2, 2, length = 41))
    truth = max.(-x, 0.5 .* x, 2 .* x .- 1.5)
    fe = FunctionEvaluations(Tuple.(x), truth)
    alg(; kwargs...) = Interpol(; optimizer, planes = 10, kwargs...)
    for metric ∈ (:l1, :l2, :max)
        sel = approx(fe, Convex(), alg(; metric, selection = :tolerance,
                                       tolerance = 1e-8))
        @test length(sel.planes) == 3
    end
    for metric ∈ (:l1, :l2), s ∈ (:aic, :bic)
        @test length(approx(fe, Convex(), alg(; metric, selection = s)).planes) == 3
    end

    # Noisy data: the selected number of planes is close to the truth, and
    # cross validation is deterministic
    z = truth .+ 0.03 .* randn(rng, length(x))
    fen = FunctionEvaluations(Tuple.(x), z)
    for s ∈ (:bic, :cv), metric ∈ (:l1, :l2)
        k = length(approx(fen, Convex(), alg(; metric, selection = s)).planes)
        @test 3 ≤ k ≤ 6
    end
    a = approx(fen, Convex(), alg(; selection = :cv))
    b = approx(fen, Convex(), alg(; selection = :cv))
    @test a.planes == b.planes

    # Tolerance that cannot be reached gives the maximum number of planes
    pwa = @test_logs (:warn,) approx(
        fen,
        Convex(),
        alg(; planes = 2, selection = :tolerance, tolerance = 1e-6),
    )
    @test length(pwa.planes) == 2

    # The maximum number of planes is capped by the number of points
    few = FunctionEvaluations([(0.0,), (1.0,), (3.0,)], [1.0, 0.0, 2.0])
    pwa = approx(few, Convex(), alg(; planes = 20, selection = :bic))
    @test length(pwa.planes) ≤ 2

    # Invalid options
    @test_throws ArgumentError approx(fe, Convex(), alg(; selection = :tolerance))
    @test_throws ArgumentError approx(
        fe,
        Convex(),
        alg(; metric = :max, selection = :aic),
    )
    @test_throws ArgumentError approx(fe, Convex(), alg(; selection = :cv, folds = 1))
    @test_throws ArgumentError approx(fe, Convex(), alg(; selection = :foo))
end

@testset "Interpolation adjusting knots" begin
    rng = StableRNG(9)

    # Nonnegative least squares satisfies the optimality conditions
    for _ ∈ 1:5
        A = randn(rng, 30, 6)
        b = randn(rng, 30)
        c = PWA._nnls(A, b)
        w = A' * (b - A * c)
        @test all(c .≥ 0)
        @test all(w[c.>0] .≈ 0 .|| abs.(w[c.>0]) .< 1e-8)
        @test all(w[c.==0] .< 1e-8)
    end

    # Refitted values at fixed breakpoints are convex
    x = collect(range(0, 3, length = 60))
    z = sin.(3 .* x)
    t = [0.0, 0.7, 1.5, 2.2, 3.0]
    for metric ∈ (:l1, :l2, :max)
        v = PWA._fit_values(x, z, t, metric, optimizer)
        @test PWA._isconvex(t, v)
    end

    # Adjusting never increases the deviation from the data
    fe = FunctionEvaluations(Tuple.(x), z)
    for metric ∈ (:l1, :l2, :max)
        e = map((:none, :values, :breakpoints)) do adjust
            alg = Interpol(; optimizer, planes = 4, metric, adjust)
            return PWA._ipol_error(approx(fe, Convex(), alg), x, z, metric)
        end
        @test e[2] ≤ e[1] + 1e-9
        @test e[3] ≤ e[2] + 1e-9
    end

    # For the maximum deviation from a convex function, refitting the values
    # halves the error of the interpolant (shifting it down)
    x = collect(range(-2, 2, length = 101))
    fe = FunctionEvaluations(Tuple.(x), x .^ 2)
    e = map((:none, :values)) do adjust
        alg = Interpol(; optimizer, planes = 3, metric = :max, adjust)
        return PWA._ipol_error(approx(fe, Convex(), alg), x, x .^ 2, :max)
    end
    @test e[2] ≈ e[1] / 2 rtol = 1e-6

    # Moving the breakpoints recovers kinks between the data points
    truth = max.(-x, 0.4 .* x .+ 0.3, 2 .* x .- 1.2)
    fe = FunctionEvaluations(Tuple.(x), truth)
    e = map((:none, :breakpoints)) do adjust
        alg = Interpol(; planes = 3, metric = :l2, adjust)
        return PWA._ipol_error(approx(fe, Convex(), alg), x, truth, :l2)
    end
    @test e[2] < 1e-3 * e[1]

    # Selection uses adjusted fits, which can place kinks between data points.
    # For noisy data from 3 planes, this selects fewer and more correct planes.
    # Adjusting with :l2 does not require an optimizer.
    k = map((:none, :values)) do adjust
        alg = Interpol(; planes = 10, metric = :l2, selection = :bic, adjust)
        return map(1:8) do seed
            zₙ = truth .+ 0.02 .* randn(StableRNG(seed), length(x))
            fe = FunctionEvaluations(Tuple.(x), zₙ)
            return length(approx(fe, Convex(), alg).planes)
        end
    end
    @test count(==(3), k[2]) > count(==(3), k[1])
    @test sum(k[2]) < sum(k[1])
    @test_throws ErrorException approx(
        fe,
        Convex(),
        Interpol(; planes = 3, metric = :l1, adjust = :values),
    )
    @test_throws ArgumentError approx(
        fe,
        Convex(),
        Interpol(; optimizer, planes = 3, adjust = :foo),
    )
end
