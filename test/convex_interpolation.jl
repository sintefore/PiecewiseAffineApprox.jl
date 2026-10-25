
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

    # Test convexification of the knot values
    v = PWA._convexify_values(pwa.x, pwa.z, :l1, optimizer)
    @test length(v) == 6
    @test all(diff(diff(v) ./ diff(pwa.x)) .≥ -1e-8)
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

    # Repeated identical measurements do not change the largest deviation, so
    # the fit with metric = :max is unchanged (with :l1 and :l2 they count
    # with their number)
    xd = vcat(x, x[[3, 8]])
    zd = vcat(z, z[[3, 8]])
    algₘ = Interpol(; optimizer, planes = 4, metric = :max)
    refₘ = approx(FunctionEvaluations(Tuple.(x), z), Convex(), algₘ)
    dup = approx(FunctionEvaluations(Tuple.(xd), zd), Convex(), algₘ)
    @test all(evaluate(dup, (t,)) ≈ evaluate(refₘ, (t,)) for t ∈ x)

    # Duplicates with different values are replaced by their average
    d = PWA._sort_unique([2.0, 1.0, 2.0, 0.0, 2.0], [1.0, 5.0, 2.0, 7.0, 6.0])
    @test d.x == [0.0, 1.0, 2.0]
    @test d.z ≈ [7.0, 5.0, 3.0]
    @test d.n == [1, 1, 3]
    @test d.vals[3] == [1.0, 2.0, 6.0]

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
                       (:auto, :dp, :data))
        e = Dict(c => error(pwa, metric) for (c, pwa) ∈ results)
        @test all(length(pwa.planes) ≤ 5 for pwa ∈ values(results))

        # The dynamic program interpolates the end points of the data
        @test evaluate(results[:dp], (x[1],)) ≈ z[1]
        @test evaluate(results[:dp], (x[end],)) ≈ z[end]

        # Auto selects the best of :dp and :data
        @test e[:auto] ≈ min(e[:dp], e[:data])

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
        ) for c ∈ (:auto, :dp, :data)
    ]
    @test all(e ≈ ec[1] for e ∈ ec)

    # Options requiring an optimizer, and unknown options
    @test Interpol().optimizer === nothing
    @test Interpol().convexify == :auto
    @test_throws ErrorException approx(fe, Convex(), Interpol(; convexify = :data))
    @test_throws ArgumentError approx(fe, Convex(), Interpol(; convexify = :knots))
    @test_throws ArgumentError approx(fe, Convex(), Interpol(; convexify = :foo))
end

@testset "Interpolation selection of planes" begin
    rng = StableRNG(5)
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

@testset "Interpolation with repeated measurements" begin
    # Unequal numbers of measurements per load, as in a test campaign
    rng = StableRNG(1)
    counts = [1, 1, 1, 2, 3, 6, 12, 12, 12]
    loads = reduce(vcat, fill(l, c) for (l, c) ∈ zip(0.2:0.1:1.0, counts))
    z = [155 + 12 / l + 60(l - 0.7)^2 + 3randn(rng) for l ∈ loads]
    d = PWA._sort_unique(loads, z)
    @test d.n == counts

    # Penalties and errors are measured against every measurement
    for metric ∈ (:l1, :l2, :max)
        p = PWA._segment_penalties(d, metric)
        for (i, j) ∈ ((1, 9), (3, 7), (6, 9))
            c = (d.z[j] - d.z[i]) / (d.x[j] - d.x[i])
            r = [abs(d.z[i] + c * (d.x[g] - d.x[i]) - v) for g ∈ i:j-1 for v ∈ d.vals[g]]
            @test p[i, j] ≈ PWA._aggregate(r, metric)
        end
        knots = PWA.PWAFunction1D(d.x[[1, 5, 9]], d.z[[1, 5, 9]])
        f = PWA._interpolate(knots.x, knots.z, loads)
        @test PWA._knots_error(knots, d, metric) ≈ PWA._aggregate(abs.(f .- z), metric)
    end

    # The fit is better on the measurements than a fit to the means per load
    fe = FunctionEvaluations(tuple.(loads), z)
    means = FunctionEvaluations(tuple.(d.x), d.z)
    sse(pwa) = sum((evaluate(pwa, x) - v)^2 for (x, v) ∈ fe)
    maxdev(pwa) = maximum(abs(evaluate(pwa, x) - v) for (x, v) ∈ fe)
    for (metric, err) ∈ ((:l2, sse), (:max, maxdev))
        alg = Interpol(; optimizer, planes = 3, metric, convexify = :dp, adjust = :values)
        @test err(approx(fe, Convex(), alg)) < err(approx(means, Convex(), alg))
    end

    # Convexifying the data and refitting values account for all measurements
    v = PWA._convexify_values(d, :l2, optimizer)
    @test all(diff(diff(v) ./ diff(d.x)) .≥ -1e-8)
    vfit = PWA._fit_values(d, d.x[[1, 5, 9]], :l1, optimizer)
    @test length(vfit) == 3
end

@testset "Interpolation selection at the maximum" begin
    x = collect(range(-1, 1, length = 40))
    truth = max.(-x, 0.3x .+ 0.2, 1.5x .- 0.8, 3x .- 2.5, -2x .- 1.2)
    fe = FunctionEvaluations(tuple.(x), truth .+ 0.005 .* randn(StableRNG(3), 40))
    alg(planes) = Interpol(; optimizer, planes, metric = :l2, selection = :bic)
    @test_logs (:warn, r"maximum number of planes") approx(fe, Convex(), alg(4))
    pwa = @test_logs min_level = Logging.Warn approx(fe, Convex(), alg(10))
    @test length(pwa.planes) > 4
end

@testset "Interpolation validation and near-duplicate points" begin
    x = collect(range(-1, 1, length = 20))
    fe = FunctionEvaluations(Tuple.(x), x .^ 2)
    for kw ∈ (
        (metric = :foo,),
        (convexify = :foo,),
        (selection = :foo,),
        (adjust = :foo,),
        (planes = 0,),
        (selection = :cv, folds = 1),
        (selection = :tolerance,),
    )
        @test_throws ArgumentError approx(fe, Convex(), Interpol(; optimizer, kw...))
    end
    # :rms is accepted as :l2
    a = approx(fe, Convex(), Interpol(; planes = 3, metric = :rms))
    b = approx(fe, Convex(), Interpol(; planes = 3, metric = :l2))
    @test a.planes == b.planes

    # Measurements at indistinguishably close points (within 1e-6 of the
    # range) are grouped, giving the same result as exact repetitions, and
    # points above the grouping tolerance are solved with well scaled
    # constraints
    rng = StableRNG(1)
    loads = repeat(collect(0.2:0.1:1.0), inner = 5)
    zl = [155 + 12 / l + 60(l - 0.7)^2 + 2randn(rng) for l ∈ loads]
    exact = FunctionEvaluations(tuple.(loads), zl)
    for metric ∈ (:l2, :l1, :max)
        alg = Interpol(; optimizer, planes = 3, metric, convexify = :data, adjust = :values)
        p1 = approx(exact, Convex(), alg)
        close = FunctionEvaluations(tuple.(loads .+ 1e-9 .* repeat(0:4, outer = 9)), zl)
        p2 = approx(close, Convex(), alg)
        @test all(evaluate(p1, (t,)) ≈ evaluate(p2, (t,)) for t ∈ 0.2:0.1:1.0)
        apart = FunctionEvaluations(tuple.(loads .+ 1e-5 .* repeat(0:4, outer = 9)), zl)
        p3 = approx(apart, Convex(), alg)
        @test length(p3.planes) ≤ 3
    end
end

@testset "Interpolation scale invariance, size guard and candidate reuse" begin
    x = collect(range(-1, 1, length = 40))
    z = exp.(x) .+ 0.01 .* randn(StableRNG(2), 40)
    fe = FunctionEvaluations(Tuple.(x), z)
    scaled(sx, ox, sz) = FunctionEvaluations(Tuple.(sx .* x .+ ox), sz .* z)
    for kw ∈ (
        (convexify = :dp,),
        (convexify = :dp, metric = :l2, adjust = :values),
        (convexify = :data, metric = :l1, adjust = :breakpoints),
        (metric = :l2, selection = :bic, planes = 8),
    )
        alg = Interpol(; optimizer, planes = 5, kw...)
        base = approx(fe, Convex(), alg)
        sse(p, pts, scale) = sum((evaluate(p, (pt,)) / scale - v)^2 for (pt, v) ∈ zip(pts, z))
        for (sx, ox, sz) ∈ ((1e6, 3e7, 1e-13), (1e-4, -2.0, 1e9))
            pwa = approx(scaled(sx, ox, sz), Convex(), alg)
            @test length(pwa.planes) == length(base.planes)
            if get(kw, :adjust, :none) == :breakpoints
                # The search over the breakpoints is a local search that can
                # end in different local optima, so only the error is compared
                @test sse(pwa, sx .* x .+ ox, sz) ≈ sse(base, x, 1.0) rtol = 1e-2
            else
                @test all(
                    isapprox(evaluate(pwa, (sx * t + ox,)) / sz, evaluate(base, (t,));
                        rtol = 1e-6) for t ∈ x
                )
            end
        end
    end

    # The tolerance is in the units of the data (here scaled by 1e-13)
    tol = 3e-15
    alg = Interpol(; planes = 8, metric = :max, selection = :tolerance, tolerance = tol)
    pwa = approx(scaled(1.0, 0.0, 1e-13), Convex(), alg)
    @test maximum(abs(evaluate(pwa, (t,)) - 1e-13 * v) for (t, v) ∈ zip(x, z)) ≤ tol * (1 + 1e-6)
    @test length(pwa.planes) < 8

    # Data too large for the dynamic program is rejected before any work
    big = FunctionEvaluations(Tuple.(range(0, 1, length = 30_000)), zeros(30_000))
    @test_throws ArgumentError approx(big, Convex(), Interpol(; planes = 50))
    @test_throws ArgumentError approx(big, Convex(), Interpol(; planes = 50, selection = :bic))

    # After selection, the result is the fit with the selected number of planes
    sel = approx(fe, Convex(), Interpol(; planes = 8, metric = :l2, selection = :bic))
    fixed = approx(fe, Convex(), Interpol(; planes = length(sel.planes), metric = :l2))
    @test all(evaluate(sel, (t,)) ≈ evaluate(fixed, (t,)) for t ∈ x)
end
