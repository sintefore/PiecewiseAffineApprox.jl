# The contract of the algorithms, see the documentation of ApproximationAlgorithm: every
# algorithm is run on the same cases, checking the properties that the rest of
# the package relies on.

# An algorithm without _fit, and a minimal one: the tangent planes at a few
# points of one-dimensional data
struct NoFit <: PWA.ApproximationAlgorithm end
struct Tangents <: PWA.ApproximationAlgorithm
    planes::Int
end
function PWA._fit(a::Tangents, data::FunctionEvaluations{1})
    o = sortperm(first.(data.points))
    idx = round.(Int, range(1, length(o); length = a.planes + 2))[2:(end-1)]
    planes = map(idx) do i
        l, r = o[max(i - 1, 1)], o[min(i + 1, length(o))]
        slope =
            (data.values[r] - data.values[l]) /
            (data.points[r][1] - data.points[l][1])
        return Plane((slope,), data.values[o[i]] - slope * data.points[o[i]][1])
    end
    return PWAFunc{Convex,1}(planes)
end

@testset "ApproximationAlgorithm contract" begin
    x = collect(range(-1, 1; length = 20))
    f1 = FunctionEvaluations(x, x .^ 2)
    s(x) = x[1]^2 + x[2]^2
    f2 = PWA._sample_uniform(s, [(-1, 1), (-1, 1)], 6)
    err(f, pwa) = approximation_error(f, pwa, :max)

    algorithms = [
        Cluster(planes = 3, rng = StableRNG(1)),
        Cluster(
            planes = 3,
            metric = :l1,
            strict = :outer;
            optimizer,
            rng = StableRNG(1),
        ),
        Interpol(planes = 3),
        Interpol(planes = 3; optimizer, metric = :max),
        MILP(planes = 3; optimizer),
        MILP(planes = 3, strict = :inner; optimizer),
        Progressive(; optimizer, tolerance = 0.1),
        FullOrder(; optimizer),
    ]
    for a ∈ algorithms
        @testset "$(nameof(typeof(a)))" begin
            # Traits read the options
            @test PWA.planes(a) ===
                  (hasfield(typeof(a), :planes) ? a.planes : nothing)
            @test PWA.concave_planes(a) === PWA.planes(a)
            @test PWA.metric(a) === a.metric
            @test PWA.strict(a) ===
                  (hasfield(typeof(a), :strict) ? a.strict : :none)
            @test PWA.optimizer(a) === a.optimizer
            @test PWA.requires_convex_data(a) ==
                  (a isa Union{Progressive,FullOrder})
            @test (1 ∈ PWA.supported_dimensions(a))
            @test (2 ∈ PWA.supported_dimensions(a)) == !(a isa Interpol)

            # A convex approximation with at most planes planes
            pwa = approx(f1, Convex(), a)
            @test pwa isa PWAFunc{Convex,1}
            p = PWA.planes(a)
            isnothing(p) || @test PWA._planes(pwa) ≤ p
            @test err(pwa, f1) < 0.3
            # The strict guarantee at every point
            if PWA.strict(a) == :outer
                @test all(evaluate(pwa, x) ≤ z + 1e-6 for (x, z) ∈ f1)
            elseif PWA.strict(a) == :inner
                @test all(evaluate(pwa, x) ≥ z - 1e-6 for (x, z) ∈ f1)
            end
            # Concave data by negation (the rng of Cluster is advanced by
            # each call, so a fresh one is used for the same result)
            fresh(a) =
                hasfield(typeof(a), :rng) ? PWA.with(a; rng = StableRNG(1)) : a
            cc = approx(FunctionEvaluations(x, -x .^ 2), Concave(), fresh(a))
            @test cc isa PWAFunc{Concave,1}
            @test sort(cc.planes; by = p -> p.α) ==
                  sort(pwa.planes; by = p -> p.α)
            # Deterministic
            @test approx(f1, Convex(), fresh(a)) == pwa

            # Two variables, or a clear error
            if 2 ∈ PWA.supported_dimensions(a)
                pwa2 = approx(f2, Convex(), a)
                @test pwa2 isa PWAFunc{Convex,2}
                @test err(pwa2, f2) < 1.0
            else
                @test_throws ArgumentError approx(f2, Convex(), a)
            end
        end
    end

    # with changes options and rejects unknown ones
    a = Cluster(planes = 3, metric = :l1)
    b = PWA.with(a; planes = 5, strict = :outer)
    @test b.planes == 5 && b.strict == :outer && b.metric == :l1
    @test PWA.with(a) == a
    @test_throws ArgumentError PWA.with(a; folds = 3)
    @test PWA.with(Interpol(planes = 4); selection = :aic).selection == :aic

    # Empty data is rejected before the algorithm is called
    @test_throws ArgumentError approx(
        FunctionEvaluations(Float64[], Float64[]),
        Convex(),
        Cluster(),
    )

    # A new algorithm needs _fit; the trait defaults work for any struct
    @test_throws ArgumentError approx(f1, Convex(), NoFit())
    @test PWA.planes(NoFit()) === nothing
    @test PWA.strict(NoFit()) == :none
    @test PWA.optimizer(NoFit()) === nothing
    @test PWA._default_samples(NoFit(), 1) == 100

    # A minimal algorithm, see the top of the file
    t = approx(f1, Convex(), Tangents(3))
    @test PWA._planes(t) == 3
    @test err(t, f1) < 0.3
    # ... which then also works for concave and difference of convex data,
    # and for sampling a function
    @test approx(FunctionEvaluations(x, -x .^ 2), Concave(), Tangents(3)) isa
          PWAFunc{Concave,1}
    dc = approx(FunctionEvaluations(x, sin.(3x)), DC(), Tangents(3))
    @test dc isa DCFunc{1}
    @test PWA._planes(dc.convex) == 3 && PWA._planes(dc.concave) == 3
    @test approx(x -> x[1]^2, [(-1, 1)], Convex(), Tangents(2)) isa
          PWAFunc{Convex,1}
    @test_throws ArgumentError approx(f2, Convex(), Tangents(3))
end

@testset "Select" begin
    # Noisy samples of a quadratic: the information criteria pick a small
    # number of planes, the tolerance the fewest that reach it
    rng = StableRNG(5)
    x = collect(range(-1, 1; length = 60))
    z = x .^ 2 .+ 0.05 .* randn(rng, 60)
    f = FunctionEvaluations(x, z)
    err(pwa, metric = :l1) = approximation_error(pwa, f, metric)

    s = Select(Interpol(planes = 8); criterion = :bic)
    @test s isa Select{Interpol{Nothing}}
    @test PWA.planes(s) == 8 && PWA.metric(s) == :l1 && PWA.strict(s) == :none
    @test PWA.supported_dimensions(s) == 1:1
    @test occursin("criterion = :bic", repr(s))
    pwa = approx(f, Convex(), s)
    @test 2 ≤ PWA._planes(pwa) < 8
    # The same choice as the selection built into Interpol
    @test pwa == approx(f, Convex(), Interpol(planes = 8, selection = :bic))
    pwa = approx(f, Convex(), Select(Interpol(planes = 8); criterion = :aic))
    @test pwa == approx(f, Convex(), Interpol(planes = 8, selection = :aic))
    pwa = approx(
        f,
        Convex(),
        Select(Interpol(planes = 8; optimizer, metric = :l2); criterion = :aic),
    )
    @test pwa == approx(
        f,
        Convex(),
        Interpol(planes = 8, selection = :aic; optimizer, metric = :l2),
    )

    # Works with any algorithm with a planes option
    s = Select(
        Cluster(planes = 6, rng = StableRNG(1));
        criterion = :tolerance,
        tolerance = 0.05,
    )
    pwa = approx(f, Convex(), s)
    @test err(pwa) ≤ 0.05
    @test PWA._planes(pwa) ≤ 6
    pwa = approx(
        f,
        Convex(),
        Select(MILP(planes = 4; optimizer); criterion = :cv, folds = 4),
    )
    @test 1 ≤ PWA._planes(pwa) ≤ 4
    @test err(pwa) < 0.05
    # Two variables
    s2(x) = x[1]^2 + x[2]^2
    f2 = PWA._sample_uniform(s2, [(-1, 1), (-1, 1)], 7)
    pwa2 = approx(
        f2,
        Convex(),
        Select(Cluster(planes = 6, rng = StableRNG(1)); criterion = :bic),
    )
    @test pwa2 isa PWAFunc{Convex,2}
    @test 2 ≤ PWA._planes(pwa2) ≤ 6

    # Options are changed on the wrapper or forwarded to the inner algorithm
    s = Select(Cluster(planes = 6); criterion = :aic)
    @test PWA.with(s; planes = 3).inner.planes == 3
    @test PWA.with(s; criterion = :bic).criterion == :bic
    @test PWA.with(s; strict = :outer).inner.strict == :outer
    @test_throws ArgumentError PWA.with(s; foo = 1)
    @test PWA.concave_planes(Select(Cluster(planes = 6, concave_planes = 2))) ==
          2
    # The tolerance is warned about when not reached
    @test_logs (:warn, r"not reached") approx(
        f,
        Convex(),
        Select(Interpol(planes = 2); criterion = :tolerance, tolerance = 1e-6),
    )
    # Concave and difference of convex data through the generic paths
    @test approx(
        FunctionEvaluations(x, -z),
        Concave(),
        Select(Interpol(planes = 8)),
    ) isa PWAFunc{Concave,1}
    dc = approx(
        FunctionEvaluations(x, sin.(3x)),
        DC(),
        Select(Interpol(planes = 6); criterion = :bic),
    )
    @test dc isa DCFunc{1}
    @test PWA._planes(dc.convex) ≤ 6 && PWA._planes(dc.concave) ≤ 6

    # The selection option of Interpol is the wrapper, including repeated
    # measurements in cross validation and the adjustment of the knots
    xd = repeat(collect(range(-1, 1; length = 20)), 2)
    fd = FunctionEvaluations(xd, xd .^ 2 .+ 0.05 .* randn(StableRNG(7), 40))
    for (alg, sel) ∈ (
        (
            Interpol(planes = 8, selection = :cv),
            Select(Interpol(planes = 8); criterion = :cv),
        ),
        (
            Interpol(
                planes = 8,
                selection = :aic,
                adjust = :values,
                metric = :l2,
            ),
            Select(
                Interpol(planes = 8, adjust = :values, metric = :l2);
                criterion = :aic,
            ),
        ),
        (
            Interpol(planes = 8, selection = :tolerance, tolerance = 0.06),
            Select(
                Interpol(planes = 8);
                criterion = :tolerance,
                tolerance = 0.06,
            ),
        ),
    )
        @test approx(fd, Convex(), alg) == approx(fd, Convex(), sel)
    end
    a = Interpol(
        planes = 8,
        selection = :aic,
        adjust = :breakpoints,
        metric = :l2,
    )
    pwa = approx(fd, Convex(), a)
    k = PWA._planes(
        approx(
            fd,
            Convex(),
            Select(
                Interpol(planes = 8, adjust = :values, metric = :l2);
                criterion = :aic,
            ),
        ),
    )
    @test PWA._planes(pwa) == k
    @test pwa == approx(
        fd,
        Convex(),
        Interpol(planes = k, adjust = :breakpoints, metric = :l2),
    )

    # Invalid options
    @test_throws ArgumentError approx(
        f,
        Convex(),
        Select(Progressive(; optimizer)),
    )
    @test_throws ArgumentError approx(
        f,
        Convex(),
        Select(Interpol(); criterion = :foo),
    )
    @test_throws ArgumentError approx(
        f,
        Convex(),
        Select(Interpol(); criterion = :tolerance),
    )
    @test_throws ArgumentError approx(
        f,
        Convex(),
        Select(Interpol(); criterion = :cv, folds = 1),
    )
    @test_throws ArgumentError approx(
        f,
        Convex(),
        Select(Interpol(); criterion = :aic, metric = :max),
    )
    @test_throws ArgumentError approx(
        f,
        Convex(),
        Select(Interpol(); metric = :foo),
    )
end

@testset "Restart" begin
    x = collect(range(-1, 1; length = 40))
    f = FunctionEvaluations(x, x .^ 2 .+ 0.05 .* randn(StableRNG(5), 40))
    err(pwa) = approximation_error(pwa, f, :l2)

    # The trials of Cluster are the wrapper, with identical results
    a = Cluster(planes = 4, trials = 8, rng = StableRNG(3))
    r = Restart(Cluster(planes = 4, trials = 1); trials = 8, rng = StableRNG(3))
    @test r isa Restart
    @test PWA.planes(r) == 4 && PWA.metric(r) == :l2 && PWA.strict(r) == :none
    @test approx(f, Convex(), r) == approx(f, Convex(), a)
    # More trials do not give a worse result with the same first seeds
    single =
        approx(f, Convex(), Cluster(planes = 4, trials = 1, rng = StableRNG(3)))
    @test err(
        approx(
            f,
            Convex(),
            Restart(
                Cluster(planes = 4, trials = 1);
                trials = 8,
                rng = StableRNG(3),
            ),
        ),
    ) ≤ err(single) + 1e-12
    @test occursin("trials = 8", repr(r))

    # Restarts of a wrapped algorithm, and of other metrics and strictness
    s = Restart(
        Select(Cluster(planes = 6, trials = 1); criterion = :bic);
        trials = 4,
        rng = StableRNG(1),
    )
    pwa = approx(f, Convex(), s)
    @test PWA._planes(pwa) ≤ 6
    @test err(pwa) < 0.1
    o = Restart(
        Cluster(
            planes = 4,
            trials = 1,
            metric = :l1,
            strict = :outer;
            optimizer,
        );
        trials = 3,
        rng = StableRNG(1),
    )
    pwa = approx(f, Convex(), o)
    @test all(evaluate(pwa, p) ≤ v + 1e-6 for (p, v) ∈ f)
    @test PWA.metric(o) == :l1 && PWA.strict(o) == :outer
    # Through the concave and difference of convex paths
    @test approx(FunctionEvaluations(x, -x .^ 2), Concave(), r) isa
          PWAFunc{Concave,1}
    @test approx(FunctionEvaluations(x, sin.(3x)), DC(), r) isa DCFunc{1}

    # Options go to the wrapper or the inner algorithm
    @test PWA.with(r; trials = 3).trials == 3
    @test PWA.with(r; planes = 2).inner.planes == 2
    @test PWA.with(r; rng = StableRNG(9)).rng isa StableRNG
    @test_throws ArgumentError PWA.with(r; foo = 1)

    # An algorithm without randomness cannot be restarted
    @test_throws ArgumentError approx(
        f,
        Convex(),
        Restart(Interpol(planes = 3)),
    )
    @test_throws ArgumentError approx(
        f,
        Convex(),
        Restart(Cluster(); trials = 0),
    )
    @test_throws ArgumentError approx(
        f,
        Convex(),
        Restart(Cluster(); metric = :foo),
    )
    # Errors in a trial are reported directly
    failing = () -> error("no solver available")
    @test_throws ErrorException approx(
        f,
        Convex(),
        Restart(
            Cluster(trials = 1, metric = :l1, optimizer = failing);
            trials = 2,
        ),
    )
end

@testset "fit and diagnostics" begin
    x = collect(range(-1, 1; length = 30))
    f = FunctionEvaluations(x, x .^ 2 .+ 0.05 .* randn(StableRNG(5), 30))

    r = fit(f, Convex(), Cluster(planes = 3, rng = StableRNG(1)))
    @test r isa FitResult
    @test r.func == approx(f, Convex(), Cluster(planes = 3, rng = StableRNG(1)))
    @test r.error == approximation_error(r.func, f, :l2)
    @test r.metric == :l2
    @test r.time > 0
    @test r.algorithm isa Cluster
    @test r.details.trials == 20 && length(r.details.errors) == 20
    @test r.details.errors[r.details.best] == minimum(r.details.errors)
    @test occursin("error = ", repr(r)) && occursin("trials", repr(r))

    r = fit(f, Convex(), MILP(planes = 2; optimizer))
    @test r.details.status == MOI.OPTIMAL
    @test r.details.objective ≥ 0 && r.details.big_M_rounds ≥ 0
    @test r.metric == :l1
    r = fit(f, Convex(), Interpol(planes = 3))
    @test isempty(r.details)
    q = FunctionEvaluations(x, x .^ 2)
    r = fit(q, Convex(), Progressive(; optimizer, tolerance = 0.05))
    @test r.details.rounds ≥ 0 && r.details.planes == PWA._planes(r.func)
    r = fit(f, Convex(), Select(Interpol(planes = 6); criterion = :bic))
    @test r.details.planes == PWA._planes(r.func)
    @test length(r.details.scores) == 6 &&
          argmin(r.details.scores) == r.details.planes
    r = fit(
        f,
        Convex(),
        Restart(
            Select(Cluster(planes = 5, trials = 1); criterion = :aic);
            trials = 3,
            rng = StableRNG(1),
        ),
    )
    @test r.details.trials == 3 && haskey(r.details.inner, :scores)

    # Concave, difference of convex, and sampled functions
    r = fit(FunctionEvaluations(x, -x .^ 2), Concave(), Interpol(planes = 3))
    @test r.func isa PWAFunc{Concave,1} &&
          r.error ≈
          approximation_error(r.func, FunctionEvaluations(x, -x .^ 2), :l1)
    r = fit(
        FunctionEvaluations(x, sin.(3x)),
        DC(),
        Interpol(planes = 3; optimizer),
    )
    @test r.func isa DCFunc{1}
    @test r.details.refits ≥ 0 &&
          length(r.details.errors) == r.details.refits + 1
    @test r.details.errors[end] ≥ r.error
    r = fit(
        sin,
        [(0, 2π)],
        DC(),
        Alternating(
            convex = Interpol(planes = 3),
            concave = Interpol(planes = 2),
        ),
    )
    @test r.algorithm isa Alternating && haskey(r.details, :refits)
    r = fit(x -> x[1]^2, [(-1, 1)], Convex(), Interpol(planes = 2))
    @test domain(r.func) == [(-1.0, 1.0)]
end

@testset "Option validation" begin
    x = collect(range(-1, 1; length = 20))
    f = FunctionEvaluations(x, x .^ 2)
    @test PWA.supported_metrics(Cluster()) == (:l1, :l2, :rms, :max)
    @test PWA.supported_metrics(MILP(; optimizer)) == (:l1, :l2, :max)
    @test :none ∈ PWA.supported_metrics(FullOrder(; optimizer))
    @test PWA.supported_strict(Cluster()) == (:none, :inner, :outer)
    # Checked before the algorithm is called, for every algorithm
    for a ∈ (
        Cluster(metric = :foo),
        Cluster(strict = :bar; optimizer),
        MILP(metric = :rms; optimizer),
        MILP(strict = :foo; optimizer),
        Interpol(metric = :foo),
        Progressive(metric = :foo; optimizer),
        Progressive(metric = :foo, strict = :none; optimizer),
        Progressive(strict = :foo; optimizer),
        FullOrder(metric = :foo; optimizer),
        Select(Interpol(); metric = :foo),
        Restart(Cluster(); metric = :foo),
    )
        @test_throws ArgumentError approx(f, Convex(), a)
        @test_throws ArgumentError approx(f, Concave(), a)
        @test_throws ArgumentError approx(f, DC(), a)
    end
    @test approx(f, Convex(), FullOrder(metric = :none; optimizer)) isa PWAFunc
end

@testset "Normalization" begin
    rng = StableRNG(11)
    X = [3.0 .+ 2.0 .* rand(rng, 30)'; 100.0 .+ 50.0 .* rand(rng, 30)']
    z = 1e4 .+ 1e3 .* vec(sum(X .^ 2; dims = 1))
    data = FunctionEvaluations(X, z; dims = 2)
    for whiten ∈ (false, true), values ∈ (true, false)
        nd, n = PWA._normalize(data; whiten, values)
        Xn = PWA._matrix(nd)
        if whiten
            @test PWA.cov(Xn; dims = 2) ≈ I atol = 1e-8
        else
            @test all(-1 .≤ Xn .≤ 1)
            @test vec(minimum(Xn; dims = 2)) ≈ [-1.0, -1.0]
            @test vec(maximum(Xn; dims = 2)) ≈ [1.0, 1.0]
        end
        if values
            @test extrema(nd.values) == (-1.0, 1.0)
        else
            @test nd.values == z
        end
        # A plane fitted to the normalized data is the same function after
        # transforming it back
        p = Plane((0.3, -0.7), 0.2)
        q = PWA._denormalize(p, n)
        @test all(
            isapprox(evaluate(q, x), n.ν + n.s * evaluate(p, x̃); rtol = 1e-10)
            for (x, x̃) ∈ zip(data.points, nd.points)
        )
        pwa =
            PWA._denormalize(PWAFunc{Convex,2}([p, Plane((1.0, 1.0), 0.0)]), n)
        @test pwa isa PWAFunc{Convex,2}
    end
    # Constant coordinates and values are shifted only, or scaled by the floor
    c = FunctionEvaluations(
        [(1.0, 5.0), (2.0, 5.0), (3.0, 5.0)],
        [4.0, 4.0, 4.0],
    )
    nd, n = PWA._normalize(c)
    @test PWA._matrix(nd)[2, :] == zeros(3) && nd.values == zeros(3) && n.s == 1
    nd, n = PWA._normalize(c; floor = 0.5)
    @test n.s == 0.5
    # The algorithms are invariant to the scale and offset of the data
    x = collect(range(0, 1; length = 25))
    f = FunctionEvaluations(x, x .^ 2 .+ 0.02 .* randn(StableRNG(2), 25))
    g = FunctionEvaluations(1e3 .* x .- 5e2, 1e-3 .* f.values .+ 7.0)
    for a ∈ (
        Interpol(planes = 3),
        MILP(planes = 3; optimizer),
        Cluster(planes = 3, rng = StableRNG(1)),
        Progressive(; optimizer, tolerance = 0.08, strict = :none),
    )
        pf = approx(f, Convex(), a)
        pg = approx(
            g,
            Convex(),
            PWA.with(
                a;
                (hasfield(typeof(a), :rng) ? (rng = StableRNG(1),) : (;))...,
                (
                    hasfield(typeof(a), :tolerance) ? (tolerance = 8e-5,) : (;)
                )...,
            ),
        )
        @test all(
            isapprox(
                evaluate(pg, (1e3 * t - 5e2,)),
                1e-3 * evaluate(pf, (t,)) + 7.0;
                atol = 1e-6,
            ) for t ∈ x
        )
    end
end
