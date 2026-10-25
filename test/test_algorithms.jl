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
