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
        pts = vec([(a, b, c) for a ∈ g, b ∈ g, c ∈ (-1.0, 1.0)])
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

@testset "Quad Approx" begin
    points = [i => i^2 for i ∈ 0:0.1:1]
    data = FunctionEvaluations([(p.first,) for p ∈ points], [p.second for p ∈ points])
    for alg ∈ (
        MILP(optimizer = optimizer, metric = :l1, planes = 5),
        Cluster(metric = :l2, planes = 5, rng = StableRNG(1)),
    )
        pwa = approx(data, Convex(), alg)
        # Added to a model with a new output variable
        m = Model(optimizer)
        @variable(m, x)
        y = pwa(tuple(x))
        @objective(m, Min, y)
        @constraint(m, x >= 0.3)
        optimize!(m)
        @test isapprox(value(y), 0.09, atol = 0.015)
        @constraint(m, x >= 0.9)
        optimize!(m)
        @test isapprox(value(y), 0.81, atol = 0.1)

        # Added using an existing output variable
        m = Model(optimizer)
        @variable(m, x)
        @variable(m, test_y)
        pwa(tuple(x); z = test_y)
        @objective(m, Min, test_y)
        @constraint(m, x >= 0.3)
        optimize!(m)
        @test isapprox(value(test_y), 0.09, atol = 0.015)
    end
end

@testset "2D Convex Approx" begin
    xg = -1:0.5:1
    X = [repeat(xg, inner = length(xg)) repeat(xg, outer = length(xg))]
    points = Tuple.(eachrow(X))
    z = X[:, 1] .^ 2 + X[:, 2] .^ 2
    alg = MILP(optimizer = optimizer, planes = 17, strict = :outer, metric = :l1)
    pwa1 = approx(FunctionEvaluations(points, z), Convex(), alg)
    pwa2 = approx(FunctionEvaluations(points, -z), Concave(), alg)

    @test isapprox(PWA.evaluate(pwa1, (0.5, 0.5)), 0.5, atol = 0.1)
    @test isapprox(
        PWA.evaluate(pwa1, (-0.3, 0.4)),
        -PWA.evaluate(pwa2, (-0.3, 0.4)),
        atol = 0.01,
    )

    # Added to models with existing variables: a tuple of inputs and an output
    function solve(pwa, sense)
        m = Model(optimizer)
        @variable(m, xvar)
        @variable(m, yvar)
        @variable(m, f)
        pwa((xvar, yvar); z = f)
        @objective(m, sense, f)
        @constraint(m, xvar == √0.5)
        @constraint(m, yvar == √0.5)
        optimize!(m)
        return value(xvar), value(yvar), value(f)
    end
    x1, y1, f1 = solve(pwa1, MIN_SENSE)
    x2, y2, f2 = solve(pwa2, MAX_SENSE)
    @test isapprox(f1, 1.0, rtol = 0.12)
    @test isapprox(f1, -f2, atol = 0.01)
    @test isapprox(x1, x2, atol = 0.01)
    @test isapprox(y1, y2, atol = 0.01)
end

@testset "Big-M calculations" begin
    function mat2tuples(x)
        return collect(Tuple(x'[:, i]) for i ∈ 1:size(x', 2))
    end

    # big-M for a 1D function
    x = collect(range(-1, 1; length = 10))
    z = x .^ 2

    M = PWA._linear_big_M(FunctionEvaluations(tuple.(x), z))
    @test M ≈ 3.91 atol = 0.01

    # big-M for a 2D function
    xg = collect(range(-1, 1; length = 4))
    yg = collect(range(-1, 1; length = 4))
    X = [repeat(xg, inner = [size(yg, 1)]) repeat(yg, outer = [size(xg, 1)])]
    z = X[:, 1] .^ 2 + X[:, 2] .^ 2
    np = 5

    M = PWA._linear_big_M(FunctionEvaluations(mat2tuples(X), z))
    @test M ≈ 14.67 atol = 0.01

    @testset "Duplicate inputs" begin
        # 1D
        x = [302.3633312955, 300.6680553405, 274.281686739, 274.281686739]
        z = [-184.9368247, -183.9371662, -167.94263, -167.94263]

        M¹ = PWA._linear_big_M(FunctionEvaluations(tuple.(x), z))
        @test M¹ ≈ 0.479 atol = 0.001

        M² =
            PWA._linear_big_M(FunctionEvaluations(tuple.(unique(x)), unique(z)))
        @test M¹ == M²

        # Duplicate with different function value
        z[4] = -166.3
        M³ = PWA._linear_big_M(FunctionEvaluations(tuple.(x), z))
        @test M³ ≈ 2.286 atol = 0.001

        # 2D
        xg = collect(range(-1, 1; length = 4))
        yg = collect(range(-1, 1; length = 4))

        X =
            [repeat(xg, inner = [size(yg, 1)]) repeat(yg, outer = [size(xg, 1)])]
        z = X[:, 1] .^ 2 + X[:, 2] .^ 2

        # Duplicate point with different function value
        X = vcat(X, [-1 1 / 3])
        z = vcat(z, 1.5)

        M = PWA._linear_big_M(FunctionEvaluations(mat2tuples(X), z))
        @test M ≈ 15.95 atol = 0.01
    end

    @testset "Close points" begin
        # 1D
        x = collect(range(-1, 1; length = 10))
        z = x .^ 2

        # Not close enough
        x1 = vcat(x, 0.55)
        z1 = vcat(z, 0.41)
        M¹ = PWA._linear_big_M(FunctionEvaluations(tuple.(x1), z1))
        @test M¹ ≈ 9.68 atol = 0.001

        # Close points
        x2 = vcat(x, 0.5555)
        z2 = vcat(z, 0.41)
        M² = PWA._linear_big_M(FunctionEvaluations(tuple.(x2), z2))
        @test M² ≈ 3.911 atol = 0.001
    end

    @testset "Scale invariance 2D" begin
        # big-M should scale with the data, not drop to zero for small domains
        g = range(-1, 1; length = 4)
        pts = vec([(a, b) for a ∈ g, b ∈ g])
        z = [a^2 + b^2 for (a, b) ∈ pts]
        M = PWA._linear_big_M(FunctionEvaluations(pts, z))
        for s ∈ (10.0, 0.1, 0.01)
            Mˢ = PWA._linear_big_M(
                FunctionEvaluations([s .* p for p ∈ pts], s^2 .* z),
            )
            @test Mˢ ≈ s^2 * M
        end
    end

    @testset "Same as brute force" begin
        # Direct implementation of the definitions: all segments (1D) or
        # triangles (2D) without data points in the interior
        function brute_force_1d(x, z)
            δ = (maximum(x) - minimum(x)) / (100 * length(x))
            M = zeros(length(x))
            for i ∈ eachindex(x), j ∈ eachindex(x)
                x[j] - x[i] > δ || continue
                any(p -> x[i] + δ < p < x[j] - δ, x) && continue
                a = (z[j] - z[i]) / (x[j] - x[i])
                M .= max.(M, z .- a .* x .- (z[i] - a * x[i]))
            end
            return 1.1 .* M
        end
        function brute_force_2d(pts, z)
            lo, hi = minimum.(zip(pts...)), maximum.(zip(pts...))
            x = [(p .- lo) ./ (hi .- lo) for p ∈ pts]
            M = zeros(length(x))
            for i ∈ eachindex(x), j ∈ i+1:length(x), k ∈ j+1:length(x)
                v = (x[i], x[j], x[k])
                any(p -> p ∉ v && PWA._point_in_triangle(p, v...), x) && continue
                PWA._sliver_parameter(v...) > 0.025 || continue
                t = [(x[l]..., z[l]) for l ∈ (i, j, k)]
                n = cross(collect(t[1] .- t[2]), collect(t[1] .- t[3]))
                M .= max.(M, z .- [PWA._plane_f(p, n, dot(n, t[1])) for p ∈ x])
            end
            return 1.1 .* M
        end

        rng = StableRNG(1)
        for n ∈ (2, 5, 20)
            x = rand(rng, n)
            x = vcat(x, x[1], x[2] + 1e-6)
            z = x .^ 2 .+ 0.1 .* randn(rng, length(x))
            M = PWA._linear_big_M_points(FunctionEvaluations(tuple.(x), z))
            @test M ≈ brute_force_1d(x, z)

            pts = [(rand(rng), rand(rng)) for _ ∈ 1:n+1]
            pts = vcat(pts, pts[1:1], [(0.5, 0.5), (0.5, 0.25), (0.5, 0.75)])
            z = [a^2 + b^2 + 0.1 * randn(rng) for (a, b) ∈ pts]
            M = PWA._linear_big_M_points(FunctionEvaluations(pts, z))
            @test M ≈ brute_force_2d(pts, z)
        end
    end
end
