# The largest deviations on a fine grid, as (min, max) of pwa - f
function grid_deviation(pwa, f, box, n)
    lo, hi = Inf, -Inf
    for x ∈ Iterators.product((range(b[1], b[2], n) for b ∈ box)...)
        d = evaluate(pwa, x) - f(x)
        lo, hi = min(lo, d), max(hi, d)
    end
    return lo, hi
end

@testset "Error bounds over the domain" begin
    # x² on [-1, 1] with the tangents at ±0.5 and the chord-like plane z = 1
    f = x -> x[1]^2
    pwa = PWAFunc{Convex,1}(
        [Plane((1.0,), -0.25), Plane((-1.0,), -0.25)];
        domain = [(-1.0, 1.0)],
    )
    b = error_bounds(pwa, f)
    # Tangents never exceed the function, and are furthest below at x = ±1, 0
    @test b.above ≈ 0 atol = 1e-9
    @test b.below ≈ 0.25
    @test abs(b.x_below[1]) ∈ (0.0, 1.0)
    # Shifted up by 0.1: above by 0.1 at the points of tangency
    up = PWAFunc{Convex,1}([Plane(p.α, p.β + 0.1) for p ∈ pwa.planes])
    b = error_bounds(up, f, [(-1, 1)])
    @test b.above ≈ 0.1 atol = 1e-9
    @test abs(b.x_above[1]) ≈ 0.5 atol = 1e-6
    @test b.below ≈ 0.15
    # Adding the same constant to both functions gives the same bounds. Before,
    # with a large constant, the vertex at 0 was lost and below was -0.5.
    v = PWAFunc{Convex,1}([Plane((1.0,), 0.0), Plane((-1.0,), 0.0)]; domain = [(-1.0, 1.0)])
    for c ∈ (0.0, 1e6, 1e12)
        lifted = PWAFunc{Convex,1}([Plane(p.α, p.β + c) for p ∈ v.planes]; domain = [(-1.0, 1.0)])
        b = error_bounds(lifted, x -> c + 0.5)
        @test b.below ≈ 0.5
        @test b.x_below[1] ≈ 0 atol = 1e-9
        @test b.above ≈ 0.5 atol = 1e-3     # found by a search, near 10¹² the spacing is 10⁻⁴
    end
    @test_throws ArgumentError error_bounds(up, f)
    # The same for the concave mirror image
    down = PWAFunc{Concave,1}(up.planes; domain = [(-1.0, 1.0)])
    b = error_bounds(down, x -> -x[1]^2)
    @test b.below ≈ 0.1 atol = 1e-9
    @test b.above ≈ 0.15

    # Two variables, against a fine grid
    g = (x, y) -> x^2 + 2y^2
    box = [(-1.0, 1.0), (-0.5, 1.0)]
    pwa = approx(g, box, Convex(), Cluster(planes = 6, rng = StableRNG(1)))
    b = error_bounds(pwa, g)
    lo, hi = grid_deviation(pwa, x -> g(x...), box, 301)
    @test b.above ≥ hi - 1e-9 && b.above ≤ hi + 0.01
    @test b.below ≥ -lo - 1e-9 && b.below ≤ -lo + 0.01
end

@testset "One-sided approximations over the whole domain" begin
    cases = (
        (x -> exp(x[1]), [(-1.0, 2.0)], Convex()),
        (x -> log(x[1]), [(0.1, 3.0)], Concave()),
        ((x, y) -> x^2 + 2y^2, [(-1.0, 1.0), (-1.0, 1.0)], Convex()),
        ((x, y) -> sqrt(x * y), [(0.5, 3.0), (0.5, 3.0)], Concave()),
    )
    for (i, (f, box, c)) ∈ enumerate(cases), s ∈ (:outer, :inner)
        fx = x -> f(x...)
        n = length(box) == 1 ? 2001 : 151
        # A new algorithm for each fit, as the fits advance its rng
        alg() = Cluster(trials = 4, planes = 5, metric = :l1, strict = s; optimizer, rng = StableRNG(2))
        # At the samples only, the side is violated between them
        plain = approx(f, box, c, alg())
        r = fit(f, box, c, alg(); everywhere = true, refinements = length(box) == 1 ? 3 : 1)
        lo, hi = grid_deviation(r.func, fx, box, n)
        below = (s == :outer) == (c isa Convex)
        if below
            @test hi ≤ 0
            @test grid_deviation(plain, fx, box, n)[2] > 1e-6
        else
            @test lo ≥ 0
            @test grid_deviation(plain, fx, box, n)[1] < -1e-6
        end
        # The reported bounds are those of the result
        b = r.details.bounds
        @test b.above ≥ hi - 1e-9 && b.below ≥ -lo - 1e-9
        @test (below ? b.above : b.below) ≤ 0
        @test r.details.shift ≥ 0 && 0 ≤ r.details.refinements ≤ 3
        @test domain(r.func) == box
        @test r.func isa PWAFunc{typeof(c)}
        # approx passes the option on (checked once, as it fits again)
        if i == 1 && s == :outer
            @test approx(f, box, c, alg(); everywhere = true, refinements = r.details.refinements) == r.func
        end
    end

    # Without refits, the planes are only shifted
    f = x -> x[1]^4
    a = Progressive(; optimizer, tolerance = 0.01, metric = :max, strict = :inner)
    r = fit(f, [(-1.0, 1.0)], Convex(), a; everywhere = true)
    @test r.details.refinements == 0 && r.details.shift > 0
    @test grid_deviation(r.func, f, [(-1.0, 1.0)], 5001)[1] ≥ 0
    plain = approx(f, [(-1.0, 1.0)], Convex(), a)
    @test [p.α for p ∈ r.func.planes] == [p.α for p ∈ plain.planes]

    # Without a side, only the bounds are added
    cl() = Cluster(planes = 4, rng = StableRNG(1))
    r = fit(f, [(-1.0, 1.0)], Convex(), cl(); everywhere = true)
    @test r.func == approx(f, [(-1.0, 1.0)], Convex(), cl())
    @test r.details.shift == 0 && r.details.refinements == 0
    lo, hi = grid_deviation(r.func, f, [(-1.0, 1.0)], 5001)
    @test r.details.bounds.above ≈ hi atol = 1e-3
    @test r.details.bounds.below ≈ -lo atol = 1e-3

    @test_throws ArgumentError fit(sin, [(0.0, 6.0)], DC(), cl(); everywhere = true)
    @test_throws ArgumentError fit(f, [(-1.0, 1.0)], Convex(), cl(); everywhere = true, refinements = -1)
end
