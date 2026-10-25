@testset "PWA vertex enumeration 1D" begin
    # Test 1: Simple convex 1D function with 3 planes
    # f(x) = max(x + 1, -3x - 2, -2x - 1)
    pwa = PWAFunc{Convex,1}([Plane([1], 1), Plane([-3], -2), Plane([-2], -1)])
    verts = PiecewiseAffineApprox.vertices_by_subsets(pwa)

    # Verify we found vertices
    @test length(verts) >= 2

    # Verify each vertex is on the epigraph boundary
    for (x, z) ∈ verts
        @test abs(z - evaluate(pwa, x)) < 1e-8
        # Verify at least 2 planes are active at each vertex (1D needs 2 planes)
        active_planes = count(p -> abs(evaluate(p, x) - z) < 1e-8, pwa.planes)
        @test active_planes >= 2
    end

    # Test 2: Concave 1D function
    # f(x) = min(-x - 1, 3x + 2, 2x + 1)
    pwa_concave =
        PWAFunc{Concave,1}([Plane([1], 1), Plane([-3], -2), Plane([-2], -1)])
    verts_concave = PWA.vertices_by_subsets(pwa_concave)

    @test length(verts_concave) >= 2
    for (x, z) ∈ verts_concave
        @test abs(z - evaluate(pwa_concave, x)) < 1e-8
    end
end

@testset "PWA vertex enumeration 2D" begin
    # Test 1: Simple convex 2D function
    # f(x,y) = max(x + 2y + 1, -3x + y - 2, x - y + 3, 0.5x + 4y - 2)
    pwa = PWAFunc{Convex,2}([
        Plane([1, 2], 1),
        Plane([-3, 1], -2),
        Plane([1, -1], 3),
        Plane([0.5, 4], -2),
    ])

    verts = PiecewiseAffineApprox.vertices_by_subsets(pwa)

    # Verify we found vertices
    @test length(verts) == 2

    # Verify each vertex lies on the epigraph boundary
    for (x, z) ∈ verts
        @test abs(z - evaluate(pwa, x)) < 1e-8
        # In 2D, at least 3 planes should be active at each vertex
        active_planes = count(
            p -> abs(PiecewiseAffineApprox.evaluate(p, x) - z) < 1e-8,
            pwa.planes,
        )
        @test active_planes >= 3
    end

    # Test 2: Known simple case - 3 planes forming a pyramid
    # f(x,y) = max(x, y, 0)
    pwa_pyramid = PWAFunc{Convex,2}([
        Plane([1, 0], 0),   # z >= x
        Plane([0, 1], 0),   # z >= y
        Plane([0, 0], 0),    # z >= 0
    ])

    verts_pyramid = PiecewiseAffineApprox.vertices_by_subsets(pwa_pyramid)

    # The only vertex should be at the origin (0, 0) with value 0
    @test length(verts_pyramid) == 1
    @test all(abs.(verts_pyramid[1][1]) .< 1e-8)  # x ≈ [0, 0]
    @test abs(verts_pyramid[1][2]) < 1e-8         # z ≈ 0

    # Test 3: Quadratic approximation test
    pwa_quad = JSON.parsefile(
        joinpath(@__DIR__, "pwa_ex", "pwa_quad_2d.json"),
        PWAFunc,
    )
    vertices = PiecewiseAffineApprox.vertices_by_subsets(pwa_quad)

    # Should find multiple vertices around the domain boundary
    @test length(vertices) >= 3

    # All vertices should be on the epigraph boundary
    for (x, z) ∈ vertices
        @test abs(z - evaluate(pwa_quad, x)) < 1e-8
    end
end

@testset "PWA vertex enumeration - edge cases" begin
    # Test 1: Single plane (no vertices - unbounded)
    pwa_single = PWAFunc{Convex,2}([Plane([1, 0], 0)])
    verts_single = PiecewiseAffineApprox.vertices_by_subsets(pwa_single)
    @test length(verts_single) == 0  # No vertices for single plane

    # Test 2: Parallel planes (no vertices - the kink is a line)
    pwa_parallel = PWAFunc{Convex,2}([
        Plane([1, 1], 0),
        Plane([1, 1], 1),
        Plane([0, 0], 0),
    ])
    verts_parallel = PiecewiseAffineApprox.vertices_by_subsets(pwa_parallel)
    @test length(verts_parallel) == 0  # No vertices for parallel planes

    # Test 3: Concave 2D function
    pwa_concave = PWAFunc{Concave,2}([
        Plane([1, 0], 0),
        Plane([0, 1], 0),
        Plane([-1, -1], 2),
    ])

    verts_concave = PiecewiseAffineApprox.vertices_by_subsets(pwa_concave)
    @test length(verts_concave) >= 1

    # Verify vertices are on the hypograph boundary
    for (x, z) ∈ verts_concave
        @test abs(z - evaluate(pwa_concave, x)) < 1e-8
    end

    # Test 4: Three planes meeting in a point give a single vertex
    pwa_conc = PWAFunc{Convex,1}([
        Plane((1,), 0),
        Plane((-1,), 0),
        Plane((0,), 0),
        Plane((2,), -1),
    ])
    verts_conc = PiecewiseAffineApprox.vertices_by_subsets(pwa_conc)
    @test sort(only.(first.(verts_conc))) ≈ [0.0, 1.0]
end

@testset "PWA vertex enumeration with bounds" begin
    # 1D: tangents to x² at ±0.25 and ±0.75, with breakpoints at -0.5, 0 and 0.5
    pwa = PWAFunc{Convex,1}([
        Plane((2t,), -t^2) for t ∈ (-0.75, -0.25, 0.25, 0.75)
    ])
    verts = PiecewiseAffineApprox.vertices_by_subsets(pwa, [(-1, 1)])
    @test sort(only.(first.(verts))) ≈ [-1.0, -0.5, 0.0, 0.5, 1.0]
    for (x, z) ∈ verts
        @test z ≈ evaluate(pwa, x)
    end
    # Bounds between two breakpoints only give the end points
    verts = PiecewiseAffineApprox.vertices_by_subsets(pwa, [(0.2, 0.3)])
    @test sort(only.(first.(verts))) ≈ [0.2, 0.3]
    # Single plane
    verts = PiecewiseAffineApprox.vertices_by_subsets(
        PWAFunc{Convex,1}([Plane((1,), 0)]),
        [(-1, 2)],
    )
    sort!(verts)
    @test only.(first.(verts)) ≈ [-1.0, 2.0]
    @test last.(verts) ≈ [-1.0, 2.0]
    # Concave
    verts = PiecewiseAffineApprox.vertices_by_subsets(
        PWAFunc{Concave,1}(pwa.planes),
        [(-1, 1)],
    )
    @test length(verts) == 5
    for (x, z) ∈ verts
        @test z ≈ -evaluate(pwa, x)
    end

    # 2D: f(x,y) = max(x, y, 0) on [-1, 1]². Vertices are the origin, the 4 corners
    # and where the kinks x = 0 (y < 0) and y = 0 (x < 0) cross the boundary
    pwa_pyramid = PWAFunc{Convex,2}([
        Plane([1, 0], 0),
        Plane([0, 1], 0),
        Plane([0, 0], 0),
    ])
    verts = PiecewiseAffineApprox.vertices_by_subsets(
        pwa_pyramid,
        [(-1, 1), (-1, 1)],
    )
    expected = [[0, 0], [-1, -1], [-1, 1], [1, -1], [1, 1], [0, -1], [-1, 0]]
    @test length(verts) == length(expected)
    for p ∈ expected
        @test any(v -> v[1] ≈ p, verts)
    end
    for (x, z) ∈ verts
        @test z ≈ evaluate(pwa_pyramid, x)
    end

    # Invalid bounds
    @test_throws ArgumentError PiecewiseAffineApprox.vertices_by_subsets(
        pwa,
        [(1, -1)],
    )
    @test_throws ArgumentError PiecewiseAffineApprox.vertices_by_subsets(
        pwa,
        [(-1, 1), (-1, 1)],
    )
    @test_throws ArgumentError PiecewiseAffineApprox.vertices_by_subsets(
        pwa,
        [(-Inf, 1)],
    )
end
