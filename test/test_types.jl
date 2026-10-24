@testset "PWA printing" begin
    pwa = PWAFunc{Convex,2}([Plane([1, 2], 1), Plane([-3, 1], -2)])

    @test summary(pwa) == "PWAFunc{Convex, 2, Int64} with 2 planes"
    @test repr(MIME("text/plain"), pwa) ==
          "PWAFunc{Convex, 2, Int64} with 2 planes:\n z ≥ 1 x₁ + 2 x₂ + 1\n z ≥ -3 x₁ + 1 x₂ + -2\n"

    pwa = PWAFunc{Concave,2}([Plane([1, 2], 1), Plane([-3, 1], -2)])
    @test summary(pwa) == "PWAFunc{Concave, 2, Int64} with 2 planes"
    @test repr(MIME("text/plain"), pwa) ==
          "PWAFunc{Concave, 2, Int64} with 2 planes:\n z ≤ -1 x₁ + -2 x₂ + -1\n z ≤ 3 x₁ + -1 x₂ + 2\n"
end

@testset "Concrete plane types" begin
    @test isconcretetype(Plane{2,Float64})
    @test Plane([1, 2], 1) isa Plane{2,Int}
    @test Plane((1, 2.5), 1) isa Plane{2,Float64}
    @test Plane(2, 1.0) isa Plane{1,Float64}
    @test Plane{2}((1, 2), 3.0) isa Plane{2,Float64}

    # Mixed number types are promoted
    pwa = PWAFunc{Convex,2}([Plane([1, 2], 1), Plane([-3.0, 1.0], -2.0)])
    @test pwa isa PWAFunc{Convex,2,Float64}
    @test PWAFunc{Concave,3}() isa PWAFunc{Concave,3,Float64}
    @test PWAFunc([Plane([1, 2], 1)], Concave()) isa PWAFunc{Concave,2,Int}

    # Evaluation is type stable
    @test (@inferred evaluate(pwa, (0.5, 1.0))) == 3.5
    cv = PWAFunc{Concave,2}(pwa.planes)
    @test (@inferred evaluate(cv, (0.5, 1.0))) == -3.5

    # Adding planes of another number type converts them
    PiecewiseAffineApprox._addplane!(pwa, [1, 1], 0)
    @test pwa.planes[end] === Plane((1.0, 1.0), 0.0)
end

@testset "FunctionEvaluations construction" begin
    f = FunctionEvaluations([(1.0,), (2.0,)], [1.0, 4.0])
    @test FunctionEvaluations([1.0, 2.0], [1.0, 4.0]) == f
    @test FunctionEvaluations(1:2, [1, 4]) isa
          FunctionEvaluations{1,Tuple{Int},Int}
    @test FunctionEvaluations(1:2, (1:2) .^ 2).values == [1, 4]
    # Points as rows (default) or columns of a matrix
    X = [1.0 2.0; 3.0 4.0; 5.0 6.0]
    z = [1.0, 2.0, 3.0]
    @test FunctionEvaluations(X, z).points ==
          [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)]
    @test FunctionEvaluations(X', z; dims = 2).points ==
          [(1.0, 2.0), (3.0, 4.0), (5.0, 6.0)]
    @test FunctionEvaluations(X, z) isa
          FunctionEvaluations{2,Tuple{Float64,Float64},Float64}
    @test_throws ArgumentError FunctionEvaluations(X, z; dims = 3)
    @test_throws ArgumentError FunctionEvaluations(X, [1.0, 2.0])
    @test_throws ArgumentError FunctionEvaluations([1.0, 2.0], [1.0])
    @test_throws ArgumentError FunctionEvaluations([(1.0,), (2.0,)], [1.0])
end

@testset "Function conventions" begin
    # The function is called with the point as a tuple if it accepts one, and
    # otherwise with the coordinates as separate arguments
    s = PiecewiseAffineApprox._sample_uniform
    @test s(sin, [(0, 1)], 5).values ≈ s(x -> sin(x[1]), [(0, 1),], 5).values
    @test s(sin, [(0, 1)], 5).values ≈ sin.(range(0, 1; length = 5))
    g(x, y) = x^2 + 2y
    @test s(g, [(0, 1), (0, 1)], 3).values ≈
          s(x -> x[1]^2 + 2x[2], [(0, 1), (0, 1)], 3).values
    @test s(g, [(0, 1), (0, 1)], 3).values ==
          [g(p...) for p ∈ s(g, [(0, 1), (0, 1)], 3).points]
    pwa = approx(sin, [(0, 2π)], DC(), Interpol(planes = 3); nsamples = 20)
    @test pwa isa DCFunc{1,Float64}
end
