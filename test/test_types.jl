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
