@testset "PWA printing" begin
    pwa = PWAFunc{Convex,2}([Plane([1, 2], 1), Plane([-3, 1], -2)])

    @test summary(pwa) == "PWAFunc{Convex, 2} with 2 planes"
    @test repr(MIME("text/plain"), pwa) ==
          "PWAFunc{Convex, 2} with 2 planes:\n z ≥ 1 x₁ + 2 x₂ + 1\n z ≥ -3 x₁ + 1 x₂ + -2\n"

    pwa = PWAFunc{Concave,2}([Plane([1, 2], 1), Plane([-3, 1], -2)])
    @test summary(pwa) == "PWAFunc{Concave, 2} with 2 planes"
    @test repr(MIME("text/plain"), pwa) ==
          "PWAFunc{Concave, 2} with 2 planes:\n z ≤ -1 x₁ + -2 x₂ + -1\n z ≤ 3 x₁ + -1 x₂ + 2\n"
end
