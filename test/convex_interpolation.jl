
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
