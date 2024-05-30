
@testset "Convex Interpolation" begin
    # Points to approximate
    x = [i for i ∈ -2:0.1:2]
    z = x .^ 2

    # Test interpolation routine
    pwl = PWA._interpolatepw(x, z, Interpol(optimizer = optimizer, planes = 5))

    @test length(pwl.x) == 6
    @test length(pwl.z) == 6

    @test pwl.z[2] ≈ 1.44
    @test pwl.z[3] ≈ 0.16

    # Test convexification
    cpwl = PWA._convexify1D(pwl, Interpol(optimizer = optimizer))

    @test cpwl.x[2] ≈ -1.2
    @test cpwl.x[3] ≈ -0.4
end
