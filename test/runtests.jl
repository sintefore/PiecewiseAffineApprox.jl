using Cbc
using JuMP
using PiecewiseLinearApprox
using Test
using Xpress

const Opt = optimizer_with_attributes(Cbc.Optimizer, MOI.Silent() => true)
const QuadOpt = Xpress.Optimizer

@testset "Test piece-wise linear approximation" begin
    @testset "Quad Approx" begin
        include("testquadapprox.jl")
    end
    @testset "Convex Interpolation" begin
        include("convex_interpolation.jl")
    end
    @testset "1D Convex Approx" begin
        include("testconvexapprox.jl")
    end
    @testset "2D Convex Approx" begin
        include("test_2D_convexapprox.jl")
    end
end
