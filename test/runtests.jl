using JuMP
using PiecewiseLinearApprox
using Test
using Xpress

const PWL = PiecewiseLinearApprox
const optimizer = optimizer_with_attributes(Xpress.Optimizer, MOI.Silent()=>true)
quadopt(absgap=nothing) = isnothing(absgap) ? optimizer : optimizer_with_attributes(Xpress.Optimizer, MOI.Silent()=>true, MOI.RawOptimizerAttribute("MIPABSSTOP") => absgap)

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
    @testset "Big-M calculations" begin
         include("test_big_M.jl")
    end
end
