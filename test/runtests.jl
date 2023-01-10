using HiGHS
using JuMP
using PiecewiseLinearApprox
using StableRNGs
using Test

rng = StableRNG(123)

const PWL = PiecewiseLinearApprox
const optimizer =
    optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
function quadopt(absgap = nothing)
    return isnothing(absgap) ? optimizer :
           optimizer_with_attributes(
        HiGHS.Optimizer,
        MOI.Silent() => true,
        MOI.RawOptimizerAttribute("MIPABSSTOP") => absgap,
    )
end

@testset "Test piece-wise linear approximation" begin
    # @testset "Quad Approx" begin
    #     include("testquadapprox.jl")
    # end
    @testset "Convex Interpolation" begin
        include("convex_interpolation.jl")
    end
    # @testset "1D Convex Approx" begin
    #     include("testconvexapprox.jl")
    # end
    # @testset "2D Convex Approx" begin
    #     include("test_2D_convexapprox.jl")
    # end
    # @testset "Big-M calculations" begin
    #      include("test_big_M.jl")
    # end
end
