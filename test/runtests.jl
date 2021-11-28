using JuMP
using PiecewiseLinearApprox
using SCIP
using Test

const optimizer = optimizer_with_attributes(SCIP.Optimizer, MOI.Silent()=>true)
quadopt(absgap=nothing) = isnothing(absgap) ? optimizer : optimizer_with_attributes(SCIP.Optimizer, MOI.Silent()=>true, MOI.RawOptimizerAttribute("limits/absgap") => absgap)

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
