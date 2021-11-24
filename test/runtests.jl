using Cbc
using Xpress
using JuMP
using PiecewiseLinearApprox
using Test

@testset "Test piece-wise linear approximation" begin
    include("testquadapprox.jl")
    include("convex_interpolation.jl")
    include("testconvexapprox.jl")
    include("test_2D_convexapprox.jl")
end
