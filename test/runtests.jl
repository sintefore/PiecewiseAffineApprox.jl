using Cbc
using JuMP
using PiecewiseLinearApprox
using Test

@testset "Test piece-wise linear approximation" begin
    include("testquadapprox.jl")
    include("convex_interpolation.jl")
    include("convex_approx.jl")
end

