using Cbc
using JuMP
using PiecewiseLinearApprox
using Test

@testset "Test piece-wise linear approximation" begin
    include("testquadapprox.jl")
end

