module PiecewiseAffineApprox

using JuMP
using Printf
using Requires
using LinearAlgebra
using Statistics
using Distributions
using Combinatorics

include("types.jl")
include("convexapprox.jl")
include("linopt.jl")
include("magnani_boyd.jl")
include("interpolation.jl")

function __init__()
    @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include(
        "plotting.jl",
    )
    @require GLMakie = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a" include(
        "plotting_makie.jl",
    )
end

export Convex, Concave
export Optimized, Heuristic, Interpol
export Plane
export FunctionEvaluations
export PWLFunc
export approx
export pwlinear
export evaluate



end # module
