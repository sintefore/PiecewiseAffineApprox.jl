module PiecewiseAffineApprox

using Combinatorics
using Distributions
using JuMP
using LinearAlgebra
using Printf
using Statistics

include("types.jl")
include("convexapprox.jl")
include("linopt.jl")
include("magnani_boyd.jl")
include("interpolation.jl")

function plotconv2D end

export Convex, Concave
export Optimized, Heuristic, Interpol
export Plane
export FunctionEvaluations
export PWLFunc
export approx
export pwlinear
export evaluate
export plotconv2D

end # module
