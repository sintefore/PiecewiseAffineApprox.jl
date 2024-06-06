module PiecewiseAffineApprox

using Combinatorics
using Distributions
using JuMP
using LinearAlgebra
using Printf
using Statistics
using Base.Threads: @spawn

include("types.jl")
include("convexapprox.jl")
include("linopt.jl")
include("magnani_boyd.jl")
include("kazda_li.jl")
include("interpolation.jl")

export Convex, Concave
export Optimized, Heuristic, Interpol
export ProgressiveFitting, FullOrderFitting
export Plane
export FunctionEvaluations
export PWAFunc
export approx
export pwaffine
export evaluate
export convexify

end # module
