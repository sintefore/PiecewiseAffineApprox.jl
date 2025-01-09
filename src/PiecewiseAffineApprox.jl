module PiecewiseAffineApprox

using Combinatorics
using Distributions
using JuMP
using LinearAlgebra
using Printf
using Statistics
using StructTypes 
using JSON3 # TODO: move to extension
using Base.Threads: @spawn

include("types.jl")
include("convexapprox.jl")
include("linopt.jl")
include("magnani_boyd.jl")
include("kazda_li.jl")
include("interpolation.jl")
include("struct_types.jl")

export Convex, Concave
export MILP, Cluster, Interpol
export Progressive, FullOrder
export Plane
export FunctionEvaluations
export PWAFunc
export approx
export pwaffine
export evaluate
export enforce_curvature

end # module
