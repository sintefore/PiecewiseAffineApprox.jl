module PiecewiseAffineApprox

using Combinatorics
using JuMP
using LinearAlgebra
using Printf
using Random
using StableRNGs
using Statistics
using StructTypes
using Base.Threads: @spawn

include("types.jl")
include("convexapprox.jl")
include("milp.jl")
include("linopt.jl")
include("magnani_boyd.jl")
include("kazda_li.jl")
include("interpolation.jl")
include("dcapprox.jl")
include("struct_types.jl")

export Convex, Concave, DC
export MILP, Cluster, Interpol, Alternating
export Progressive, FullOrder
export Plane
export FunctionEvaluations
export PWAFunc, DCFunc
export approx
export pwaffine
export evaluate
export enforce_curvature
export dc_decomposition

end # module
