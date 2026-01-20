module PiecewiseAffineApprox

using Combinatorics
using Distributions
using JuMP
using LinearAlgebra
using Printf
using Statistics
using StructTypes
using Base.Threads: @spawn

include("types.jl")
include("convexapprox.jl")
include("vertex_enum.jl")
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
export λ_Formulation, Ψ_Formulation, Δ_Formulation, Φ_Formulation
export evaluate
export enforce_curvature

end # module
