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
include("fitting.jl")
include("convexapprox.jl")
include("milp.jl")
include("vertex_enum.jl")
include("linopt.jl")
include("magnani_boyd.jl")
include("kazda_li.jl")
include("interpolation.jl")
include("dcapprox.jl")
include("struct_types.jl")

export Convex, Concave, DC
export MILP, Cluster, Interpol, Alternating, Select, Restart
export Progressive, FullOrder
export Plane
export FunctionEvaluations
export PWAFunc, DCFunc
export approx
export pwaffine
export PWAFormulation, PsiFormulation, LambdaFormulation, DeltaFormulation, PhiFormulation
export Ψ_Formulation, λ_Formulation, Δ_Formulation, Φ_Formulation
export Direction, Graph, Epigraph, Hypograph
export domain
export evaluate
export enforce_curvature
export approximation_error
export dc_decomposition

end # module
