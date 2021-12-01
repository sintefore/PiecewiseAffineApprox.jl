module PiecewiseLinearApprox

using JuMP
using Printf

using Requires
using LinearAlgebra

include("types.jl")
include("convexapprox.jl")
include("linopt.jl")

function __init__()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plotting.jl")
    @require GLMakie="e9467ef8-e4e7-5192-8a1a-b1aee30e663a" include("plotting_makie.jl")
end

# New exports
# Types:
export Convex, Concave
export Optimized, Heuristic
export Plane
export FunctionEvaluations
export PWLFunc

# Methods:
export approx

# Old exports (to be revised)
export ConvexPWLFunction
export ConcavePWLFunction
export convex_linearization
export concave_linearization
export convex_pwlinear
export concave_pwlinear

export plot
export plotconvND

end # module
