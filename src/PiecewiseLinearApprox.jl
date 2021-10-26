module PiecewiseLinearApprox

using JuMP
using Printf

using Requires


include("types.jl")
include("convexapprox.jl")
include("linopt.jl")

function __init__()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("plotting.jl")
end

export ConvexPWLFunction
export ConcavePWLFunction
export convex_linearization
export concave_linearization
export convex_pwlinear
export concave_pwlinear
export plot


end # module
