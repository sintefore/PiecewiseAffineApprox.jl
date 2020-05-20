module PiecewiseLinearApprox

using JuMP
using Memoize
using Plots
using Printf

include("computelinearapprox.jl")

export bestlinearization
export convexlinearization

include("types.jl")
include("convexapprox.jl")
include("linopt.jl")


export ConvexPWLFunction
export ConcavePWLFunction
export convex_linearization
export concave_linearization
export convex_pwlinear
export concave_pwlinear


end # module
