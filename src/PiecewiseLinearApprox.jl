module PiecewiseLinearApprox

using Cbc
using JuMP
using Memoize

include("computelinearapprox.jl")

export bestlinearization
export convexlinearization


end # module
