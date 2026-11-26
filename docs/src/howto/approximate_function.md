# Approximate a function
Piecewise affine approximations can be obtained either for a given function 𝒻(x) or for a set of sampled points 𝒳 with their corresponding function values 𝒵 = {f(x): ∀x ∈ 𝒳}. Both methods have their applicability depending on the context and the specific application requirements. Methods that use the functions directly can be updated iteratively in real-time applications whereas approximations that use a collection of points may be more relevant when a dataset is provided for off-line construction of the approximate functions prior to optimization.


For piecewise-affine approximations obtained from a pre-defined function, PWA provides a function called approx() which receives as the first arguments the sampling function and the region in the domain for the sampling. The following code creates a convex piecewise affine approximation to 𝓍² on the interval [-1,1] with 5 planes, or segments for this unidimensional case. Notice that the domain [-1,1] will be sampled uniformly thereby splitting the feasible region into 5 segments with corresponding 6 points.
 
```jldoctest

using JuMP, PiecewiseAffineApprox, HiGHS
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)

m = Model(optimizer)
@variable(m, x)
pwa = approx(x -> sum(x.^2), [(-1, 1)], Convex(), MILP(optimizer = optimizer, planes=5))
z = pwa(x)

typeof(z)

# output
VariableRef (alias for GenericVariableRef{Float64})
```

## Approximations on one side of the function

With `strict = :outer`, the approximation of a convex function lies below the
sampled values, and with `strict = :inner` above them (interchanged for concave
functions). Between the samples it can still cross the function. With
`everywhere = true` it stays on its side over the whole box, as needed for a
valid relaxation, a bound or a conservative constraint:

```jldoctest
using PiecewiseAffineApprox, HiGHS, JuMP
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)

algorithm = Progressive(; optimizer, tolerance = 0.05, metric = :max, strict = :inner)
r = fit(x -> x[1]^2, [(-1, 1)], Convex(), algorithm; everywhere = true)
b = r.details.bounds
b.below ≤ 0, round(b.above; digits = 3)

# output
(true, 0.054)
```

The function must be convex or concave on the box as stated, which is used to
find where the approximation is furthest on the wrong side, and its planes are
shifted by that amount. Here the approximation is above ``x^2`` everywhere, by
at most 0.054: the tolerance at the samples plus the shift. `bounds` are the
largest deviations over the box, see [`error_bounds`](@ref), which also gives
them for any convex or concave approximation of a function.

## Non-convex functions

A function that is neither convex nor concave is approximated with `DC()` as a
difference of convex functions: the sum of a convex and a concave piecewise
affine function. Any of the algorithms can be used, and fits both parts, with
`planes` planes in the convex part and `concave_planes` in the concave part. The
concave part is the expensive one to model, so it usually gets fewer planes.

```jldoctest
using PiecewiseAffineApprox, HiGHS, JuMP
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)

dc = approx(sin, [(0, 2π)], DC(), Interpol(; optimizer, planes = 4, concave_planes = 2); nsamples = 30)
length(dc.convex.planes), length(dc.concave.planes)

# output
(4, 2)
```

Use `Alternating` to fit the two parts with different algorithms.

A tolerance, as in `Progressive` or `Select` with `criterion = :tolerance`,
applies to the approximation as a whole: it is divided between the two parts
so that their sum is within the tolerance of the data.
