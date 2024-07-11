# Approximate a function
Piecewise affine approximations can be obtained either for a given function 𝒻(x) or for a set of sampled points 𝒳 with their corresponding function values 𝒵 = {f(x): ∀x ∈ 𝒳}. Both methods have their applicability depending on the context and the specific application requirements. Methods that use the functions directly can be updated iteratively in real-time applications whereas approximations that use a collection of points may be more relevant when a dataset is provided for off-line construction of the approximate functions prior to optimization.


For piecewise-affine approximations obtained from a pre-defined function, PWA provides a function called approx() which receives as the first arguments the sampling function and the region in the domain for the sampling. The following code creates a convex piecewise affine approximation to 𝓍² on the interval [-1,1] with 5 planes, or segments for this unidimensional case. Notice that the domain [-1,1] will be sampled uniformly thereby splitting the feasible region into 5 segments with corresponding 6 points.
 
```jldoctest

using JuMP, PiecewiseAffineApprox, HiGHS
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)

m = Model(optimizer)
@variable(m, x)
pwa = approx(x -> sum(x.^2), [(-1, 1)], Convex(), MILP(optimizer = optimizer, planes=5))
z = pwaffine(m, x, pwa)

typeof(z)

# output
VariableRef (alias for GenericVariableRef{Float64})
```