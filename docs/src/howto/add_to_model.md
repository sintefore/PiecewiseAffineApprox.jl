# Add approximation to model

 A piecewise-affine approximation ̃𝓏 ≈ 𝒻(𝔁) can be added to a JuMP model for optimization purposes. In that case, the optimization variables 𝔁 will denote the coordinates in the domain whereas the approximate value of function 𝒻 will be represented by ̃𝓏.
 
The code below creates a piecewise affine approximate to 𝓍² on the interval [-1,1] using the JuMP variable 𝓍. The piecewise affine approximation is linked the other variable in the model 𝓏. The objective is set as the minimization of this variable suject to 𝓍 being greater or equal than 0.5. After the optimization of the constrained problem, the optimal value of 𝓏 is 0.2653.

```jldoctest; output = false
using JuMP, PiecewiseAffineApprox, HiGHS


optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)
m = Model(optimizer)
@variable(m, x)
# Create a piecewise linear approximation to x^2 on the interval [-1, 1]
pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), MILP(optimizer = optimizer, planes=5))
# Add the pwa function to the model
z = pwaffine(m, x, pwa)
# Minimize
@objective(m, Min, z)
# Check approximation/solution at x = 0.5
@constraint(m, x >= 0.5)
optimize!(m);
round(value(z); digits=4)

# output
0.2653
```

A piecewise affine approximation can be also be linked to existing variables in an existing JuMP model. In the following code, a 3D piecewise affine approximation is obtained for existing variables of a model (xvar and yvar). The approximation is constructed for a given dataset (X,z) and then linked to a JuMP model, which is optimized in the sequence.


```jldoctest; output=false
using JuMP, PiecewiseAffineApprox, HiGHS
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)

xg = [i for i ∈ -1:0.5:1]
yg = [j for j ∈ -1:0.5:1]
X = [repeat(xg, inner = [size(yg, 1)]) repeat(yg, outer = [size(xg, 1)])]
z = X[:, 1] .^ 2 + X[:, 2] .^ 2

np = 17

m = Model()
@variable(m, xvar)
@variable(m, yvar)
@variable(m, approx_f)

tuple_var = (xvar, yvar)

y = pwaffine(
    m,
    tuple_var,
    approx(
        FunctionEvaluations(tuple.(eachcol(X)...), z),
        Convex(),
        MILP(
            optimizer = optimizer,
            planes = np,
            strict = :outer,
            metric = :l1,
        ),
    );
    z = approx_f,
)

@objective(m, Min, y)
set_optimizer(m, optimizer)
@constraint(m, xvar == √0.5)
@constraint(m, yvar == √0.5)
optimize!(m)
round(value.(y))
# output
1.0
```

## Non-convex functions

A difference of convex approximation, see [`DC`](@ref), is added to a model in
the same way. The function is modelled as the sum of its convex and concave
parts, and the part that is bounded against the direction of its curvature
needs one binary variable per plane to select the active plane. The binary
variables select the active plane with big-M constraints computed from the
bounds on the variables, so the variables must be bounded, or the bounds given
with the `bounds` keyword. With `formulation = :indicator`, indicator
constraints are used instead for solvers that support them.

By default `z` equals the approximation. If `z` is only minimized, pass
`direction = :epigraph`, and if it is only maximized `direction = :hypograph`.
Then only the concave or the convex part, respectively, needs binary variables.

```jldoctest; output = false
using JuMP, PiecewiseAffineApprox, HiGHS
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)

m = Model(optimizer)
@variable(m, 0 ≤ x ≤ 2π)
# A difference of convex approximation of sin(x) on [0, 2π]
dc = approx(sin, [(0, 2π)], DC(), Interpol(; optimizer, planes = 3))
# Only minimized, so binary variables are needed for the concave part only
z = pwaffine(m, (x,), dc; direction = :epigraph)
@constraint(m, x ≥ 4)
@objective(m, Min, z)
optimize!(m)
round(value(z); digits = 2)

# output
-1.0
```
