# Add approximation to model

 A piecewise-affine approximation ̃𝓏 ≈ 𝒻(𝔁) can be added to a JuMP model for optimization purposes. In that case, the optimization variables 𝔁 will denote the coordinates in the domain whereas the approximate value of function 𝒻 will be represented by ̃𝓏.
 
The code below creates a piecewise affine approximate to 𝓍² on the interval [-1,1] using the JuMP variable 𝓍. The piecewise affine approximation is linked the other variable in the model 𝓏. The objective is set as the minimization of this variable suject to 𝓍 being greater or equal than 0.5. After the optimization of the constrained problem, the optimal value of 𝓏 is 0.2653.

```julia
using JuMP, PiecewiseAffineApprox, HiGHS

m = Model(HiGHS.Optimizer)
@variable(m, x)
# Create a piecewise linear approximation to x^2 on the interval [-1, 1]
pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), Optimized(optimizer = HiGHS.Optimizer, planes=5))
# Add the pwa function to the model
z = pwaffine(m, x, pwa)
# Minimize
@objective(m, Min, z)
# Check approximation/solution at x = 0.5
@constraint(m, x >= 0.5)
optimize!(m)
value(z) # 0.2653
```

A piecewise affine approximation can be also be linked to existing variables in an existing JuMP model. In the following code, a 3D piecewise affine approximation is obtained for existing variables of a model (xvar and yvar). The approximation is constructed for a given dataset (X,z) and then linked to a JuMP model, which is optimized in the sequence.


```julia
using JuMP, PiecewiseAffineApprox, HiGHS

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
        FunctionEvaluations(mat2tuples(X), z),
        Convex(),
        Optimized(
            optimizer = optimizer,
            planes = np,
            strict = :outer,
            pen = :l1,
        ),
    );
    z = approx_f,
)

@objective(m, Min, y)
set_optimizer(m, optimizer)
@constraint(m, xvar == √0.5)
@constraint(m, yvar == √0.5)
optimize!(m)
```