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


## Formulations

PiecewiseAffineApprox.jl provides four different linear programming formulations for adding convex (or concave) piecewise affine functions to optimization models. All formulations are pure linear programs (no binary variables required) and differ in their computational structure and preprocessing requirements. Within the bounds of the input variables they model the same function. The formulations and their names follow [Coffrin et al. (2021)](https://doi.org/10.1016/j.epsr.2021.107191) (preprint: [arXiv:2005.14087](https://arxiv.org/abs/2005.14087)), which considers one-dimensional piecewise linear cost functions. The λ-formulation for functions of two or more variables is a natural extension of the one-dimensional formulation.

Note that all formulations model the epigraph `z ≥ f(x)` of a convex function (or the hypograph `z ≤ f(x)` of a concave function). The value of `z` is therefore only equal to `f(x)` when the objective pushes `z` towards the function, i.e. when minimizing a convex function or maximizing a concave function.

The formulations are selected with the keyword argument `formulation` to `pwaffine`, using the types `PsiFormulation`, `LambdaFormulation`, `DeltaFormulation` and `PhiFormulation`. The Greek names `Ψ_Formulation`, `λ_Formulation`, `Δ_Formulation` and `Φ_Formulation` are available as aliases.

**Ψ-Formulation (`PsiFormulation`, default):** Uses the plane representation directly without requiring vertex enumeration, with one constraint for each plane. This formulation is the most intuitive, works for any dimension and gives the smallest model. Bounds on the input variables are optional. If bounds are given explicitly, the output `z` is also bounded by the smallest and largest value of the function over the bounds, as in Coffrin et al. This requires vertex enumeration.

**λ-Formulation (`LambdaFormulation`):** Represents the piecewise affine function as a convex combination of its vertices. This formulation works for any dimension. It requires vertex enumeration as a preprocessing step, which computes all vertices of the piecewise affine function's epigraph (for convex functions) or hypograph (for concave functions) restricted to the bounds of the input variables.

**Δ-Formulation (`DeltaFormulation`):** A specialized formulation that exploits the ordered structure of breakpoints in one-dimensional piecewise affine functions. Each segment between two consecutive breakpoints is a "bin" with a variable for how much of the segment is used, i.e. the formulation integrates the slope of the function along the x-axis. It requires breakpoint enumeration as a preprocessing step. This formulation only works for 1D problems and will error for higher dimensions.

**Φ-Formulation (`PhiFormulation`):** An alternative specialized formulation for one-dimensional problems. It extends the first segment over the whole domain and adds the change in slope at each later breakpoint, and is the mathematical dual of the Δ-formulation. It is uncommon in practice, and is included mainly for completeness. Like the Δ-formulation, it requires breakpoint enumeration and only works for 1D problems.

The λ, Δ and Φ formulations need a bounded domain for the input variables. The bounds are taken from the JuMP variables if they have finite bounds, or can be given explicitly with the keyword argument `bounds = [(lo₁, hi₁), …, (lo_D, hi_D)]`. Typically, these should be the bounds of the domain used to compute the approximation. For the one-dimensional formulations, adjacent segments with nearly equal slopes are merged to avoid redundant variables.

```jldoctest; output = false
using JuMP, PiecewiseAffineApprox, HiGHS

optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)
pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), MILP(optimizer = optimizer, planes=5))

function solve(formulation)
    m = Model(optimizer)
    @variable(m, -1 <= x <= 1)
    z = pwaffine(m, x, pwa; formulation = formulation)
    @objective(m, Min, z)
    @constraint(m, x >= 0.5)
    optimize!(m)
    return round(value(z); digits=4)
end

[solve(f) for f in (LambdaFormulation(), PsiFormulation(), DeltaFormulation(), PhiFormulation())]

# output
4-element Vector{Float64}:
 0.2653
 0.2653
 0.2653
 0.2653
```

The bounds can also be given explicitly, e.g. when the input is an expression:

```julia
z = pwaffine(m, x + y, pwa; formulation = DeltaFormulation(), bounds = [(-1, 1)])
```

**Choosing a formulation:**
- **Linear programs** (e.g. HiGHS): The choice of formulation has limited impact. In the experiments of Coffrin et al. (DC optimal power flow), λ and Δ were fastest on average, but all formulations were within a factor of 1.7 of the best.
- **Nonlinear models with interior point solvers** (e.g. Ipopt, Knitro): Avoid the Ψ-formulation. Coffrin et al. found it on average 5 times, and up to 10 times, slower than the other formulations for AC optimal power flow. The λ and Δ formulations performed best, with λ best for convex (second order cone) relaxations.
- For **1D problems**: All four formulations are available.
- For **2D and higher**: Only λ-Formulation and Ψ-Formulation are available (Δ and Φ will error).
- **Preprocessing considerations**: λ, Δ, and Φ formulations require enumeration of vertices or breakpoints, which may be computationally expensive for complex approximations. The Ψ-formulation works directly with planes and avoids this preprocessing step.
- **Default**: Ψ-Formulation is used when no formulation is specified. It requires no preprocessing and no bounds, and works well for linear and mixed-integer linear programs.

## Non-convex functions

A difference of convex approximation, see [`DC`](@ref), is added to a model in
the same way. The function is modelled as the sum of its convex and concave
parts, and the part that is bounded against the direction of its curvature
needs one binary variable per plane to select the active plane. The binary
variables select the active plane with big-M constraints computed from the
bounds on the variables, so the variables must be bounded, or the bounds given
with the `bounds` keyword. With `selection = :indicator`, indicator
constraints are used instead for solvers that support them.

For functions of one variable, the Δ- and Φ-formulations are also available.
They model the function by its breakpoints over the bounds, with one binary
variable at each kink of the wrong curvature for the direction, i.e. at the
concave kinks when `z` is minimized. These formulations are ideal and need no
big-M values, so they are preferable to the Ψ-formulation in one dimension.

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
# Only minimized, so binary variables are needed at the concave kinks only
z = pwaffine(m, x, dc; direction = :epigraph, formulation = DeltaFormulation())
@constraint(m, x ≥ 4)
@objective(m, Min, z)
optimize!(m)
round(value(z); digits = 2)

# output
-1.0
```
