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

PiecewiseAffineApprox.jl provides four different linear programming formulations for adding convex (or concave) piecewise affine functions to optimization models. All formulations are pure linear programs (no binary variables required) and produce mathematically equivalent results but differ in their computational structure and preprocessing requirements. For more details see [Coffrin et al. (2021)](https://doi.org/10.1016/j.epsr.2021.107191).

**Ψ-Formulation (default):** Uses the plane representation directly without requiring vertex enumeration. The formulation directly models the piecewise affine function as a pointwise maximum of the planes. This formulation works for any dimension and can be more efficient when vertex enumeration is expensive (large number of planes).

**λ-Formulation:** Represents the piecewise affine function as a convex combination of its vertices. This formulation works for any dimension and is the default choice. It requires vertex enumeration as a preprocessing step, which computes all vertices of the piecewise affine function's epigraph (for convex functions) or hypograph (for concave functions).

**Δ-Formulation:** A specialized formulation that exploits the ordered structure of breakpoints in one-dimensional piecewise affine functions. It requires breakpoint enumeration as a preprocessing step. This formulation only works for 1D problems and will error for higher dimensions.

**Φ-Formulation:** An alternative specialized formulation for one-dimensional problems. Like the Δ-formulation, it requires breakpoint enumeration and only works for 1D problems.

```jldoctest; output = false
using JuMP, PiecewiseAffineApprox, HiGHS

optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)
m = Model(optimizer)
@variable(m, x)

pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), MILP(optimizer = optimizer, planes=5))

z1 = pwaffine(m, x, pwa; formulation = λ_Formulation())
z2 = pwaffine(m, x, pwa; formulation = Ψ_Formulation())
z3 = pwaffine(m, x, pwa; formulation = Δ_Formulation())
z4 = pwaffine(m, x, pwa; formulation = Φ_Formulation())

@objective(m, Min, z1)
@constraint(m, x >= 0.5)
optimize!(m)
round(value(z1); digits=4)

# output
0.2653
```

**Choosing a formulation:**
- For **1D problems**: All four formulations are available. The Δ and Φ formulations may offer computational advantages for problems with many breakpoints.
- For **2D and higher**: Only λ-Formulation and Ψ-Formulation are available (Δ and Φ will error).
- **Preprocessing considerations**: λ, Δ, and Φ formulations require enumeration of vertices or breakpoints, which may be computationally expensive for complex approximations. The Ψ-formulation works directly with planes and avoids this preprocessing step.
- **Default**: λ-Formulation is used when no formulation is specified and works well in most cases.

