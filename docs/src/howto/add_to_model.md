# Add to a JuMP model

An approximation is used like a function in a JuMP model. Called with numbers
it is evaluated, and called with JuMP variables or expressions it adds its
constraints to their model and returns the output, which can be used in
expressions, constraints and objectives:

```jldoctest model; output = false
using JuMP, PiecewiseAffineApprox, HiGHS
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)

cost = approx(x -> x[1]^2 + 1 / x[1], [(0.2, 2)], Convex(), Interpol(planes = 5))

m = Model(optimizer)
@variable(m, 0.2 ≤ x[1:3] ≤ 2)
@constraint(m, sum(x) ≥ 3)
@objective(m, Min, sum(cost(x[i]) for i ∈ 1:3))
optimize!(m)
round(objective_value(m); digits = 1), round(cost(1.0); digits = 1)

# output

(6.1, 2.0)
```

Each call adds a new formulation, so an approximation of the same inputs that
is used several times should be added once, e.g. with
`@expression(m, c, cost(x[1]))`.

## Inputs and output

A function of several variables takes the inputs as separate arguments, or as
one vector or tuple. The output is a new variable, or an existing variable or
expression given with the keyword `z`:

```jldoctest model; output = false
using StableRNGs
paraboloid = approx((x, y) -> x^2 + y^2, [(-1, 1), (-1, 1)], Convex(),
    Cluster(planes = 8, rng = StableRNG(1)))

m = Model(optimizer)
@variable(m, -1 ≤ y[1:2] ≤ 1)
@variable(m, w)
paraboloid(y; z = w)
@constraint(m, y[1] + y[2] ≥ 1)
@objective(m, Min, w)
optimize!(m)
round(value(w); digits = 2)

# output

0.47
```

## Bounds of the inputs and the output

The approximation is only valid on the box it was fitted on, which it records
as its [`domain`](@ref). The formulations that need a bounded input take the
box from the keyword `bounds = [(lo₁, hi₁), …]`, the domain, or the bounds of
the variables, in that order. Give `bounds` explicitly when the input is an
expression:

```julia
z = cost(x[1] + x[2]; bounds = [(0.2, 2)])
```

Over the box, the output is bounded by the smallest and largest value of the
function, except in the direction it is pushed: no upper bound for a convex
function that is minimized, and no lower bound for a concave function that is
maximized.

## Record and rollback

[`add_pwaformulation`](@ref) does the same as calling the function, and also
returns a [`PWAFormulation`](@ref), a record of the variables and constraints
that were added, which `delete` removes again. If a function cannot be added,
the model is left unchanged.

```jldoctest model
m = Model()
@variable(m, 0.2 ≤ u ≤ 2)
z, record = add_pwaformulation(m, u, cost; name = "cost")
n = num_variables(m)
delete(m, record)
n, num_variables(m)

# output

(2, 1)
```

## Direction and formulation

A convex function is modelled by its epigraph, `z ≥ f(x)`, and a concave
function by its hypograph, `z ≤ f(x)`, as pure linear programs. The output
then only equals the function when the objective pushes it towards the
function: when a convex function is minimized or a concave function is
maximized. The keyword `direction` changes this:

| Direction | Relation | Use when | Binary variables for |
|:--|:--|:--|:--|
| `Epigraph()` | `z ≥ f(x)` | `z` is only minimized | the concave part |
| `Hypograph()` | `z ≤ f(x)` | `z` is only maximized | the convex part |
| `Graph()` | `z = f(x)` | anything else | both parts |

The directions are not exported, and are given as e.g.
`PiecewiseAffineApprox.Graph()`. The default is the direction that needs no
binary variables for a convex or concave function, and `Graph()` for a
[`DCFunc`](@ref).

The keyword `formulation` selects how the function is modelled: by its planes
([`PsiFormulation`](@ref), the default for convex and concave functions), its
vertices ([`LambdaFormulation`](@ref)), its segments
([`DeltaFormulation`](@ref), the default in one variable when binary variables
are needed) or its slope changes ([`PhiFormulation`](@ref)). All model the same
function within the bounds; see [Formulations](@ref) for when each is best.

```jldoctest model
pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), Interpol(planes = 4))

function solve(formulation)
    m = Model(optimizer)
    @variable(m, -1 ≤ x ≤ 1)
    z = pwa(x; formulation)
    @constraint(m, x ≥ 0.5)
    @objective(m, Min, z)
    optimize!(m)
    return round(value(z); digits = 4)
end

[solve(F) for F ∈ (PsiFormulation(), LambdaFormulation(), DeltaFormulation(), PhiFormulation())]

# output

4-element Vector{Float64}:
 0.2727
 0.2727
 0.2727
 0.2727
```

## Store the options with the function

[`Formulated`](@ref) keeps the options of `add_pwaformulation` with the
function, so that every call uses them:

```julia
cost = Formulated(approx(sin, [(0, 2π)], DC(), Interpol(planes = 4));
    formulation = DeltaFormulation(), direction = PiecewiseAffineApprox.Epigraph())
@objective(m, Min, cost(x))
```

For difference of convex approximations in a model, see
[Non-convex functions](@ref).
