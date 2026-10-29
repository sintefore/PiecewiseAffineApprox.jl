```@meta
CurrentModule = PiecewiseAffineApprox
```

# Add a formulation

A formulation method for adding a piecewise affine function to a JuMP model is a
subtype of [`PWAMethod`](@ref), and a method of [`formulate!`](@ref) that adds
its variables and constraints. Everything else is provided by
[`add_pwaformulation`](@ref): the inputs are checked, the output variable is
created and bounded, the variables and constraints are named and recorded, and
the model is left unchanged if adding the function fails.

As an example, the following method models a convex function by one container
of constraints for its planes, `z ≥ αₖ⋅x + βₖ`:

```@example formulation
using JuMP, HiGHS, PiecewiseAffineApprox
using PiecewiseAffineApprox: Epigraph, record!

struct PlanesFormulation <: PWAMethod end

function PiecewiseAffineApprox.formulate!(
    m, x, z, f::PWAFunc{Convex,D}, ::PlanesFormulation, ::Epigraph;
    name, kwargs...,
) where {D}
    P = f.planes
    record!(m, @constraint(
        m,
        [k ∈ eachindex(P)],
        z ≥ sum(P[k].α[i] * x[i] for i ∈ 1:D) + P[k].β,
        base_name = "planes_$(name)",
    ))
    return nothing
end
nothing # hide
```

The method is then used like the built-in ones:

```@example formulation
f = approx(x -> x[1]^2, [(-1, 1)], Convex(), Interpol(planes = 4))
m = Model(HiGHS.Optimizer)
set_silent(m)
@variable(m, -1 ≤ x ≤ 1)
z, formulation = add_pwaformulation(m, x, f; formulation = PlanesFormulation())
@objective(m, Min, z - x / 2)
optimize!(m)
value(x), value(z), formulation
```

It can also be stored with the function, `Formulated(f; formulation =
PlanesFormulation())`, see [`Formulated`](@ref). It is never chosen by
default.

## What a method receives

The method dispatches on the function, the formulation and the direction, and
implements the combinations it supports. For any other combination,
`add_pwaformulation` throws an `ArgumentError`.

- A convex or concave [`PWAFunc`](@ref) is passed in its natural direction,
  `Epigraph()` for a convex and `Hypograph()` for a concave function, which
  needs no binary variables.
- A [`DCFunc`](@ref) is passed in any direction. A convex or concave function
  in another direction is passed as a `DCFunc` with an empty part, so a method
  for `DCFunc` also covers that case.
- The inputs `x` are a vector, and the output `z` is a variable or an
  expression.
- The keyword `bounds` is the box given by the user, or `nothing`. A method
  that needs a box can fall back on the [`domain`](@ref) of the function, or
  on the bounds of the variables.
- The keyword `name` is a label for the names of the variables and
  constraints. Other keywords, such as `selection`, are accepted with
  `kwargs...`.

## What a method must do

Every variable and constraint is added through
[`record!`](@ref PiecewiseAffineApprox.record!), so that it belongs to the
[`PWAFormulation`](@ref) and is removed by `JuMP.delete(m, formulation)`.
Variable bounds and integrality need no recording, as they are removed with
the variables. If the number of variables or constraints in the model changes
by more than was recorded, `add_pwaformulation` warns.

The method returns the box `[(lo₁, hi₁), …]` that it restricts `x` to, or
`nothing` if it does not restrict `x`. The output is bounded over that box, or
over the bounds of the variables if `nothing` is returned. In the example, `z`
gets a lower bound from the bounds of `x`.
