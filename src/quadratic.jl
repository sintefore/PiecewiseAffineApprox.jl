# Piecewise quadratic functions in a JuMP model. A convex PWQFunc is added in
# the epigraph direction and a concave one in the hypograph direction, with
# one convex quadratic constraint for each piece, as the Ψ-formulation does
# for planes. The other directions would need binary variables and a
# nonconvex model, and are not available.

"""
    add_pwaformulation(m::JuMP.Model, x, f::PWQFunc{C,D}; z, formulation, bounds, direction, name)

Add the piecewise quadratic function `f` of the inputs `x` to the model `m`
with one constraint for each piece: `z ≥ xᵀQx + αᵀx + β` for a convex function
and `z ≤ -(xᵀQx + αᵀx + β)` for a concave one (with the stored pieces, see
[`PWQFunc`](@ref)). Returns `(z, formulation)` as for a [`PWAFunc`](@ref), and
can likewise be called as `f(x)`.

Each constraint is a rotated second order cone, `u ≥ ‖Lx‖²` with `Q = LᵀL`
and `u` the output minus the affine part, or a linear constraint for a piece
without curvature. Conic solvers such as Clarabel, Mosek or SCS take these
directly, and JuMP reformulates them as convex quadratic constraints for
solvers such as Gurobi or Ipopt. As a quadratic constraint, a piece whose
curvature is singular (zero in some direction) would not be accepted by the
reformulation for conic solvers. As for the Ψ-formulation of a `PWAFunc`, the
relation is only tight if the objective pushes `z` towards the function.

# Keyword arguments
- `z`: the output, a JuMP variable or expression, or `nothing` (default) to
  create one. The output is not bounded.
- `formulation`: `nothing` or `PsiFormulation()`, the only one available.
- `direction`: `Epigraph()` for a convex and `Hypograph()` for a concave
  function, the default and the only one available.
- `bounds`: a box `[(lo₁, hi₁), …]` that is added as constraints on `x`, or
  `nothing` (default).
- `name`: the label in the names of the variables and constraints added.
"""
function add_pwaformulation(m::JuMP.Model, x, f::PWQFunc; kwargs...)
    return _add_pwaformulation(m, x, f; kwargs...)
end

function _pwaffine!(
    m::JuMP.Model,
    x,
    f::PWQFunc{C,D},
    name;
    z = nothing,
    formulation::Union{Nothing,PWAMethod} = nothing,
    bounds = nothing,
    direction::Direction = _natural_direction(C),
    selection = :bigm,
) where {C,D}
    xs = _inputs(x)
    length(xs) == D || throw(
        ArgumentError("The function has $D input(s), but x has $(length(xs)) element(s)"),
    )
    isnothing(formulation) || formulation isa PsiFormulation || throw(
        ArgumentError(
            "$(nameof(typeof(formulation))) is not available for a PWQFunc, " *
            "which is added with one quadratic constraint for each piece",
        ),
    )
    direction == _natural_direction(C) || throw(
        ArgumentError(
            "A PWQFunc{$C} can only be added in the direction " *
            "$(nameof(typeof(_natural_direction(C))))(), which needs no binary " *
            "variables",
        ),
    )
    if isnothing(z)
        z = record!(m, JuMP.@variable(m, base_name = "z_$(name)"))
    end
    for (k, q) ∈ enumerate(f.pieces)
        # The slack u = z - (αᵀx + β) for a convex function, or -z - (αᵀx + β)
        # for a concave one, must be at least xᵀQx = ‖Lx‖²
        affine = dot(q.α, xs) + q.β
        u = C === Convex ? z - affine : -z - affine
        L = _factor(q.Q)
        con = if isempty(L)
            JuMP.@constraint(m, u ≥ 0)
        else
            # 2 u (1/2) ≥ ‖Lx‖²
            JuMP.@constraint(m, [u; 0.5; L * xs] ∈ JuMP.RotatedSecondOrderCone())
        end
        record!(m, con)
        JuMP.set_name(con, "pwq_$(name)_$(k)")
    end
    if !isnothing(bounds)
        lo, hi = _check_bounds(bounds, D)
        record!(m, JuMP.@constraint(
            m,
            [i = 1:D],
            lo[i] ≤ xs[i] ≤ hi[i],
            base_name = "pwq_$(name)_bounds"
        ))
    end
    return z
end

# A factor L with Q = LᵀL for a positive semidefinite Q, with one row for each
# eigenvalue that is not zero up to rounding
function _factor(Q::AbstractMatrix)
    isempty(Q) && return zeros(0, 0)
    λ, V = eigen(Symmetric(Matrix{Float64}(Q)))
    keep = λ .> 1e-12 * max(maximum(abs, λ), floatmin(Float64))
    return Diagonal(sqrt.(λ[keep])) * V[:, keep]'
end
