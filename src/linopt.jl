
mutable struct PWLData
    counter::Int
    PWLData() = new(0)
end

function initPWL!(m::JuMP.Model)
    if !haskey(m.ext, :PWL)
        m.ext[:PWL] = PWLData()
    end
    return nothing
end

const VarOrAff = Union{JuMP.VariableRef,JuMP.AffExpr}

constr(::Type{Convex}, m, z, p, x) = JuMP.@constraint(m, z ≥ dot(p.α, x) + p.β)
function constr(::Type{Concave}, m, z, p, x)
    JuMP.@constraint(m, z ≤ dot(-1 .* p.α, x) - p.β)
end

"""
    pwaffine(m::JuMP.Model, x::Tuple, pwa::PWAFunc{C,D}; z=nothing, kwargs...) where {C,D}

Add constraints to JuMP-model `m` for JuMP-variable `z` as a  
piecewise linear function/approximation `pwa` of JuMP-variables `x`    
"""
function pwaffine(m::JuMP.Model, x, pwa::PWAFunc{C,D}; z = nothing) where {C,D}
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter

    if isnothing(z)
        # @warn "NB: Skipping bounds for now"
        z = JuMP.@variable(m, base_name = "z_$(counter)")
    end
    for (k, p) ∈ enumerate(pwa.planes)
        con = constr(C, m, z, p, x)
        JuMP.set_name(con, "pwa_$(counter)_$(k)")
    end
    return z
end

function pwaffine(
    m::JuMP.Model,
    x,
    fevals::FunctionEvaluations,
    curvature::Curvature,
    a::Algorithm;
    z = nothing,
    kwargs...,
)
    return pwaffine(m, x, approx(fevals, curvature, a); z = z, kwargs...)
end

"""
    pwaffine(m::JuMP.Model, x, f::DCFunc{D}; z = nothing, direction = :graph,
             formulation = :bigm, bounds = nothing)

Add constraints to the JuMP model `m` for the variable `z` as the difference of
convex piecewise linear function `f` of the `D` variables `x`, and return `z`.
The function is modelled as z = z₁ + z₂ with z₁ its convex and z₂ its concave
part. A part that is the maximum (minimum) of planes needs binary variables to
select the active plane for the constraints that bound it from above (below).

- `direction = :graph` (default): z equals the function, with binary variables
  for both parts.
- `direction = :epigraph`: z is at least the function, which is sufficient when
  z is minimized or only bounded from below. The convex part needs no binary
  variables, as for a convex [`PWAFunc`](@ref).
- `direction = :hypograph`: z is at most the function, which is sufficient
  when z is maximized. The concave part needs no binary variables.

The active plane is selected by one binary variable per plane, and the
constraints for the other planes are relaxed with `formulation = :bigm` by
big-M values computed from the bounds on `x`, or with
`formulation = :indicator` by indicator constraints, which do not need
bounds but require solver support. The bounds are the lower and upper bounds
of the variables `x`, or given with `bounds` as a vector of `(lower, upper)`
tuples, e.g. when `x` are expressions.
"""
function pwaffine(
    m::JuMP.Model,
    x,
    f::DCFunc{D};
    z = nothing,
    direction = :graph,
    formulation = :bigm,
    bounds = nothing,
) where {D}
    direction ∈ (:graph, :epigraph, :hypograph) || throw(
        ArgumentError(
            "Unsupported direction :$direction, use :graph, :epigraph or :hypograph",
        ),
    )
    formulation ∈ (:bigm, :indicator) || throw(
        ArgumentError(
            "Unsupported formulation :$formulation, use :bigm or :indicator",
        ),
    )
    length(x) == D || throw(
        ArgumentError("The function has $D inputs, got $(length(x)) variables"),
    )
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter

    if isnothing(z)
        z = JuMP.@variable(m, base_name = "z_$(counter)")
    end
    # The box containing x, needed for the big-M values when a part has
    # binary variables
    exact₁ = direction != :epigraph && _planes(f.convex) > 1
    exact₂ = direction != :hypograph && _planes(f.concave) > 1
    box = if formulation == :bigm && (exact₁ || exact₂)
        _variable_box(x, bounds)
    else
        nothing
    end
    z₁ = _dc_part!(
        m,
        x,
        f.convex,
        direction != :epigraph,
        formulation,
        box,
        "$(counter)_cvx",
    )
    z₂ = _dc_part!(
        m,
        x,
        f.concave,
        direction != :hypograph,
        formulation,
        box,
        "$(counter)_ccv",
    )
    con = JuMP.@constraint(m, z == z₁ + z₂)
    JuMP.set_name(con, "pwa_$(counter)")
    return z
end

# The value of the plane in the sense of the curvature, as an affine expression
_plane_expr(::Type{Convex}, p::Plane, x) = dot(p.α, x) + p.β
_plane_expr(::Type{Concave}, p::Plane, x) = dot(-1 .* p.α, x) - p.β

# Add the variable zᵢ for a convex (concave) part of a DCFunc: zᵢ is at least
# (at most) every plane, and with exact = true also at most (at least) the
# active plane, selected by binary variables. Returns zᵢ, or zero for an
# empty part.
function _dc_part!(
    m::JuMP.Model,
    x,
    pwa::PWAFunc{C,D},
    exact,
    formulation,
    box,
    name,
) where {C,D}
    K = _planes(pwa)
    K == 0 && return zero(JuMP.AffExpr)
    zᵢ = JuMP.@variable(m, base_name = "z_$(name)")
    v = [_plane_expr(C, p, x) for p ∈ pwa.planes]
    # The constraints giving the maximum (minimum) of the planes
    for k ∈ 1:K
        con =
            C == Convex ? JuMP.@constraint(m, zᵢ ≥ v[k]) :
            JuMP.@constraint(m, zᵢ ≤ v[k])
        JuMP.set_name(con, "pwa_$(name)_$(k)")
    end
    exact || return zᵢ
    if K == 1
        con =
            C == Convex ? JuMP.@constraint(m, zᵢ ≤ v[1]) :
            JuMP.@constraint(m, zᵢ ≥ v[1])
        JuMP.set_name(con, "pwa_$(name)_active_1")
        return zᵢ
    end
    u = JuMP.@variable(m, [1:K], Bin, base_name = "u_$(name)")
    con = JuMP.@constraint(m, sum(u) == 1)
    JuMP.set_name(con, "pwa_$(name)_select")
    M = formulation == :bigm ? _dc_big_M(pwa, box) : nothing
    for k ∈ 1:K
        con = if formulation == :indicator && C == Convex
            JuMP.@constraint(m, u[k] --> {zᵢ ≤ v[k]})
        elseif formulation == :indicator
            JuMP.@constraint(m, u[k] --> {zᵢ ≥ v[k]})
        elseif C == Convex
            JuMP.@constraint(m, zᵢ ≤ v[k] + M[k] * (1 - u[k]))
        else
            JuMP.@constraint(m, zᵢ ≥ v[k] - M[k] * (1 - u[k]))
        end
        JuMP.set_name(con, "pwa_$(name)_active_$(k)")
    end
    return zᵢ
end

# The bounds (lower, upper) of each of the variables x, from `bounds` if given
# and otherwise from the bounds of the variables
function _variable_box(x, bounds)
    if !isnothing(bounds)
        length(bounds) == length(x) || throw(
            ArgumentError(
                "bounds must have one (lower, upper) tuple for each variable",
            ),
        )
        box = [(Float64(b[1]), Float64(b[2])) for b ∈ bounds]
    else
        box = map(x) do xᵢ
            if xᵢ isa JuMP.VariableRef && JuMP.is_fixed(xᵢ)
                return (JuMP.fix_value(xᵢ), JuMP.fix_value(xᵢ))
            elseif xᵢ isa JuMP.VariableRef &&
                   JuMP.has_lower_bound(xᵢ) &&
                   JuMP.has_upper_bound(xᵢ)
                return (JuMP.lower_bound(xᵢ), JuMP.upper_bound(xᵢ))
            end
            return throw(
                ArgumentError(
                    "pwaffine with formulation = :bigm requires lower and " *
                    "upper bounds on the variables x, as variable bounds or " *
                    "with the bounds keyword. Use formulation = :indicator " *
                    "for unbounded variables.",
                ),
            )
        end
    end
    all(b -> isfinite(b[1]) && isfinite(b[2]) && b[1] ≤ b[2], box) || throw(
        ArgumentError("The bounds on x must be finite with lower ≤ upper"),
    )
    return collect(box)
end

# Big-M for each plane k of a part: an upper bound on the distance from the
# part to plane k over the box, max over planes i and x in the box of
# (αᵢ - αₖ)⋅x + βᵢ - βₖ. The same expression applies to the concave part,
# where the stored planes are negated.
function _dc_big_M(pwa::PWAFunc{C,D}, box) where {C,D}
    return map(pwa.planes) do pₖ
        return maximum(pwa.planes) do pᵢ
            return pᵢ.β - pₖ.β + sum(
                max(
                    (pᵢ.α[j] - pₖ.α[j]) * box[j][1],
                    (pᵢ.α[j] - pₖ.α[j]) * box[j][2],
                ) for j ∈ 1:D
            )
        end
    end
end
