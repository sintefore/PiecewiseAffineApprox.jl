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

"""
    PWAFormulation

Abstract type for the formulations used by [`pwaffine`](@ref) to add a piecewise affine
function to a JuMP model. The formulations follow Coffrin et al. (2021), and are also
available under their Greek names `Ψ_Formulation`, `λ_Formulation`, `Δ_Formulation` and
`Φ_Formulation`.
"""
abstract type PWAFormulation end

"Ψ-formulation: one constraint for each plane of the piecewise affine function."
struct PsiFormulation <: PWAFormulation end
"λ-formulation: convex combination of the vertices of the piecewise affine function."
struct LambdaFormulation <: PWAFormulation end
"Δ-formulation (1D only): a bin for each segment between consecutive breakpoints."
struct DeltaFormulation <: PWAFormulation end
"Φ-formulation (1D only): first segment extended, plus the change in slope at each breakpoint."
struct PhiFormulation <: PWAFormulation end

# Greek aliases, following the notation of Coffrin et al. (2021)
const Ψ_Formulation = PsiFormulation
const λ_Formulation = LambdaFormulation
const Δ_Formulation = DeltaFormulation
const Φ_Formulation = PhiFormulation

const VarOrAff = Union{JuMP.VariableRef,JuMP.AffExpr}

"""
    pwaffine(m::JuMP.Model, x, pwa::PWAFunc{C,D}; z, formulation::PWAFormulation, bounds) where {C,D}

Add constraints to JuMP-model `m` for JuMP-variable `z` as a
piecewise linear function/approximation `pwa` of JuMP-variables `x`.

# Arguments
- `m::JuMP.Model`: The JuMP model to add constraints to
- `x`: JuMP variable(s) representing the input to the piecewise affine function
- `pwa::PWAFunc{C,D}`: The piecewise affine function to model

# Keyword Arguments
- `z`: JuMP variable or expression representing the output. If `nothing`, a new variable is created
- `formulation::PWAFormulation`: The formulation to use. Options are:
  - `PsiFormulation()` (default): Direct plane representation. Works for any dimension. No preprocessing required.
  - `LambdaFormulation()`: Convex combination of vertices. Works for any dimension. Requires vertex enumeration.
  - `DeltaFormulation()`: Specialized for 1D. Uses ordered breakpoint structure. Requires breakpoint enumeration.
  - `PhiFormulation()`: Alternative 1D formulation. Requires breakpoint enumeration.
- `bounds`: Box `[(lo₁, hi₁), …, (lo_D, hi_D)]` restricting the domain of `x`. Required by the
  λ, Δ and Φ formulations unless all of `x` are JuMP variables with finite bounds (or fixed
  values), in which case these are used. For the Ψ formulation the bounds are optional. If
  given, they are added as constraints on `x`, and `z` is bounded by the smallest and
  largest value of the function over the box.

# Returns
- The JuMP variable `z` representing the piecewise affine function value

# Note
All formulations are pure linear programs (no binary variables). They model z ≥ f(x) for
convex and z ≤ f(x) for concave functions, i.e. the relation is only tight if the
objective pushes `z` towards the function (minimizing a convex or maximizing a concave
function). The Δ and Φ formulations only work for 1D problems and will error for higher
dimensions. With interior point solvers for nonlinear models, the λ and Δ formulations are
typically much faster than Ψ (see Coffrin et al. (2021), arXiv:2005.14087).

# Example
```julia
m = Model(HiGHS.Optimizer)
@variable(m, -1 ≤ x ≤ 1)
pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), MILP(optimizer = optimizer, planes=5))
z = pwaffine(m, x, pwa; formulation = LambdaFormulation())
```
"""
function pwaffine(
    m::JuMP.Model,
    x,
    pwa::PWAFunc{C,D};
    z = nothing,
    formulation::PWAFormulation = PsiFormulation(),
    bounds = nothing,
) where {C,D}
    n = length(_inputs(x))
    n == D || throw(
        ArgumentError(
            "The pwa function has dimension $D, but x has $n element(s)",
        ),
    )
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter
    if isnothing(z)
        z = JuMP.@variable(m, base_name = "z_$(counter)")
    end
    pwaffine_formulation(m, x, z, counter, pwa, formulation, bounds)
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

# The input variables as an indexable collection
_inputs(x::JuMP.AbstractJuMPScalar) = [x]
_inputs(x) = collect(x)

function _var_bounds(v)
    JuMP.is_fixed(v) && return (JuMP.fix_value(v), JuMP.fix_value(v))
    return (JuMP.lower_bound(v), JuMP.upper_bound(v))
end

function _has_finite_bounds(v)
    v isa JuMP.VariableRef || return false
    JuMP.is_fixed(v) && return true
    return JuMP.has_lower_bound(v) && JuMP.has_upper_bound(v)
end

# Use the bounds given explicitly, otherwise the bounds of the JuMP variables
function _resolve_bounds(x, bounds, D)
    isnothing(bounds) || return bounds
    xs = _inputs(x)
    if length(xs) == D && all(_has_finite_bounds, xs)
        return [_var_bounds(v) for v ∈ xs]
    end
    return throw(
        ArgumentError(
            "Formulation requires finite bounds on x: pass `bounds = [(lo, hi), ...]` " *
            "or set bounds on the variables",
        ),
    )
end

# Model the piecewise affine function as a convex combination of its vertices
function pwaffine_formulation(
    m,
    x,
    z,
    counter,
    pwa::PWAFunc{C,D},
    ::LambdaFormulation,
    bounds,
) where {C,D}
    vertices = vertices_by_subsets(pwa, _resolve_bounds(x, bounds, D))
    xs = _inputs(x)
    V = eachindex(vertices)
    λ = @variable(m, [V], lower_bound = 0, base_name = "_λ_$(counter)")
    @constraint(
        m,
        [i = 1:D],
        xs[i] == sum(λ[v] * vertices[v][1][i] for v ∈ V),
        base_name = "pwa_$(counter)_x"
    )
    @constraint(
        m,
        z == sum(λ[v] * vertices[v][2] for v ∈ V),
        base_name = "pwa_$(counter)_z"
    )
    @constraint(m, sum(λ[v] for v ∈ V) == 1, base_name = "pwa_$(counter)_sum")
    return
end

constr(::Type{Convex}, m, z, p, x) = JuMP.@constraint(m, z ≥ dot(p.α, x) + p.β)
function constr(::Type{Concave}, m, z, p, x)
    JuMP.@constraint(m, z ≤ dot(-1 .* p.α, x) - p.β)
end

# Model the piecewise affine function directly as the maximum of its hyperplanes
function pwaffine_formulation(
    m,
    x,
    z,
    counter,
    pwa::PWAFunc{C,D},
    ::PsiFormulation,
    bounds,
) where {C,D}
    for (k, p) ∈ enumerate(pwa.planes)
        con = constr(C, m, z, p, x)
        JuMP.set_name(con, "pwa_$(counter)_$(k)")
    end
    if !isnothing(bounds)
        lo, hi = _check_bounds(bounds, D)
        xs = _inputs(x)
        @constraint(
            m,
            [i = 1:D],
            lo[i] ≤ xs[i] ≤ hi[i],
            base_name = "pwa_$(counter)_bounds"
        )
        # The function values over the box are bounded by the values at the vertices
        zmin, zmax = extrema(v[2] for v ∈ vertices_by_subsets(pwa, bounds))
        _bound_z!(m, z, zmin, zmax, counter)
    end
    return
end

# Bound z to [zmin, zmax]: as variable bounds for a free variable, otherwise as a constraint
function _bound_z!(m, z, zmin, zmax, counter)
    if z isa JuMP.VariableRef &&
       !JuMP.has_lower_bound(z) &&
       !JuMP.has_upper_bound(z) &&
       !JuMP.is_fixed(z)
        JuMP.set_lower_bound(z, zmin)
        JuMP.set_upper_bound(z, zmax)
    else
        @constraint(m, zmin ≤ z ≤ zmax, base_name = "pwa_$(counter)_zbounds")
    end
    return
end

# Sorted breakpoints (including the end points of the domain) of a 1D pwa function
function _breakpoints(x, pwa::PWAFunc{C,1}, bounds) where {C}
    verts = vertices_by_subsets(pwa, _resolve_bounds(x, bounds, 1))
    sort!(verts, by = v -> v[1][1])
    return _merge_segments([v[1][1] for v ∈ verts], [v[2] for v ∈ verts])
end

# Remove breakpoints where the slopes of the adjacent segments are (nearly) equal, to avoid
# redundant variables in the formulations
function _merge_segments(xp, zp; rtol = 1e-6)
    keep = [1]
    for i ∈ 2:(length(xp)-1)
        j = keep[end]
        s₁ = (zp[i] - zp[j]) / (xp[i] - xp[j])
        s₂ = (zp[i+1] - zp[i]) / (xp[i+1] - xp[i])
        if abs(s₂ - s₁) > rtol * max(1.0, abs(s₁), abs(s₂))
            push!(keep, i)
        end
    end
    length(xp) > 1 && push!(keep, length(xp))
    return xp[keep], zp[keep]
end

function pwaffine_formulation(
    m,
    x,
    z,
    counter,
    pwa::PWAFunc,
    ::DeltaFormulation,
    bounds,
)
    return throw(
        ArgumentError(
            "DeltaFormulation is only available for functions of one variable",
        ),
    )
end

# Model the 1D piecewise affine function by having a "bin" for each segment between
# consecutive breakpoints
function pwaffine_formulation(
    m,
    x,
    z,
    counter,
    pwa::PWAFunc{C,1},
    ::DeltaFormulation,
    bounds,
) where {C}
    xp, zp = _breakpoints(x, pwa, bounds)
    x1 = only(_inputs(x))
    S = 1:(length(xp)-1)
    width = diff(xp)
    slope = diff(zp) ./ width
    Δ = @variable(
        m,
        [s ∈ S],
        lower_bound = 0,
        upper_bound = width[s],
        base_name = "_Δ_$(counter)"
    )
    @constraint(
        m,
        x1 == xp[1] + sum(Δ[s] for s ∈ S),
        base_name = "pwa_$(counter)_x"
    )
    @constraint(
        m,
        z == zp[1] + sum(slope[s] * Δ[s] for s ∈ S),
        base_name = "pwa_$(counter)_z"
    )
    return
end

function pwaffine_formulation(
    m,
    x,
    z,
    counter,
    pwa::PWAFunc,
    ::PhiFormulation,
    bounds,
)
    return throw(
        ArgumentError(
            "PhiFormulation is only available for functions of one variable",
        ),
    )
end

# Model the 1D piecewise affine function as the first segment extended over the domain,
# corrected by the change in slope at each breakpoint
function pwaffine_formulation(
    m,
    x,
    z,
    counter,
    pwa::PWAFunc{C,1},
    ::PhiFormulation,
    bounds,
) where {C}
    xp, zp = _breakpoints(x, pwa, bounds)
    x1 = only(_inputs(x))
    @constraint(m, xp[1] ≤ x1 ≤ xp[end], base_name = "pwa_$(counter)_bounds")
    if length(xp) == 1
        @constraint(m, z == zp[1], base_name = "pwa_$(counter)_z")
        return
    end
    slope = diff(zp) ./ diff(xp)
    # Φ[s] measures how far x extends beyond breakpoint s
    S̃ = 2:(length(xp)-1)
    Φ = @variable(
        m,
        [s ∈ S̃],
        lower_bound = 0,
        upper_bound = xp[end] - xp[s],
        base_name = "_Φ_$(counter)"
    )
    @constraint(m, [s ∈ S̃], Φ[s] ≥ x1 - xp[s], base_name = "pwa_$(counter)_Φ")
    @constraint(
        m,
        z ==
        zp[1] +
        slope[1] * (x1 - xp[1]) +
        sum((slope[s] - slope[s-1]) * Φ[s] for s ∈ S̃),
        base_name = "pwa_$(counter)_z"
    )
    return
end

"""
    pwaffine(m::JuMP.Model, x, f::DCFunc{D}; z = nothing, direction = :graph,
             formulation = PsiFormulation(), selection = :bigm, bounds = nothing)

Add constraints to the JuMP model `m` for the variable `z` as the difference of
convex piecewise linear function `f` of the `D` variables `x`, and return `z`.
The function is modelled as z = z₁ + z₂ with z₁ its convex and z₂ its concave
part, each with the Ψ-formulation ([`PsiFormulation`](@ref), the only
`formulation` available for these functions). A part that is the maximum
(minimum) of planes needs binary variables to select the active plane for the
constraints that bound it from above (below).

- `direction = :graph` (default): z equals the function, with binary variables
  for both parts.
- `direction = :epigraph`: z is at least the function, which is sufficient when
  z is minimized or only bounded from below. The convex part needs no binary
  variables, as for a convex [`PWAFunc`](@ref).
- `direction = :hypograph`: z is at most the function, which is sufficient
  when z is maximized. The concave part needs no binary variables.

The active plane is selected by one binary variable per plane, and the
constraints for the other planes are relaxed with `selection = :bigm` by
big-M values computed from the bounds on `x`, or with
`selection = :indicator` by indicator constraints, which do not need
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
    formulation::PWAFormulation = PsiFormulation(),
    selection = :bigm,
    bounds = nothing,
) where {D}
    direction ∈ (:graph, :epigraph, :hypograph) || throw(
        ArgumentError(
            "Unsupported direction :$direction, use :graph, :epigraph or :hypograph",
        ),
    )
    formulation isa PsiFormulation || throw(
        ArgumentError(
            "Only PsiFormulation is available for difference of convex functions",
        ),
    )
    selection ∈ (:bigm, :indicator) || throw(
        ArgumentError(
            "Unsupported selection :$selection, use :bigm or :indicator",
        ),
    )
    xs = _inputs(x)
    length(xs) == D || throw(
        ArgumentError(
            "The function has $D inputs, got $(length(xs)) variables",
        ),
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
    box = if selection == :bigm && (exact₁ || exact₂)
        _check_bounds(_resolve_bounds(x, bounds, D), D)
    else
        nothing
    end
    z₁ = _dc_part!(
        m,
        xs,
        f.convex,
        direction != :epigraph,
        selection,
        box,
        "$(counter)_cvx",
    )
    z₂ = _dc_part!(
        m,
        xs,
        f.concave,
        direction != :hypograph,
        selection,
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
    selection,
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
    M = selection == :bigm ? _dc_big_M(pwa, box) : nothing
    for k ∈ 1:K
        con = if selection == :indicator && C == Convex
            JuMP.@constraint(m, u[k] --> {zᵢ ≤ v[k]})
        elseif selection == :indicator
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

# Big-M for each plane k of a part: an upper bound on the distance from the
# part to plane k over the box (lo, hi), max over planes i and x in the box of
# (αᵢ - αₖ)⋅x + βᵢ - βₖ. The same expression applies to the concave part,
# where the stored planes are negated.
function _dc_big_M(pwa::PWAFunc{C,D}, (lo, hi)) where {C,D}
    return map(pwa.planes) do pₖ
        return maximum(pwa.planes) do pᵢ
            return pᵢ.β - pₖ.β + sum(
                max((pᵢ.α[j] - pₖ.α[j]) * lo[j], (pᵢ.α[j] - pₖ.α[j]) * hi[j])
                for j ∈ 1:D
            )
        end
    end
end
