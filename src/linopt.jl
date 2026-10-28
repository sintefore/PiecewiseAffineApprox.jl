"""
    PWARecord

The variables and constraints that [`add_pwaffine`](@ref) added to a model for
a piecewise affine function: fields `model`, `name`, `output` (the output `z`),
`variables` and `constraints`. The variables include the output if it was
created, and the constraints exclude variable bounds and integrality, which
are removed with the variables. Use `JuMP.delete(model, record)` to remove the
function from the model again; bounds that were set on an output given by
the user are constraints, and are removed as well.
"""
struct PWARecord{M<:JuMP.AbstractModel}
    model::M
    name::String
    output::Base.RefValue{Any}
    variables::Vector{JuMP.VariableRef}
    constraints::Vector{JuMP.ConstraintRef}
end

function Base.show(io::IO, r::PWARecord)
    return print(
        io,
        "PWARecord \"$(r.name)\" with $(length(r.variables)) variables and ",
        "$(length(r.constraints)) constraints",
    )
end

function JuMP.delete(model::JuMP.AbstractModel, r::PWARecord)
    model === r.model ||
        throw(ArgumentError("The record does not belong to the model"))
    for c ∈ r.constraints
        JuMP.is_valid(model, c) && JuMP.delete(model, c)
    end
    vars = filter(v -> JuMP.is_valid(model, v), r.variables)
    isempty(vars) || JuMP.delete(model, vars)
    empty!(r.constraints)
    empty!(r.variables)
    return
end

mutable struct PWLData
    counter::Int
    recording::Union{Nothing,PWARecord}
    PWLData() = new(0, nothing)
end

function initPWL!(m::JuMP.Model)
    if !haskey(m.ext, :PWL)
        m.ext[:PWL] = PWLData()
    end
    return nothing
end

# Record the variables or constraints x in the record being built, if any,
# and return x
function _rec!(m, x)
    data = get(m.ext, :PWL, nothing)
    if data isa PWLData && !isnothing(data.recording)
        _push_record!(data.recording, x)
    end
    return x
end
_push_record!(r, x::JuMP.VariableRef) = push!(r.variables, x)
_push_record!(r, x::JuMP.ConstraintRef) = push!(r.constraints, x)
_push_record!(r, x::AbstractArray) = foreach(y -> _push_record!(r, y), x)
_push_record!(r, x) = nothing

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
    Direction

How the output `z` of [`pwaffine`](@ref) relates to the function value:
[`Graph`](@ref), [`Epigraph`](@ref) or [`Hypograph`](@ref).
"""
abstract type Direction end

"""
    Graph()

The output equals the function, `z == f(x)`. The default for a
[`DCFunc`](@ref). For a convex or concave [`PWAFunc`](@ref), this needs
binary variables, as for a `DCFunc`.
"""
struct Graph <: Direction end

"""
    Epigraph()

The output is at least the function, `z ≥ f(x)`, which suffices when `z` is
minimized. The default for a convex [`PWAFunc`](@ref), where it needs no
binary variables.
"""
struct Epigraph <: Direction end

"""
    Hypograph()

The output is at most the function, `z ≤ f(x)`, which suffices when `z` is
maximized. The default for a concave [`PWAFunc`](@ref), where it needs no
binary variables.
"""
struct Hypograph <: Direction end

# The direction that needs no binary variables for the curvature
_natural_direction(::Type{Convex}) = Epigraph()
_natural_direction(::Type{Concave}) = Hypograph()

"""
    pwaffine(m::JuMP.Model, x, pwa::PWAFunc{C,D}; z, formulation::PWAFormulation, bounds, direction) where {C,D}

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
- `name`: the label in the names of the variables and constraints added, e.g. `z_name`
  (default: a counter).
- `bounds`: Box `[(lo₁, hi₁), …, (lo_D, hi_D)]` restricting the domain of `x`. Required by the
  λ, Δ and Φ formulations, and taken from the [`domain`](@ref) of `pwa` if not given, or
  else from the bounds of the JuMP variables `x` if all are finite (or fixed). For the Ψ
  formulation the bounds are optional. If given, they are added as constraints on `x`.
  Over the box that `x` is restricted to, the output `z` is bounded by the smallest value
  of the function unless the direction is `Hypograph()`, and by the largest value unless it
  is `Epigraph()`: as variable bounds for a new `z`, and as a constraint for a given one.
- `direction::Direction`: [`Epigraph`](@ref) for a convex and [`Hypograph`](@ref) for a
  concave function by default, which the formulations above model without binary
  variables. With [`Graph`](@ref), or the other one-sided direction, the function is
  modelled as a [`DCFunc`](@ref) with binary variables, see `pwaffine` for `DCFunc`,
  which also takes the `selection` keyword.

# Returns
- The JuMP variable `z` representing the piecewise affine function value. Use
  [`add_pwaffine`](@ref) to also get a record of what was added, and see the call of a
  [`PWAFunc`](@ref) for a shorter form, e.g. `pwa(x)`.

# Note
In the default direction, all formulations are pure linear programs (no binary
variables). They model z ≥ f(x) for convex and z ≤ f(x) for concave functions, i.e. the
relation is only tight if the objective pushes `z` towards the function (minimizing a
convex or maximizing a concave function). The Δ and Φ formulations only work for 1D
problems and will error for higher dimensions. With interior point solvers for nonlinear
models, the λ and Δ formulations are typically much faster than Ψ (see Coffrin et al.
(2021), arXiv:2005.14087).

# Example
```julia
m = Model(HiGHS.Optimizer)
@variable(m, -1 ≤ x ≤ 1)
pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), MILP(optimizer = optimizer, planes=5))
z = pwaffine(m, x, pwa; formulation = LambdaFormulation())
```
"""
function pwaffine(m::JuMP.Model, x, pwa::PWAFunc; kwargs...)
    return first(add_pwaffine(m, x, pwa; kwargs...))
end

"""
    add_pwaffine(m::JuMP.Model, x, f::AbstractPWAFunc; name = nothing, kwargs...)

As [`pwaffine`](@ref), but returning the tuple `(z, record)` with the
[`PWARecord`](@ref) of the variables and constraints that were added, which
can be removed again with `JuMP.delete(m, record)`. The variables and
constraints are named with `name` (by default a counter), e.g. `z_name`. If
the function cannot be added, an error is thrown and the model is left
unchanged.
"""
function add_pwaffine(
    m::JuMP.Model,
    x,
    f::AbstractPWAFunc;
    name::Union{Nothing,AbstractString} = nothing,
    kwargs...,
)
    initPWL!(m)
    data = m.ext[:PWL]
    isnothing(data.recording) ||
        error("add_pwaffine cannot be called while another function is added")
    data.counter += 1
    label = isnothing(name) ? string(data.counter) : String(name)
    record = PWARecord(
        m,
        label,
        Ref{Any}(nothing),
        JuMP.VariableRef[],
        JuMP.ConstraintRef[],
    )
    data.recording = record
    try
        record.output[] = _pwaffine!(m, x, f, label; kwargs...)
    catch
        data.recording = nothing
        JuMP.delete(m, record)
        rethrow()
    finally
        data.recording = nothing
    end
    return record.output[], record
end

function _pwaffine!(
    m::JuMP.Model,
    x,
    pwa::PWAFunc{C,D},
    counter;
    z = nothing,
    formulation::PWAFormulation = PsiFormulation(),
    bounds = nothing,
    direction::Direction = _natural_direction(C),
    selection = :bigm,
) where {C,D}
    if direction != _natural_direction(C)
        # Binary variables are needed, as for a difference of convex function
        return _pwaffine!(
            m,
            x,
            _as_dc(pwa),
            counter;
            z,
            direction,
            formulation,
            selection,
            bounds,
        )
    end
    n = length(_inputs(x))
    n == D || throw(
        ArgumentError(
            "The pwa function has dimension $D, but x has $n element(s)",
        ),
    )
    created = isnothing(z)
    if created
        z = _rec!(m, JuMP.@variable(m, base_name = "z_$(counter)"))
    end
    pwaffine_formulation(m, x, z, counter, pwa, formulation, bounds)
    # Bound the output over the box that x is restricted to, by the explicit
    # bounds (constraints of the Ψ-formulation), the domain (for the
    # formulations that restrict x to it) or the bounds of the variables
    domain = formulation isa PsiFormulation ? nothing : pwa.domain
    box = _known_box(x, bounds, D, domain)
    if !isnothing(box)
        zmin, zmax = extrema(v[2] for v ∈ vertices_by_subsets(pwa, box))
        _bound_output!(m, z, created, zmin, zmax, direction, counter)
    end
    return z
end

function pwaffine(
    m::JuMP.Model,
    x,
    fevals::FunctionEvaluations,
    curvature::Curvature,
    a::ApproximationAlgorithm;
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

# A convex or concave function as a difference of convex function with an
# empty other part
function _as_dc(pwa::PWAFunc{Convex,D,T}) where {D,T}
    return DCFunc{D,T}(pwa, PWAFunc{Concave,D,T}(Plane{D,T}[]), pwa.domain)
end
function _as_dc(pwa::PWAFunc{Concave,D,T}) where {D,T}
    return DCFunc{D,T}(PWAFunc{Convex,D,T}(Plane{D,T}[]), pwa, pwa.domain)
end

# Bound the output z by the smallest and largest function value, from the
# sides that the direction implies: as variable bounds for an output created
# here, and otherwise as a constraint
function _bound_output!(m, z, created, zmin, zmax, direction, counter)
    lower = !(direction isa Hypograph)
    upper = !(direction isa Epigraph)
    if created
        lower && JuMP.set_lower_bound(z, zmin)
        upper && JuMP.set_upper_bound(z, zmax)
    elseif lower && upper
        _rec!(m, @constraint(m, zmin ≤ z ≤ zmax, base_name = "pwa_$(counter)_zbounds"))
    elseif lower
        _rec!(m, @constraint(m, z ≥ zmin, base_name = "pwa_$(counter)_zbounds"))
    elseif upper
        _rec!(m, @constraint(m, z ≤ zmax, base_name = "pwa_$(counter)_zbounds"))
    end
    return
end

# The box that x is known to lie in: the explicit bounds, otherwise the domain
# (if given), otherwise the bounds of the variables, or nothing
function _known_box(x, bounds, D, domain)
    isnothing(bounds) || return bounds
    isnothing(domain) || return domain
    xs = _inputs(x)
    if length(xs) == D && all(_has_finite_bounds, xs)
        return [_var_bounds(v) for v ∈ xs]
    end
    return nothing
end

# Use the bounds given explicitly, otherwise the domain of the function, and
# otherwise the bounds of the JuMP variables
function _resolve_bounds(x, bounds, D, domain = nothing)
    isnothing(bounds) || return bounds
    isnothing(domain) || return domain
    xs = _inputs(x)
    if length(xs) == D && all(_has_finite_bounds, xs)
        return [_var_bounds(v) for v ∈ xs]
    end
    return throw(
        ArgumentError(
            "Formulation requires finite bounds on x: pass `bounds = [(lo, hi), ...]`, " *
            "set bounds on the variables, or fit the function with approx, which " *
            "records its domain",
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
    vertices = vertices_by_subsets(pwa, _resolve_bounds(x, bounds, D, pwa.domain))
    xs = _inputs(x)
    V = eachindex(vertices)
    λ = _rec!(m, @variable(m, [V], lower_bound = 0, base_name = "_λ_$(counter)"))
    _rec!(m, @constraint(
        m,
        [i = 1:D],
        xs[i] == sum(λ[v] * vertices[v][1][i] for v ∈ V),
        base_name = "pwa_$(counter)_x"
    ))
    _rec!(m, @constraint(
        m,
        z == sum(λ[v] * vertices[v][2] for v ∈ V),
        base_name = "pwa_$(counter)_z"
    ))
    _rec!(m, @constraint(m, sum(λ[v] for v ∈ V) == 1, base_name = "pwa_$(counter)_sum"))
    return
end

constr(::Type{Convex}, m, z, p, x) = _rec!(m, JuMP.@constraint(m, z ≥ dot(p.α, x) + p.β))
function constr(::Type{Concave}, m, z, p, x)
    _rec!(m, JuMP.@constraint(m, z ≤ dot(-1 .* p.α, x) - p.β))
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
        _rec!(m, @constraint(
            m,
            [i = 1:D],
            lo[i] ≤ xs[i] ≤ hi[i],
            base_name = "pwa_$(counter)_bounds"
        ))
    end
    return
end

# Sorted breakpoints (including the end points of the domain) of a 1D pwa function
function _breakpoints(x, pwa::PWAFunc{C,1}, bounds) where {C}
    verts = vertices_by_subsets(pwa, _resolve_bounds(x, bounds, 1, pwa.domain))
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
    Δ = _rec!(m, @variable(
        m,
        [s ∈ S],
        lower_bound = 0,
        upper_bound = width[s],
        base_name = "_Δ_$(counter)"
    ))
    _rec!(m, @constraint(
        m,
        x1 == xp[1] + sum(Δ[s] for s ∈ S),
        base_name = "pwa_$(counter)_x"
    ))
    _rec!(m, @constraint(
        m,
        z == zp[1] + sum(slope[s] * Δ[s] for s ∈ S),
        base_name = "pwa_$(counter)_z"
    ))
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
    _rec!(m, @constraint(m, xp[1] ≤ x1 ≤ xp[end], base_name = "pwa_$(counter)_bounds"))
    if length(xp) == 1
        _rec!(m, @constraint(m, z == zp[1], base_name = "pwa_$(counter)_z"))
        return
    end
    slope = diff(zp) ./ diff(xp)
    # Φ[s] measures how far x extends beyond breakpoint s
    S̃ = 2:(length(xp)-1)
    Φ = _rec!(m, @variable(
        m,
        [s ∈ S̃],
        lower_bound = 0,
        upper_bound = xp[end] - xp[s],
        base_name = "_Φ_$(counter)"
    ))
    _rec!(m, @constraint(m, [s ∈ S̃], Φ[s] ≥ x1 - xp[s], base_name = "pwa_$(counter)_Φ"))
    _rec!(m, @constraint(
        m,
        z ==
        zp[1] +
        slope[1] * (x1 - xp[1]) +
        sum((slope[s] - slope[s-1]) * Φ[s] for s ∈ S̃),
        base_name = "pwa_$(counter)_z"
    ))
    return
end

"""
    pwaffine(m::JuMP.Model, x, f::DCFunc{D}; z = nothing, direction = Graph(),
             formulation = PsiFormulation(), selection = :bigm, bounds = nothing)

Add constraints to the JuMP model `m` for the variable `z` as the difference of
convex piecewise linear function `f` of the `D` variables `x`, and return `z`.
Binary variables are needed where the function is bounded against its
curvature, depending on `direction`:

- `direction = Graph()` (default): z equals the function.
- `direction = Epigraph()`: z is at least the function, which is sufficient
  when z is minimized or only bounded from below. The convex part, or the
  convex kinks, need no binary variables, as for a convex [`PWAFunc`](@ref).
- `direction = Hypograph()`: z is at most the function, which is sufficient
  when z is maximized. The concave part, or the concave kinks, need no binary
  variables.

With the Ψ-formulation ([`PsiFormulation`](@ref), default, any dimension), the
function is modelled as z = z₁ + z₂ with z₁ its convex and z₂ its concave part,
and the active plane of a part is selected by one binary variable per plane.
The constraints for the other planes are relaxed with `selection = :bigm` by
big-M values computed from the bounds on `x`, or with `selection = :indicator`
by indicator constraints, which do not need bounds but require solver support.

With the Δ- and Φ-formulations ([`DeltaFormulation`](@ref) and
[`PhiFormulation`](@ref), one variable only), the function is modelled by its
breakpoints over the bounds, as for a [`PWAFunc`](@ref), with one binary
variable at each kink of the wrong curvature for the direction: the
formulations are then ideal and no big-M values are needed. The convex and
concave kinks of the function are found from the vertices of its parts.

The bounds are given with `bounds` as a vector of `(lower, upper)` tuples, or
taken from the [`domain`](@ref) of `f`, or else from the lower and upper bounds
of the variables `x`. Explicit bounds, and the domain when used for big-M
values, are added as constraints on `x`. The output `z` is bounded over the
box that `x` is restricted to, from the sides that the direction implies, as
for a [`PWAFunc`](@ref). The keyword `name` and the record are as for
[`add_pwaffine`](@ref).
"""
function pwaffine(m::JuMP.Model, x, f::DCFunc; kwargs...)
    return first(add_pwaffine(m, x, f; kwargs...))
end

function _pwaffine!(
    m::JuMP.Model,
    x,
    f::DCFunc{D},
    counter;
    z = nothing,
    direction::Direction = Graph(),
    formulation::PWAFormulation = PsiFormulation(),
    selection = :bigm,
    bounds = nothing,
) where {D}
    if formulation isa LambdaFormulation ||
       (!(formulation isa PsiFormulation) && D != 1)
        throw(
            ArgumentError(
                "$(typeof(formulation)) is not available for difference of " *
                "convex functions of $D variables, use PsiFormulation",
            ),
        )
    end
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
    created = isnothing(z)
    if created
        z = _rec!(m, JuMP.@variable(m, base_name = "z_$(counter)"))
    end
    if !(formulation isa PsiFormulation)
        xp, zp = _dc_breakpoints(f, _resolve_bounds(x, bounds, 1, f.domain))
        _dc_formulation_1d!(
            m,
            only(xs),
            z,
            counter,
            xp,
            zp,
            formulation,
            direction,
        )
        # x is restricted to the breakpoints
        _bound_output!(m, z, created, minimum(zp), maximum(zp), direction, counter)
        return z
    end
    # The box containing x, needed for the big-M values when a part has
    # binary variables
    exact₁ = !(direction isa Epigraph) && _planes(f.convex) > 1
    exact₂ = !(direction isa Hypograph) && _planes(f.concave) > 1
    box = if selection == :bigm && (exact₁ || exact₂)
        _check_bounds(_resolve_bounds(x, bounds, D, f.domain), D)
    else
        nothing
    end
    # x is restricted to the box of the big-M values, and to explicit bounds;
    # the bounds of the variables need no constraints
    restricted = !isnothing(box) ? box :
        !isnothing(bounds) ? _check_bounds(bounds, D) : nothing
    if !isnothing(restricted) && (!isnothing(bounds) || !isnothing(f.domain))
        lo, hi = restricted
        _rec!(
            m,
            @constraint(
                m,
                [i = 1:D],
                lo[i] ≤ xs[i] ≤ hi[i],
                base_name = "pwa_$(counter)_bounds"
            )
        )
    elseif isnothing(restricted)
        known = _known_box(x, nothing, D, nothing)
        restricted = isnothing(known) ? nothing : _check_bounds(known, D)
    end
    z₁ = _dc_part!(
        m,
        xs,
        f.convex,
        !(direction isa Epigraph),
        selection,
        box,
        "$(counter)_cvx",
    )
    z₂ = _dc_part!(
        m,
        xs,
        f.concave,
        !(direction isa Hypograph),
        selection,
        box,
        "$(counter)_ccv",
    )
    con = _rec!(m, JuMP.@constraint(m, z == z₁ + z₂))
    JuMP.set_name(con, "pwa_$(counter)")
    if !isnothing(restricted)
        zmin, zmax = _dc_range(f, collect(zip(restricted...)))
        _bound_output!(m, z, created, zmin, zmax, direction, counter)
    end
    return z
end

# Bounds on the values of a difference of convex function over a box: the sums
# of the smallest and largest values of the parts, which are attained at the
# vertices of the parts over the box
function _dc_range(f::DCFunc{D}, box) where {D}
    zmin, zmax = 0.0, 0.0
    for part ∈ (f.convex, f.concave)
        isempty(part.planes) && continue
        v = [w[2] for w ∈ vertices_by_subsets(part, box)]
        zmin += minimum(v)
        zmax += maximum(v)
    end
    return zmin, zmax
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
    zᵢ = _rec!(m, JuMP.@variable(m, base_name = "z_$(name)"))
    v = [_plane_expr(C, p, x) for p ∈ pwa.planes]
    # The constraints giving the maximum (minimum) of the planes
    for k ∈ 1:K
        con =
            C == Convex ? _rec!(m, JuMP.@constraint(m, zᵢ ≥ v[k])) :
            _rec!(m, JuMP.@constraint(m, zᵢ ≤ v[k]))
        JuMP.set_name(con, "pwa_$(name)_$(k)")
    end
    exact || return zᵢ
    if K == 1
        con =
            C == Convex ? _rec!(m, JuMP.@constraint(m, zᵢ ≤ v[1])) :
            _rec!(m, JuMP.@constraint(m, zᵢ ≥ v[1]))
        JuMP.set_name(con, "pwa_$(name)_active_1")
        return zᵢ
    end
    u = _rec!(m, JuMP.@variable(m, [1:K], Bin, base_name = "u_$(name)"))
    con = _rec!(m, JuMP.@constraint(m, sum(u) == 1))
    JuMP.set_name(con, "pwa_$(name)_select")
    M = selection == :bigm ? _dc_big_M(pwa, box) : nothing
    for k ∈ 1:K
        con = if selection == :indicator && C == Convex
            _rec!(m, JuMP.@constraint(m, u[k] --> {zᵢ ≤ v[k]}))
        elseif selection == :indicator
            _rec!(m, JuMP.@constraint(m, u[k] --> {zᵢ ≥ v[k]}))
        elseif C == Convex
            _rec!(m, JuMP.@constraint(m, zᵢ ≤ v[k] + M[k] * (1 - u[k])))
        else
            _rec!(m, JuMP.@constraint(m, zᵢ ≥ v[k] - M[k] * (1 - u[k])))
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

# The breakpoints of a one-dimensional difference of convex function over the
# bounds: the end points and the kinks of both parts, with the values of the
# function, after merging the breakpoints where the slope does not change
function _dc_breakpoints(f::DCFunc{1}, bounds)
    lo, hi = _check_bounds(bounds, 1)
    xp = [lo[1], hi[1]]
    for part ∈ (f.convex, f.concave)
        isempty(part.planes) && continue
        append!(xp, v[1][1] for v ∈ vertices_by_subsets(part, bounds))
    end
    sort!(xp)
    τ = 1e-9 * max(1.0, hi[1] - lo[1])
    unique = [xp[1]]
    for xᵢ ∈ xp
        xᵢ - unique[end] > τ && push!(unique, xᵢ)
    end
    zp = [evaluate(f, (xᵢ,)) for xᵢ ∈ unique]
    return _merge_segments(unique, zp)
end

# The interior breakpoints that need a binary variable: the kinks of the
# wrong curvature for the direction, i.e. where the slope decreases for the
# epigraph and increases for the hypograph, and all kinks for the graph
function _dc_binary_kinks(slope, direction::Direction)
    S̃ = 2:length(slope)
    direction isa Graph && return collect(S̃)
    if direction isa Epigraph
        return [s for s ∈ S̃ if slope[s] < slope[s-1]]
    end
    return [s for s ∈ S̃ if slope[s] > slope[s-1]]
end

# Model the one-dimensional function with breakpoints xp and values zp by the
# Φ-formulation: the first segment extended over the domain, corrected by the
# change in slope at each breakpoint s times Φ[s] = max(x - xp[s], 0). A
# convex kink in the epigraph direction (and a concave one in the hypograph
# direction) only needs Φ[s] ≥ x - xp[s] and Φ[s] ≥ 0. The other kinks need a
# binary variable y[s] = 1 if x ≥ xp[s], which makes Φ[s] exact.
function _dc_formulation_1d!(
    m,
    x1,
    z,
    counter,
    xp,
    zp,
    ::PhiFormulation,
    direction,
)
    _rec!(m, @constraint(m, xp[1] ≤ x1 ≤ xp[end], base_name = "pwa_$(counter)_bounds"))
    if length(xp) == 1
        _rec!(m, @constraint(m, z == zp[1], base_name = "pwa_$(counter)_z"))
        return
    end
    slope = diff(zp) ./ diff(xp)
    S̃ = 2:(length(xp)-1)
    Φ = _rec!(m, @variable(
        m,
        [s ∈ S̃],
        lower_bound = 0,
        upper_bound = xp[end] - xp[s],
        base_name = "_Φ_$(counter)"
    ))
    _rec!(m, @constraint(m, [s ∈ S̃], Φ[s] ≥ x1 - xp[s], base_name = "pwa_$(counter)_Φ"))
    _rec!(m, @constraint(
        m,
        z ==
        zp[1] +
        slope[1] * (x1 - xp[1]) +
        sum((slope[s] - slope[s-1]) * Φ[s] for s ∈ S̃),
        base_name = "pwa_$(counter)_z"
    ))
    B = _dc_binary_kinks(slope, direction)
    isempty(B) && return
    y = _rec!(m, @variable(m, [s ∈ B], Bin, base_name = "_y_$(counter)"))
    _rec!(m, @constraint(
        m,
        [s ∈ B],
        Φ[s] ≤ (xp[end] - xp[s]) * y[s],
        base_name = "pwa_$(counter)_Φup"
    ))
    _rec!(m, @constraint(
        m,
        [s ∈ B],
        Φ[s] ≤ x1 - xp[s] + (xp[s] - xp[1]) * (1 - y[s]),
        base_name = "pwa_$(counter)_Φx"
    ))
    # x ≥ xp[s+1] implies x ≥ xp[s]
    _rec!(m, @constraint(
        m,
        [i ∈ 1:(length(B)-1)],
        y[B[i]] ≥ y[B[i+1]],
        base_name = "pwa_$(counter)_order"
    ))
    return
end

# Model the one-dimensional function with breakpoints xp and values zp by the
# Δ-formulation: a variable for how much of each segment is used, with x and z
# integrating the segments. Within a run of segments with increasing slopes
# (decreasing for the hypograph direction), the segments are used in order
# when z is pushed towards the function. The other kinks split the segments
# into such runs, and a binary variable y[s] = 1 if x is beyond the kink
# forces all earlier segments to be full and all later segments to be empty.
function _dc_formulation_1d!(
    m,
    x1,
    z,
    counter,
    xp,
    zp,
    ::DeltaFormulation,
    direction,
)
    if length(xp) == 1
        _rec!(m, @constraint(m, x1 == xp[1], base_name = "pwa_$(counter)_x"))
        _rec!(m, @constraint(m, z == zp[1], base_name = "pwa_$(counter)_z"))
        return
    end
    S = 1:(length(xp)-1)
    width = diff(xp)
    slope = diff(zp) ./ width
    Δ = _rec!(m, @variable(
        m,
        [s ∈ S],
        lower_bound = 0,
        upper_bound = width[s],
        base_name = "_Δ_$(counter)"
    ))
    _rec!(m, @constraint(
        m,
        x1 == xp[1] + sum(Δ[s] for s ∈ S),
        base_name = "pwa_$(counter)_x"
    ))
    _rec!(m, @constraint(
        m,
        z == zp[1] + sum(slope[s] * Δ[s] for s ∈ S),
        base_name = "pwa_$(counter)_z"
    ))
    B = _dc_binary_kinks(slope, direction)
    isempty(B) && return
    y = _rec!(m, @variable(m, [s ∈ B], Bin, base_name = "_y_$(counter)"))
    _rec!(m, @constraint(
        m,
        [s ∈ B, j ∈ S; j ≥ s],
        Δ[j] ≤ width[j] * y[s],
        base_name = "pwa_$(counter)_Δempty"
    ))
    _rec!(m, @constraint(
        m,
        [s ∈ B, j ∈ S; j < s],
        Δ[j] ≥ width[j] * y[s],
        base_name = "pwa_$(counter)_Δfull"
    ))
    return
end

# The inputs of a call: one collection of inputs, or the inputs themselves
_call_inputs(x::Tuple{Union{Tuple,AbstractVector}}) = Tuple(only(x))
_call_inputs(x::Tuple) = x

# The model of the inputs, or nothing if they are all numbers
function _owner_model(xs)
    model = nothing
    for xi ∈ xs
        xi isa Real && continue
        mi = JuMP.owner_model(xi)
        if isnothing(model)
            model = mi
        elseif mi !== model
            throw(ArgumentError("The inputs belong to different models"))
        end
    end
    return model
end

"""
    (f::AbstractPWAFunc)(x...; kwargs...)

Evaluate `f` at a point given as numbers, e.g. `f(0.5)` or `f(1.0, 2.0)`, or,
if the inputs are JuMP variables or expressions, add `f` to their model with
[`pwaffine`](@ref) and return the output. The inputs can also be given as one
vector or tuple, e.g. `f(x)` for `x = @variable(m, [1:2])`. The keyword
arguments are passed to `pwaffine`; see [`Formulated`](@ref) to store them
with the function. Each call adds a new formulation, so a function of the same
inputs that is used several times should be added once, e.g. with
`@expression(m, c, f(x))`.

```julia
cost = approx(x -> x[1]^2, [(0, 1)], Convex(), Interpol(planes = 4))
@variable(m, 0 ≤ x[1:3] ≤ 1)
@objective(m, Min, sum(cost(x[i]) for i ∈ 1:3))
```
"""
function (f::AbstractPWAFunc)(x...; kwargs...)
    xs = _call_inputs(x)
    model = _owner_model(xs)
    if isnothing(model)
        isempty(kwargs) ||
            throw(ArgumentError("Keyword arguments are only used with JuMP inputs"))
        return evaluate(f, xs)
    end
    return pwaffine(model, length(xs) == 1 ? only(xs) : collect(xs), f; kwargs...)
end

"""
    Formulated(f::AbstractPWAFunc; formulation, direction, selection, bounds)

The function `f` with the options of [`pwaffine`](@ref) for adding it to a
model, which is called like the function: `g(x)` evaluates `g` for numbers
and adds it with these options for JuMP inputs, see the call of
[`PWAFunc`](@ref). Options given in a call override the stored ones. It can
also be passed to [`pwaffine`](@ref), [`add_pwaffine`](@ref) and
[`evaluate`](@ref).

```julia
cost = Formulated(approx(sin, [(0, 2π)], DC(), Interpol(planes = 4));
    formulation = DeltaFormulation(), direction = Epigraph())
@objective(m, Min, cost(x))
```
"""
struct Formulated{F<:AbstractPWAFunc,K<:NamedTuple}
    func::F
    options::K
end
const _FORMULATED_OPTIONS = (:formulation, :direction, :selection, :bounds)
function Formulated(f::AbstractPWAFunc; options...)
    for k ∈ keys(options)
        k ∈ _FORMULATED_OPTIONS || throw(
            ArgumentError(
                "Unknown option $k, use one of " * join(_FORMULATED_OPTIONS, ", "),
            ),
        )
    end
    return Formulated(f, NamedTuple(options))
end
function Base.show(io::IO, g::Formulated)
    print(io, "Formulated(")
    summary(io, g.func)
    for (k, v) ∈ pairs(g.options)
        print(io, "; ", k, " = ", v)
    end
    return print(io, ")")
end
function (g::Formulated)(x...; kwargs...)
    # The stored options are for adding the function to a model only
    isnothing(_owner_model(_call_inputs(x))) && return g.func(x...; kwargs...)
    return g.func(x...; g.options..., kwargs...)
end
evaluate(g::Formulated, x) = evaluate(g.func, x)
domain(g::Formulated) = domain(g.func)
function pwaffine(m::JuMP.Model, x, g::Formulated; kwargs...)
    return pwaffine(m, x, g.func; g.options..., kwargs...)
end
function add_pwaffine(m::JuMP.Model, x, g::Formulated; kwargs...)
    return add_pwaffine(m, x, g.func; g.options..., kwargs...)
end
