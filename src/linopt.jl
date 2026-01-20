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

abstract type PWAFormulation end

struct Ψ_Formulation <: PWAFormulation end
struct λ_Formulation <: PWAFormulation end
struct Δ_Formulation <: PWAFormulation end
struct Φ_Formulation <: PWAFormulation end


const VarOrAff = Union{JuMP.VariableRef,JuMP.AffExpr}

"""
    pwaffine(m::JuMP.Model, x, pwa::PWAFunc{C,D}; z, formulation::PWAFormulation) where {C,D}

Add constraints to JuMP-model `m` for JuMP-variable `z` as a
piecewise linear function/approximation `pwa` of JuMP-variables `x`.

# Arguments
- `m::JuMP.Model`: The JuMP model to add constraints to
- `x`: JuMP variable(s) representing the input to the piecewise affine function
- `pwa::PWAFunc{C,D}`: The piecewise affine function to model

# Keyword Arguments
- `z::Union{Nothing, VarOrAff}`: JuMP variable representing the output. If `nothing`, a new variable is created
- `formulation::PWAFormulation`: The formulation to use. Options are:
  - `Ψ_Formulation()` (default): Direct plane representation. Works for any dimension. No preprocessing required.
  - `λ_Formulation()`: Convex combination of vertices. Works for any dimension. Requires vertex enumeration.
  - `Δ_Formulation()`: Specialized for 1D. Uses ordered breakpoint structure. Requires breakpoint enumeration.
  - `Φ_Formulation()`: Alternative 1D formulation. Requires breakpoint enumeration.

# Returns
- The JuMP variable `z` representing the piecewise affine function value

# Note
All formulations are pure linear programs (no binary variables). The Δ and Φ formulations
only work for 1D problems and will error for higher dimensions.

# Example
```julia
m = Model(HiGHS.Optimizer)
@variable(m, x)
pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), MILP(optimizer = optimizer, planes=5))
z = pwaffine(m, x, pwa; formulation = λ_Formulation())
```
"""
function pwaffine(
    m::JuMP.Model,
    x,
    pwa::PWAFunc{C,D};
    z::Union{Nothing, VarOrAff} = nothing,
    formulation = Ψ_Formulation()
) where {C,D}

    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter
    if isnothing(z)
        z = JuMP.@variable(m, base_name = "z_$(counter)")
    end
    pwaffine_formulation(m, x, z, counter, pwa, formulation)
    return z
end

function pwaffine(
    m::JuMP.Model,
    x,
    fevals::FunctionEvaluations,
    curvature::Curvature,
    a::Algorithm;
    z = nothing,
    formulation = Ψ_Formulation()
)
    return pwaffine(m, x, approx(fevals, curvature, a); z = z, formulation = formulation)
end

# Model the piecewise affine function as a convex combination of its vertices
function pwaffine_formulation(m, x, z, counter, pwa::PWAFunc{C,D}, _::λ_Formulation) where {C,D}
    vertices = PiecewiseAffineApprox.vertices_by_subsets(pwa)
    V = 1 : length(vertices)
    λ = @variable(m, [V], lower_bound = 0, upper_bound = 1, base_name = "_λ_$(counter)")
    @constraint(m, x .== sum(λ[i] * v[1] for (i,v) in enumerate(vertices)))
    @constraint(m, z == sum(λ[i] * v[2] for (i,v) in enumerate(vertices)))
    @constraint(m, sum(λ[v] for v in V) == 1)
    return
end

constr(::Type{Convex}, m, z, p, x) = JuMP.@constraint(m, z ≥ dot(p.α, x) + p.β)
function constr(::Type{Concave}, m, z, p, x)
    JuMP.@constraint(m, z ≤ dot(-1 .* p.α, x) - p.β)
end

# Model the piecewise affine function directly as the maximum of its hyperplanes
function pwaffine_formulation(m, x, z, counter, pwa::PWAFunc{C,D}, _::Ψ_Formulation) where {C,D}
    for p in pwa.planes
        constr(C, m, z, p, x)
    end
    return
end

function pwaffine_formulation(m, x, z, counter, pwa::PWAFunc{C,D}, _::Δ_Formulation) where {C,D}
    error("Δ_Formulation for piecewise affine functions only possible for 1D functions")
end

# Model the 1D piecewise affine function by having a "bin" between consecutive breakpoints
function pwaffine_formulation(m, x, z, counter, pwa::PWAFunc{C,1}, _::Δ_Formulation) where {C}
    vp = PiecewiseAffineApprox.vertices_by_subsets(pwa)
    sort!(vp, by = v -> v[1][1])
    xp = [v[1][1] for v in vp]
    zp = [v[2] for v in vp]
    V = 2 : length(vp)
    slope = Dict(v => (zp[v] - zp[v-1]) / (xp[v] - xp[v-1]) for v in V)
    Δ = @variable(m, [v in V], lower_bound = 0, upper_bound = xp[v] - xp[v-1], base_name = "_Δ_$(counter)")
    @constraint(m, x == xp[1] + sum(Δ[v] for v in V))
    @constraint(m, z == zp[1] + sum(slope[v] * Δ[v] for v in V))
    return
end

function pwaffine_formulation(m, x, z, counter, pwa::PWAFunc{C,D}, _::Φ_Formulation) where {C,D}
    error("Φ_Formulation for piecewise affine functions only possible for 1D functions")
end

# Model the 1D piecewise affine function by having a variable that says how much each linear
# segment contributes
function pwaffine_formulation(m, x, z, counter, pwa::PWAFunc{C,1}, _::Φ_Formulation) where {C}
    vp = PiecewiseAffineApprox.vertices_by_subsets(pwa)
    sort!(vp, by = v -> v[1][1])
    xp = [v[1][1] for v in vp]
    zp = [v[2] for v in vp]
    V = 2 : length(vp)
    slope = Dict(v => (zp[v] - zp[v-1]) / (xp[v] - xp[v-1]) for v in V)
    Ṽ = 3: length(vp)
    Δ_slope = Dict(v => slope[v] - slope[v-1] for v in Ṽ)

    Φ = @variable(m, [v in Ṽ], lower_bound = 0, upper_bound = last(xp) - xp[v-1], base_name = "_Φ_$(counter)")
    @constraint(m, [v in Ṽ], Φ[v] ≥ x - xp[v-1])
    base = zp[1] + slope[2] * (x - xp[1])
    @constraint(m, z ==  base + sum(Δ_slope[v] * Φ[v] for v in Ṽ))

    return
end
