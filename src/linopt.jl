
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
    pwlinear(m::JuMP.Model, x::Tuple, pwl::PWLFunc{C,D}; z=nothing, kwargs...) where {C,D}

Add constraints to JuMP-model `m` for JuMP-variable `z` as a  
piecewise linear function/approximation `pwl` of JuMP-variables `x`    
"""
function pwlinear(m::JuMP.Model, x, pwl::PWLFunc{C,D}; z = nothing) where {C,D}
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter

    if isnothing(z)
        # @warn "NB: Skipping bounds for now"
        z = JuMP.@variable(m, base_name = "z_$(counter)")
    end
    for (k, p) ∈ enumerate(pwl.planes)
        con = constr(C, m, z, p, x)
        JuMP.set_name(con, "pwl_$(counter)_$(k)")
    end
    return z
end

function pwlinear(
    m::JuMP.Model,
    x,
    fevals::FunctionEvaluations,
    curvature::Curvature,
    a::Algorithm;
    z = nothing,
)
    return pwlinear(m, x, approx(fevals, curvature, a); z = z)
end
