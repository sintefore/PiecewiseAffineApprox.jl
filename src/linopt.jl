
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

function convex_pwlinear(m::JuMP.Model, x::VarOrAff, pwl::ConvexPWLFunction; z=nothing)
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter
    
    if isnothing(z)
         z = JuMP.@variable(m, lower_bound=minimum(pwl.z), upper_bound=maximum(pwl.z), base_name="z_$(counter)") 
    end
    for k=1:length(pwl.c)
        con = JuMP.@constraint(m, z ≥ pwl.d[k] + pwl.c[k] * x)
        JuMP.set_name(con, "pwl_$(counter)_$(k)")
    end
    return z
end

convex_pwlinear(m::JuMP.Model, x::VarOrAff, xd::Vector, zd::Vector, optimizer; z=nothing, kwargs...) =
    convex_pwlinear(m, x, convex_linearization(xd, zd, optimizer; kwargs...), z=z)

convex_pwlinear(m::JuMP.Model, x::VarOrAff, f::Function, xmin, xmax, optimizer; z=nothing, kwargs...) =
    convex_pwlinear(m, x, convex_linearization(f,xmin,xmax, optimizer; kwargs...), z=z)

function concave_pwlinear(m::JuMP.Model, x::VarOrAff, pwl::ConcavePWLFunction; z=nothing)
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter
    
    if isnothing(z)
        z = JuMP.@variable(m, lower_bound=minimum(pwl.z), upper_bound=maximum(pwl.z), base_name="z_$counter")  
    end
    for k=1:length(pwl.c)
        JuMP.@constraint(m, z ≤ pwl.d[k] + pwl.c[k] * x, name="pwl_$counter_$k")
    end
    return z
end

concave_pwlinear(m::JuMP.Model, x::VarOrAff, xd::Vector, zd::Vector, optimizer; z=nothing, kwargs...) =
    concave_pwlinear(m, x, convex_linearization(xd, fd, optimizer, kwargs...), z=z)

concave_pwlinear(m::JuMP.Model, x::VarOrAff, f::Function, xmin, xmax, optimizer; z=nothing, kwargs...) =
    concave_pwlinear(m, x, convex_linearization(f,xmin,xmax, optimizer,kwargs...), z=z)

    
convex_pwlinear(m::JuMP.Model, x::Tuple, xd::Vector, zd::Vector, optimizer; z=nothing, kwargs...) =
    convex_pwlinear(m, x, convex_linearization(xd, zd, optimizer; kwargs...), z=z)

convex_pwlinear(m::JuMP.Model, x::Tuple, xd::Matrix, zd::Vector, optimizer; z=nothing, kwargs...) =
    convex_pwlinear(m, x, convex_linearization(xd, zd, optimizer; kwargs...), z=z)


constr(::Type{Convex},m,z,p,x) = JuMP.@constraint(m, z ≥ dot(p.α, x) + p.β)
constr(::Type{Concave},m,z,p,x) = JuMP.@constraint(m, z ≤ dot(-1 .* p.α, x) - p.β)
"""
    pwlinear(m::JuMP.Model, x::Tuple, pwl::PWLFunc{C,D}; z=nothing, kwargs...) where {C,D}

Add constraints to JuMP-model `m` for JuMP-variable `z` as a  
piecewise linear function/approximation `pwl` of JuMP-variables `x`    
"""
function pwlinear(m::JuMP.Model, x, pwl::PWLFunc{C,D}; z=nothing, kwargs...) where {C,D}
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter

    if isnothing(z)
        @warn "NB: Skipping bounds for now"
        z = JuMP.@variable(m, base_name="z_$(counter)") 
    end
    for (k,p) in enumerate(pwl.planes)
        con = constr(C,m,z,p,x)
        JuMP.set_name(con, "pwl_$(counter)_$(k)")
    end
    return z
end
pwlinear(m::JuMP.Model, x, fevals::FunctionEvaluations,curvature::Curvature,a::Algorithm;kwargs...) = pwlinear(m,x,approx(fevals,curvature,a;kwargs...);kwargs...)

@deprecate convex_pwlinear pwlinear
function convex_pwlinear(m::JuMP.Model, x::Tuple, pwl::PWLFunc{C,D}; z=nothing) where {C<:Convex,D}
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter
    
    function minz(pwl)
        minimum(p.β for p ∈ pwl.planes)
    end
    function maxz(pwl)
        maximum(p.β for p ∈ pwl.planes)
    end

    if isnothing(z)
        #  z = JuMP.@variable(m, lower_bound=minz(pwl), upper_bound=maxz(pwl), base_name="z_$(counter)") 
        # TODO: Fix bounds
        z = JuMP.@variable(m, base_name="z_$(counter)") 
    end
    for (k,p) in enumerate(pwl.planes)
        con = JuMP.@constraint(m, z ≥ dot(p.α, x) + p.β)
        JuMP.set_name(con, "pwl_$(counter)_$(k)")
    end
    return z
end

@deprecate concave_pwlinear pwlinear
function concave_pwlinear(m::JuMP.Model, x::Tuple, pwl::PWLFunc{C,D}; z=nothing) where {C<:Concave, D}
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter

    # pwl = pwl_concave.pwl   
    sign = -1
    
    if isnothing(z)
        # z = JuMP.@variable(m, lower_bound=minimum(pwl.z), upper_bound=maximum(pwl.z), base_name="z_$counter")  
        z = JuMP.@variable(m, base_name="z_$counter")  
    end
    for (k,p) ∈ enumerate(pwl.planes)
        con = JuMP.@constraint(m, z ≤ dot(p.α.*sign, x) + p.β.*sign)
        JuMP.set_name(con, "pwl_$(counter)_$(k)")
    end    
    
    return z
end

concave_pwlinear(m::JuMP.Model, x::Tuple, xd::Vector, zd::Vector, optimizer; z=nothing, kwargs...) =
    concave_pwlinear(m, x, concave_linearization(xd, fd, optimizer, kwargs...), z=z)

concave_pwlinear(m::JuMP.Model, x::Tuple, f::Function, xmin, xmax, optimizer; z=nothing, kwargs...) =
    concave_pwlinear(m, x, concave_linearization(f,xmin,xmax, optimizer,kwargs...), z=z)

concave_pwlinear(m::JuMP.Model, x::Tuple, xd::Matrix, zd::Vector, optimizer; z=nothing, kwargs...) =
concave_pwlinear(m, x, concave_linearization(xd, zd, optimizer; kwargs...), z=z)
    # concave_pwlinear(m, x, approx(FunctionEvaluations(PWL.mat2tuples(xd), zd), Concave(), Optimized(); merge((;optimizer=optimizer), kwargs...), z=z)

#concave_pwlinear(m::JuMP.Model, x::Tuple, xd::Vector, zd::Vector, optimizer; z=nothing, kwargs...) =
#    concave_pwlinear(m, x, concave_linearization(xd, zd, optimizer; kwargs...), z=z)    