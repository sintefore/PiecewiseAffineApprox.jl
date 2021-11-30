
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

function convex_pwlinear(m::JuMP.Model, x::Tuple, pwl::ConvexPWLFunctionND; z=nothing)
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter
    
    if isnothing(z)
         z = JuMP.@variable(m, lower_bound=minimum(pwl.z), upper_bound=maximum(pwl.z), base_name="z_$(counter)") 
    end
    for k=1:length(pwl.a)
        con = JuMP.@constraint(m, z ≥ dot(pwl.a[k], x) + pwl.b[k])
        JuMP.set_name(con, "pwl_$(counter)_$(k)")
    end
    return z
end

function concave_pwlinear(m::JuMP.Model, x::Tuple, pwl_concave::ConcavePWLFunctionND; z=nothing)
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
    m.ext[:PWL].counter = counter

    pwl = pwl_concave.pwl   
    sign = -1
    
    if isnothing(z)
        z = JuMP.@variable(m, lower_bound=minimum(pwl.z), upper_bound=maximum(pwl.z), base_name="z_$counter")  
    end
    for k=1:length(pwl.a)
        con = JuMP.@constraint(m, z ≤ dot(pwl.a[k].*sign, x) + pwl.b[k].*sign)
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

#concave_pwlinear(m::JuMP.Model, x::Tuple, xd::Vector, zd::Vector, optimizer; z=nothing, kwargs...) =
#    concave_pwlinear(m, x, concave_linearization(xd, zd, optimizer; kwargs...), z=z)    