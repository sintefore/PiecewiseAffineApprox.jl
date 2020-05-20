defaultpenalty() = :l1
defaultseg() = 5

"""
    convex_linearization(x, z, optimizer; kwargs...)

Computes a piecewise linear function that approximates the measurements given by `x` and `z`.

# Arguments
- `method::Symbol:=fit`: the method used for approximation
- `nseg::Integer=5`: the number of segments to use 
- `pen::Symbol=:l1`: the metric used to measure deviation
"""
function convex_linearization(x::Vector, z::Vector, optimizer; kwargs...)
    @assert(length(x) == length(z))

    defaults = (;method=:fit)
    options = merge(defaults, kwargs)

    method = options.method

    if method == :fit
        return convex_linearization_fit(x, z, optimizer; kwargs...)
    elseif method == :ipol
        return convex_linearization_ipol(x, z, optimizer; kwargs...)
    else
        error("Unrecognized method $method")
    end
end

convex_linearization(x, z, optimizer; kwargs...)  = 
    convex_linearization([xx for xx in x], [zz for zz in z], optimizer; kwargs...)

function convex_linearization_fit(x::Vector, z::Vector, optimizer; kwargs...)
  
    defaults = (nseg=defaultseg(), pen=defaultpenalty(), strict=false, start_origin=false, show_res=false)
    options = merge(defaults, kwargs)
  
    N = length(x)
    𝒩 = 1:N 
    𝒦 = 1:options.nseg
    
    cᵉˢᵗ = (z[N] -z[N-1])/(x[N]-x[N-1])
    Mᵇⁱᵍ =  2 * cᵉˢᵗ * (last(x) - first(x)) - maximum(z)
    
    m = JuMP.Model()
    𝑧̂ = JuMP.@variable(m, [𝒩]) 
    𝑐 = JuMP.@variable(m, [𝒦])
    𝑑 = JuMP.@variable(m, [𝒦]) 
    𝑢 = JuMP.@variable(m, [𝒩,𝒦], Bin)

    if options.pen == :l2 
        JuMP.@objective(m, Min, sum((z[i] - 𝑧̂[i])^2 for i ∈ 𝒩))
    elseif options.pen == :max
        𝑡 = JuMP.@variable(m)
        JuMP.@objective(m, Min, 𝑡)
        for i ∈ 𝒩
            JuMP.@constraint(m,  𝑡 ≥ (z[i] - 𝑧̂[i]) )
            JuMP.@constraint(m,  𝑡 ≥ (𝑧̂[i] - z[i]) )
        end
    elseif options.pen == :l1
        𝑡 = JuMP.@variable(m, [𝒩])
        JuMP.@objective(m, Min, sum(𝑡))
        for i ∈ 𝒩
            JuMP.@constraint(m,  𝑡[i] ≥ (z[i] - 𝑧̂[i]) )
            JuMP.@constraint(m,  𝑡[i] ≥ (𝑧̂[i] - z[i]) )
        end
    else
        error("Unrecognized/unsupported penalty type $(options.pen)")
    end

    for i ∈ 𝒩, k ∈ 𝒦 
        JuMP.@constraint(m, 𝑧̂[i] ≥ 𝑐[k] * x[i] + 𝑑[k])
        JuMP.@constraint(m, 𝑧̂[i] ≤ 𝑐[k] * x[i] + 𝑑[k] + Mᵇⁱᵍ * (1-𝑢[i,k]))
    end

    if options.strict
        for i ∈ 𝒩, k ∈ 𝒦 
            JuMP.@constraint(m, z[i] ≥ 𝑐[k] * x[i] + 𝑑[k])   
        end
    end

    if options.start_origin
        JuMP.@constraint(m,𝑑[1] == 0.0)
    end

    for i ∈ 𝒩
        JuMP.@constraint(m, sum(𝑢[i,k] for k ∈ 𝒦) ≥ 1)
    end

    for k ∈ 𝒦 
        if k > 1
            JuMP.@constraint(m, 𝑐[k-1] ≤ 𝑐[k])
        end
    end


    JuMP.set_optimizer(m,optimizer)
    JuMP.optimize!(m)

    if JuMP.termination_status(m) != MOI.OPTIMAL
        error("Optimization failed")
    end

    if options.show_res
        println("Optimize succeed for $(options.pen)")
        val = JuMP.objective_value(m)
        println("Objective value = $val")
    end
    
    𝑐ᴼᵖᵗ = JuMP.value.(𝑐)
    𝑑ᴼᵖᵗ = JuMP.value.(𝑑) 

    return ConvexPWLFunction([𝑐ᴼᵖᵗ[k] for k ∈ 𝒦], [𝑑ᴼᵖᵗ[k] for k ∈ 𝒦], minimum(x), maximum(x))
end

function convex_linearization(f::Function, xmin, xmax, optimizer; kwargs...)
    @assert(xmin < xmax)

    defaults = (nsample=10, nseg=defaultseg()) 
    options = merge(defaults, kwargs)

    samples = max(options.nsample, 3*options.nseg)

    step = (xmax - xmin) / samples
    x = [i for i in xmin:step:xmax]
    y = [f(xx) for xx in x]
    return convex_linearization(x, y, optimizer; kwargs...)
end

function concave_linearization(x, z, optimizer; kwargs...)
    return ConcavePWLFunction(convex_linearization(x, -z, optimizer; kwargs...))
end

function concave_linearization(f::Function, xmin, xmax, optimizer; kwargs...)
    return ConcavePWLFunction(convex_linearization(x -> -f(x), xmin, xmax, optimizer; kwargs...))
end


function convexify(pwl::PWLFunction, optimizer)
    N = length(pwl.x)
    𝒩 = 1: N
    x = pwl.x
    z = pwl.z

    m = Model()

    δ⁺ = @variable(m, [𝒩], lower_bound=0)
    δ⁻ = @variable(m, [𝒩], lower_bound=0)

    @objective(m, Min, sum(δ⁺[i] + δ⁻[i] for i ∈ 𝒩))
    for i=2:N-1
        @constraint(m, (z[i] + δ⁺[i] - δ⁻[i] - z[i-1] - δ⁺[i-1] + δ⁻[i-1]) / 
            (x[i] -x[i-1]) <= (z[i+1] + δ⁺[i+1] - δ⁻[i+1] - z[i] - δ⁺[i] + δ⁻[i]) / (x[i+1] -x[i]))
    end

    set_optimizer(m,optimizer)
    optimize!(m)
    if termination_status(m) != MOI.OPTIMAL
        error("Optimization failed")
    end

    dp = value.(δ⁺)
    dn = value.(δ⁻)

    return ConvexPWLFunction(x, [z[i] + dp[i] - dn[i] for i ∈ 𝒩]) 
end

convexify(x, z, optimizer) =
    convexify(PWLFunction(x,z), optimizer)

function interpolatepw(x, z, optimizer; kwargs...)
    @assert(length(x) == length(z))

    defaults = (nseg=defaultseg(), pen=defaultpenalty())
    options = merge(defaults, kwargs)
   
    N = length(x)
    𝒩 = 1:N 

    # Find slopes
    c = [(z[j] -z[i]) / (x[j]-x[i])  for i ∈ 𝒩, j ∈ 𝒩]

    # Calculate penalties 
    if options.pen == :l1
        p = [(i < j ? sum(abs(c[i,j] * (x[k] - x[i]) + z[i] - z[k]) for k ∈ i:j) : 0) for i ∈ 𝒩, j ∈ 𝒩]
    elseif options.pen == :l2
        p = [(i < j ? sum((c[i,j] * (x[k] - x[i]) + z[i] - z[k])^2 for k ∈ i:j) : 0) for i ∈ 𝒩, j ∈ 𝒩]
    elseif options.pen == :max
        p = [(i < j ? maximum(abs(c[i,j] * (x[k] - x[i]) + z[i] - z[k]) for k ∈ i:j) : 0) for i ∈ 𝒩, j ∈ 𝒩]
    else
        error("Unrecognized/unsupported penalty $(options.pen)")
    end
    
    m = Model()
    @variable(m, 𝑢[𝒩,𝒩], Bin)

    # Minimize total penalty
    @objective(m, Min, sum(p[i,j] * 𝑢[i,j] for i ∈ 𝒩, j ∈ 𝒩))

    # Number of line segments in interpolant
    @constraint(m, sum(𝑢[i,j] for i ∈ 𝒩, j ∈ 𝒩) == options.nseg )

    # Only forward segments allowed
    for i ∈ 𝒩, j ∈ 𝒩 
        if i >= j 
            @constraint(m, 𝑢[i,j] == 0)
        end
    end 

    # Path modelling
    for j ∈ 𝒩
        @constraint(m, (j > 1 ? sum(𝑢[i,j] for i ∈ 𝒩) : 1) == (j < N ? sum(𝑢[j,i] for i ∈ 𝒩 ) : 1))
    end

    set_optimizer(m,optimizer)
    optimize!(m)
    if termination_status(m) != MOI.OPTIMAL
        error("Optimization failed")
    end
    
    𝑢ᴼᵖᵗ = value.(m[:𝑢])
    xᴼᵖᵗ = collect(x[i] for i ∈ 𝒩, j ∈ 𝒩 if 𝑢ᴼᵖᵗ[i,j] == 1)
    push!(xᴼᵖᵗ, x[N])
    zᴼᵖᵗ = collect(z[i] for i ∈ 𝒩, j ∈ 𝒩 if 𝑢ᴼᵖᵗ[i,j] == 1)
    push!(zᴼᵖᵗ, z[N])

    return PWLFunction(xᴼᵖᵗ, zᴼᵖᵗ)
end

convex_linearization_ipol(x, z, optimizer; kwargs...) = 
    convexify(interpolatepw(x, z, optimizer; kwargs...),optimizer)
