defaultpenalty() = :l1
defaultpenalty2D() = :l2
defaultseg() = 5
defaultplanes() = 4

"""
    convex_linearization(x, z, optimizer; kwargs...)

Computes a piecewise linear function that approximates the measurements given by `x` and `z`.

# Arguments
- `method::Symbol:=fit`: the method used for approximation
- `dimensions::Integer:=2`: the number of dimensions of the function domain
- `nseg::Integer=5`: the number of segments to use 
- `nplanes::Integer=4`: the number of planes to use in 2D PWL functions
- `strict::Symbol=:none`: defines it is a general approximation, or an overestimation or underestimation
- `pen::Symbol=:l1`: the metric used to measure deviation


"""
function convex_linearization(x, z, optimizer; kwargs...)
    defaults = (;method=:fit, dimensions=1)
    options = merge(defaults, kwargs)

    method = options.method
    dimensions = options.dimensions

    if dimensions == 1
        @assert(length(x) == length(z))
        if method == :fit
            return convex_linearization_fit(x, z, optimizer; kwargs...)
        elseif method == :ipol
            return convex_linearization_ipol(x, z, optimizer; kwargs...)
        else
            error("Unrecognized method $method")
        end
    elseif dimensions >= 2
        if method == :fit
            return convex_ND_linearization_fit(x, z, optimizer; kwargs...)
        end    
    else
        error("Unrecognized number of dimensions $dimensions")
    end    
end

function conv_linear_big_M(x, z)
    N = length(x)
    cᵉˢᵗ = (z[N] -z[N-1])/(x[N]-x[N-1])
    return  2 * cᵉˢᵗ * (last(x) - first(x)) - maximum(z)
end

#convex_linearization(x, z, optimizer; kwargs...)  = 
#    convex_linearization([xx for xx in x], [zz for zz in z], optimizer; kwargs...)

function convex_linearization_fit(x::Vector, z::Vector, optimizer; kwargs...)
  
    defaults = (nseg=defaultseg(), pen=defaultpenalty(), strict=false, start_origin=false, show_res=false)
    options = merge(defaults, kwargs)
  
    N = length(x)
    𝒩 = 1:N 
    𝒦 = 1:options.nseg
    
    Mᵇⁱᵍ =  conv_linear_big_M(x,z)
    
    m = Model()
    @variable(m, 𝑧̂[𝒩]) 
    @variable(m, 𝑐[𝒦])
    @variable(m, 𝑑[𝒦]) 
    @variable(m, 𝑢[𝒩,𝒦], Bin)

    if options.pen == :l2 
        @objective(m, Min, sum((z[i] - 𝑧̂[i])^2 for i ∈ 𝒩))
    elseif options.pen == :max
        𝑡 = @variable(m)
        @objective(m, Min, 𝑡)
        for i ∈ 𝒩
            @constraint(m,  𝑡 ≥ (z[i] - 𝑧̂[i]) )
            @constraint(m,  𝑡 ≥ (𝑧̂[i] - z[i]) )
        end
    elseif options.pen == :l1
        𝑡 = @variable(m, [𝒩])
        @objective(m, Min, sum(𝑡))
        for i ∈ 𝒩
            @constraint(m,  𝑡[i] ≥ (z[i] - 𝑧̂[i]) )
            @constraint(m,  𝑡[i] ≥ (𝑧̂[i] - z[i]) )
        end
    else
        error("Unrecognized/unsupported penalty type $(options.pen)")
    end

    for i ∈ 𝒩, k ∈ 𝒦 
        @constraint(m, 𝑧̂[i] ≥ 𝑐[k] * x[i] + 𝑑[k])
        @constraint(m, 𝑧̂[i] ≤ 𝑐[k] * x[i] + 𝑑[k] + Mᵇⁱᵍ * (1-𝑢[i,k]))
    end

    if options.strict
        for i ∈ 𝒩, k ∈ 𝒦 
            @constraint(m, z[i] ≥ 𝑐[k] * x[i] + 𝑑[k])   
        end
    end

    if options.start_origin
        @constraint(m,𝑑[1] == 0.0)
    end

    for i ∈ 𝒩
        @constraint(m, sum(𝑢[i,k] for k ∈ 𝒦) ≥ 1)
    end

    for k ∈ 𝒦 
        if k > 1
            @constraint(m, 𝑐[k-1] ≤ 𝑐[k])
        end
    end


    set_optimizer(m,optimizer)
    optimize!(m)

    if termination_status(m) != MOI.OPTIMAL
        error("Optimization failed")
    end

    if options.show_res
        println("Optimize succeed for $(options.pen)")
        val = objective_value(m)
        println("Objective value = $val")
    end
    
    𝑐ᴼᵖᵗ = value.(𝑐)
    𝑑ᴼᵖᵗ = value.(𝑑) 

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
    if (size(x,2) < 2) && (length(first(x)) < 2)
        return ConcavePWLFunction(convex_linearization(x, -z, optimizer; kwargs...))
    else
        return ConcavePWLFunctionND(convex_linearization(x, -z, optimizer; kwargs...))
    end
end

function concave_linearization(f::Function, xmin, xmax, optimizer; kwargs...)
    if (size(x,2) < 2) && (length(first(X)) < 2)
        return ConcavePWLFunction(convex_linearization(x -> -f(x), xmin, xmax, optimizer; kwargs...))
    else
        return ConcavePWLFunctionND(convex_linearization(x -> -f(x), xmin, xmax, optimizer; kwargs...))
    end
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

function mat2tuples(x::Matrix)
    return collect(Tuple(x'[:,i]) for i in 1:size(x',2))
end
# mat2tuples2(x::Matrix) = ((Tuple(x[i,:]) for i in 1:size(x,1)))

function tuples2mat(𝒫::Vector{Tuple{Float64, Float64}})
   return reduce(hcat, getindex.(𝒫,i) for i in eachindex(𝒫[1]))
end

function convex_ND_linearization_fit(x::Matrix{Float64}, z, optimizer; kwargs...)
    return convex_ND_linearization_fit(mat2tuples(x), z, optimizer; kwargs...)
end

function convND_linear_big_M(𝒫::Vector{Tuple{Float64, Float64}}, z)
    return convND_linear_big_M(tuples2mat(𝒫),z)
end

function convND_linear_big_M(x::Matrix{Float64}, z)
    ## TODO: calculate a tighter value for the big-M    
    return 2*maximum(z)
end

function convex_ND_linearization_fit(𝒫, z, optimizer; kwargs...)

    defaults = (nsegs=defaultseg(), nplanes=defaultplanes(), pen=defaultpenalty2D(), strict=:none, show_res=false)
    options = merge(defaults, kwargs)

    zᵖ = Dict(zip(𝒫, z))
    𝒦 = 1:options.nplanes
    ℐₚ = 1:length(𝒫[1])    

    Mᵇⁱᵍ = convND_linear_big_M(𝒫, z) 

    m = Model()
    @variable(m, 𝑧̂[𝒫])
    @variable(m, a[ℐₚ, 𝒦])
    @variable(m, b[𝒦])

    @variable(m, 𝑢[𝒫, 𝒦], Bin)

    if options.pen == :l2 
        @objective(m, Min, sum((zᵖ[p] - 𝑧̂[p])^2 for p ∈ 𝒫))
    elseif options.pen == :max
        𝑡 = @variable(m)
        @objective(m, Min, 𝑡)
        for p ∈ 𝒫
            @constraint(m,  𝑡 ≥ (zᵖ[p] - 𝑧̂[p]) )
            @constraint(m,  𝑡 ≥ (𝑧̂[p] - zᵖ[p]) )
        end
    elseif options.pen == :l1
        𝑡 = @variable(m, [𝒫])
        @objective(m, Min, sum(𝑡))
        for p ∈ 𝒫
            @constraint(m,  𝑡[p] ≥ (zᵖ[p] - 𝑧̂[p]) )
            @constraint(m,  𝑡[p] ≥ (𝑧̂[p] - zᵖ[p]) )
        end
    else
        error("Unrecognized/unsupported penalty type $(options.pen)")
    end
     
    for p ∈ 𝒫, k ∈ 𝒦         
        @constraint(m, 𝑧̂[p] ≥ sum(a[j,k] * p[j] for j in ℐₚ) + b[k])
        @constraint(m, 𝑧̂[p] ≤ sum(a[j,k] * p[j] for j in ℐₚ) + b[k] + Mᵇⁱᵍ * (1-𝑢[p,k]))                
    end

    if options.strict == :above
        for p ∈ 𝒫, k ∈ 𝒦 
            @constraint(m, zᵖ[p] ≥ sum(a[j,k] * p[j] for j in ℐₚ) + b[k]) 
        end
    elseif options.strict == :below
        for p ∈ 𝒫, k ∈ 𝒦 
            @constraint(m, zᵖ[p] ≤ sum(a[j,k] * p[j] for j in ℐₚ) + b[k]) 
        end
    end
    
    for p ∈ 𝒫
        @constraint(m, sum(𝑢[p,k] for k ∈ 𝒦) ≥ 1)
    end    
    
    set_optimizer(m,optimizer)
    optimize!(m)

    if termination_status(m) != MOI.OPTIMAL
        error("Optimization failed")
    end

    if options.show_res
        println("Optimize succeed for $(options.pen)")
        val = objective_value(m)
        println("Objective value = $val")
    end   
    
    aᴼᵖᵗ = value.(a)
    bᴼᵖᵗ = value.(b)    
    
    return ConvexPWLFunctionND(collect([Tuple(aᴼᵖᵗ.data[:,k]) for k ∈ 𝒦]),  [bᴼᵖᵗ[k] for k ∈ 𝒦])
    ##TODO: how to recover the data points from the coefficients?  Check package Polyhedra.    
end    

