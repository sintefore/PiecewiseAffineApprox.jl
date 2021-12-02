
function approx_error(x, z, pwl::PWLFunc, penalty = :l1)
    err = 0.0
    for i in 1:size(x,2)
        x̄ = x[:,i]
        v = evaluate(pwl, x̄)

        if penalty == :l1
            err += abs(v - z[i])
        elseif penalty == :l2 || penalty == :rms
            err += (v - z[i])^2
        elseif penalty == :max
            err = max(err, abs(v - z[i]))
        end
    end

    if penalty == :rms
        err = err / size(x,2)
    end

    if penalty == :l2 || penalty == :rms
        err = sqrt(err)
    end

    return err
end


function convex_linearization_mb(x::Matrix, z::Vector; kwargs...)

    defaults = (nseg=defaultseg(), pen=defaultpenalty(), trials=20, itlim=50, strict=:none)
    options = merge(defaults, kwargs)

    Nᵗʳ = options.trials     # Number of trials
    lᵐᵃˣ = options.itlim    # Iteration limit
    K = options.nseg
    penalty = options.pen
    strict = options.strict
    optimizer = options.optimizer

    pwl_best = nothing
    min_error = Inf

    @info "Starting heuristic search "
    for i in 1:Nᵗʳ 
        @info "  Trial $i"     
        #@info "Creating random partition"
        𝒫 = random_partition(x, K)
      
        #@info "Searching for optimal partition"
        pwl = optimal_partition(x, z, 𝒫, lᵐᵃˣ, penalty, optimizer, strict)
      
        e = approx_error(x, z, pwl, penalty)
        @info "  Approximation error = $e ($penalty)"
        if e < min_error
            min_error = e
            pwl_best = pwl
        end
    end
    @info "Terminating search - best approximation error = $(min_error) ($penalty)"
    return pwl_best
end

dist(x,y) = sqrt(sum((x[i] - y[i])^2 for i in 1:length(x)))

function random_partition(x, K)

    N = size(x)[2]

    μ = vec(mean(x, dims=2))
    σ² = cov(x, dims=2)
    n = MvNormal(μ, σ²)
    p = rand(n, K)
    
    𝒫 = Dict(j => [] for j in 1:K)
    for i in 1:N
        # Find the nearest point amongst the p's
        jmin = argmin([dist(x[:,i], p[:,j]) for j in 1:K])
        push!(𝒫[jmin], i)  
    end
    return 𝒫
end

function optimal_partition(x, z, 𝒫, lᵐᵃˣ, penalty, optimizer, strict)

    D = size(x,1)
    pwl = nothing
    for it in 1:lᵐᵃˣ
        pwl = PWLFunc{Convex,D}()
        for j in 1:length(𝒫)
            if length(𝒫[j]) > 0 
                x̄ = x[:, 𝒫[j] ]
                z̄ = z[𝒫[j]]
                a,b = local_fit(x̄, z̄, penalty, optimizer, strict)
                addplane!(pwl, a, b)
            end
        end
        𝒫ⁿᵉʷ = update_partition(x, pwl)
        if  𝒫ⁿᵉʷ == 𝒫
            break
        end
        𝒫 = 𝒫ⁿᵉʷ
    end

    return pwl
end

function local_fit(x̄, z̄, penalty, optimizer, strict)

    M, N = size(x̄)
    
    # Create an optimization model to find the best a and b such that  ax + b ≈ y
    m = Model()

    @variable(m, a[1:M])
    @variable(m, b)
    @variable(m, ẑ[1:N]) 

    for i in 1:N
        @constraint(m, ẑ[i] == sum(a[j] * x̄[j,i] for j in 1:M) + b)
    end

    if strict == :above
        for i in 1:N
            @constraint(m, ẑ[i] ≥ z̄[i])
        end
    elseif strict == :below
        for i in 1:N
            @constraint(m, ẑ[i] ≤ z̄[i])
        end
    end

    obj = AffExpr()
    if penalty == :l2 || penalty ==:rms
        obj = sum((z̄[i] - ẑ[i])^2 for i ∈ 1:N)
    elseif penalty == :max
        @variable(m, t)
        obj = t
        for i ∈ 1:N
            @constraint(m,  t ≥ (z̄[i] - ẑ[i]) )
            @constraint(m,  t ≥ (ẑ[i] - z̄[i]) )
        end
    elseif penalty == :l1
        @variable(m, t[1:N])
        obj = sum(t[i] for i in 1:N)
        for i ∈ 1:N
            @constraint(m, t[i] ≥ (z̄[i] - ẑ[i]) )
            @constraint(m, t[i] ≥ (ẑ[i] - z̄[i]) )
        end
    else
        error("Unrecognized/unsupported penalty type $(options.pen)")
    end

    # TODO: consider adding regularization term

    @objective(m, Min, obj)

    set_optimizer(m,optimizer)
    set_silent(m)
    optimize!(m)

    if termination_status(m) != MOI.OPTIMAL
        error("Optimization failed")
    end

    ā = JuMP.value.(a)
    b̄ = JuMP.value.(b)

    return ā, b̄
end

function update_partition(x, pwl)
    𝒫 = Dict(j => [] for j in 1:nplanes(pwl))
    for i in 1:size(x,2)
       push!(𝒫[active(pwl, x[:,i])], i) 
    end
    return 𝒫
end

