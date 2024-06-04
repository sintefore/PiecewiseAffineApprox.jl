
function _full_order_pwa(f::FunctionEvaluations{D}, optimizer) where {D}

    K = length(f.points)
    m = Model(optimizer)

    @variable(m, a[1:K, 1:D])
    @variable(m, b[1:K])

    for (l, (x, z)) in enumerate(point_vals(f))
        for k in 1:K
            @constraint(m, z ≥ sum(a[k, i] * x[i] for i in 1:D) + b[k])
        end
        @constraint(m, z ≤ sum(a[l, i] * x[i] for i in 1:D) + b[l])
    end
    optimize!(m)

    if !is_solved_and_feasible(m)
        error("No full order convex approximation exists. Function evaluations should be " *
              "from a convex function.")
    end

    aᶠ = value.(a)
    bᶠ = value.(b)
    pwl = PWLFunc([Plane(aᶠ[k,:], bᶠ[k]) for k in 1:K], Convex())
    return pwl
end

function _reduced_order_pwa(f::FunctionEvaluations{D}, p, U, optimizer) where {D}
    K = length(f.points)
    m = Model(optimizer)

    @variable(m, ared[1:p, 1:D])
    @variable(m, bred[1:p])
    @variable(m, s[1:K] ≥ 0)

    for (x, z) in point_vals(f)
        for k in 1:p
            @constraint(m, z ≥ sum(ared[k, i] * x[i] for i in 1:D) + bred[k])
        end
    end

    for (l,k) in U
        @constraint(m, f.values[l] ≤ sum(ared[k, i] * f.points[l][i] for i in 1:D) + bred[k] + s[k])
    end

    @objective(m, Min, sum(s[k] for k in 1:K))

    optimize!(m)

    obj = objective_value(m)

    aʳ = value.(ared)
    bʳ = value.(bred)

    pwl_red = PWLFunc([Plane(aʳ[k,:], bʳ[k]) for k in 1:p], Convex())

    return obj, pwl_red
end

# The hyperplane being active at the point x for a convex pwl function
function _active_plane(pwl::PWLFunc{Convex,D}, x) where {D}
    return argmax(collect(evaluate(p, x) for p ∈ pwl.planes))
end

function _update_allocation(f::FunctionEvaluations, pwl_red)
    U = [(l, _active_plane(pwl_red, x)) for (l, x) in enumerate(f.points)]
    return U
end

function _allocation_improvement(f::FunctionEvaluations, p, Uᴵ, optimizer)
    U = Uᴵ
    prev_obj = Inf64
    obj, pwl_red = _reduced_order_pwa(f, p, U, optimizer)
    while obj < prev_obj
        prev_obj = obj
        U = _update_allocation(f, pwl_red)
        obj, pwl_red = _reduced_order_pwa(f, p, U, optimizer)
    end
    return pwl_red
end

function _increase_order(f::FunctionEvaluations, pwl_red, pwl)
    s = [f.values[l] - evaluate(pwl_red, f.points[l]) for l in 1:length(f.points)]
    i = argmax(s)
    push!(pwl_red.planes, pwl.planes[i])
    U = _update_allocation(f, pwl_red)
    return U
end

# Calucate the approximation error between input data points and the piecewise
# affine approximation
function _approx_error(f::FunctionEvaluations, pwl::PWLFunc, penalty = :l1)
    err = 0.0
    for (x, z) in point_vals(f)
        v = evaluate(pwl, x)
        if penalty == :l1
            err += abs(v - z)
        elseif penalty == :l2 || penalty == :rms
            err += (v - z)^2
        elseif penalty == :max
            err = max(err, abs(v - z))
        end
    end
    if penalty == :rms
        err = err / length(f.points)
    end
    if penalty == :l2 || penalty == :rms
        err = sqrt(err)
    end
    return err
end

function _progressive_pwa(f::FunctionEvaluations, optimizer, δᵗᵒˡ = 1e-3, penalty = :max)

    # Check if the function evaluations is from a convex function
    # by solving a full order approximation problem
    pwl = _full_order_pwa(f, optimizer)

    # Start with an approximation with one plane and allocate all points to that segment
    p = 1
    U = [(l, 1) for l in 1:length(f.points)]

    # Increase the number of planes until the required tolerance is met
    pwl_red = _allocation_improvement(f, p, U, optimizer)
    while _approx_error(f, pwl_red, penalty) > δᵗᵒˡ
        U = _increase_order(f, pwl_red, pwl)
        p += 1
        pwl_red = _allocation_improvement(f, p, U, optimizer)
    end
    return pwl_red
end

function _progressive_pwa(f::FunctionEvaluations, options::ProgressiveFitting)
    return _progressive_pwa(f, options.optimizer, options.tolerance, options.pen)
end
