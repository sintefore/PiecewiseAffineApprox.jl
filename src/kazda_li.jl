#=
This file holds an implementation of a progressive fitting heuristic
inspired by the algorithm in the paper "A linear programming approach to
difference-of-convex piecewise linear approximation" by Kazda and Li (2024).

The current implementation requires that the function evaluations are
samples of a convex function, i.e. it is not robust with regard to noisy
and non-convex data.

The main algorithm consists of the following steps:
1. While error < tolerance
    a. Introduce an extra plane in the pwa approximation
    b. Find local optimum for allocating points to plane (using linear programming)
=#

function convexify(f::FunctionEvaluations{D}, optimizer, pen = :l1) where {D}
    K = length(f)

    m = Model(optimizer)
    @variable(m, a[1:K, 1:D])
    @variable(m, b[1:K])
    @variable(m, s⁺[1:K] ≥ 0)
    @variable(m, s⁻[1:K] ≥ 0)

    for (l, (x, z)) ∈ enumerate(f)
        for k ∈ 1:K
            @constraint(m, z + s⁺[l] - s⁻[l] ≥ dot(a[k, :], x) + b[k])
        end
        @constraint(m, z + s⁺[l] - s⁻[l] ≤ dot(a[l, :], x) + b[l])
    end
    obj = AffExpr()
    if pen == :l2 || pen == :rms
        obj = sum(s⁺[l]^2 + s⁻[l]^2 for l ∈ 1:K)
    elseif pen == :max
        @variable(m, smax)
        obj = smax
        for l ∈ 1:K
            @constraint(m, smax ≥ s⁺[l] + s⁻[l])
        end
    elseif pen == :l1
        obj = sum(s⁺[l] + s⁻[l] for l ∈ 1:K)
    else
        error("Unrecognized/unsupported penalty type $(penalty)")
    end
    @objective(m, Min, obj)

    optimize!(m)

    if !is_solved_and_feasible(m)
        JuMP.write_to_file(m, "convexify.lp")
        error("Unable to find a feasible solution")
    end

    obj_val = objective_value(m)
    if obj_val == 0
        @info "Data points are convex"
        return f
    end

    @info "Data points are not convex, deviation = $obj_val"
    values = [f.values[l] + value(s⁺[l]) - value(s⁻[l]) for l ∈ 1:K]
    return FunctionEvaluations(f.points, values)
end

# Create a full order convex approximation for the given data points,
# i.e. using a separate plane for each data point.
function _full_order_pwa(f::FunctionEvaluations{D}, optimizer) where {D}
    K = length(f)

    m = Model(optimizer)
    @variable(m, a[1:K, 1:D])
    @variable(m, b[1:K])
    for (l, (x, z)) ∈ enumerate(f)
        for k ∈ 1:K
            @constraint(m, z ≥ dot(a[k, :], x) + b[k])
        end
        @constraint(m, z ≤ dot(a[l, :], x) + b[l])
    end
    optimize!(m)

    if !is_solved_and_feasible(m)
        error(
            "No full order convex approximation exists. Function evaluations should be " *
            "from a convex function.",
        )
    end

    aᶠ = value.(a)
    bᶠ = value.(b)
    pwa = PWAFunc([Plane(aᶠ[k, :], bᶠ[k]) for k ∈ 1:K], Convex())
    return pwa
end

# Create a pwa approximation using `p` planes where `U` provides a mapping from
# data points indices to the indices of planes. The planes are required to lie below
# all data points. Points are allowed to lie above their associated plane with
# a penalty that is minimized in the optimization.
function _reduced_order_pwa(
    f::FunctionEvaluations{D},
    p::Int,
    U,
    optimizer,
) where {D}
    K = length(f)

    m = Model(optimizer)
    @variable(m, ared[1:p, 1:D])
    @variable(m, bred[1:p])
    @variable(m, s[1:K] ≥ 0)
    for (x, z) ∈ f
        for k ∈ 1:p
            @constraint(m, z ≥ dot(ared[k, :], x) + bred[k])
        end
    end
    for (l, k) ∈ U
        z = f.values[l]
        x = f.points[l]
        @constraint(m, z ≤ dot(ared[k, :], x) + bred[k] + s[l])
    end
    @objective(m, Min, sum(s[l] for l ∈ 1:K))
    optimize!(m)
    obj = objective_value(m)

    aʳ = value.(ared)
    bʳ = value.(bred)
    pwa_red = PWAFunc([Plane(aʳ[k, :], bʳ[k]) for k ∈ 1:p], Convex())

    return obj, pwa_red
end

# The hyperplane being active at the point x for a convex pwa function
function _active_plane(pwl::PWAFunc{Convex,D}, x) where {D}
    return argmax(collect(evaluate(p, x) for p ∈ pwl.planes))
end

# Map data points to the active plane for the piecewise approximation
function _update_allocation(f::FunctionEvaluations, pwa_red)
    U = [(l, _active_plane(pwa_red, f.points[l])) for l ∈ 1:length(f)]
    return U
end

# Do a local improvement of point allocation to planes until error
# is not further reduced
function _allocation_improvement(f::FunctionEvaluations, p, Uᴵ, optimizer)
    U = Uᴵ
    prev_obj = Inf64
    obj, pwa_red = _reduced_order_pwa(f, p, U, optimizer)
    while obj < prev_obj
        prev_obj = obj
        U = _update_allocation(f, pwa_red)
        obj, pwa_red = _reduced_order_pwa(f, p, U, optimizer)
    end
    return pwa_red
end

# Add an extra plane to the pwa by identifying the point lying
# further away and use the corresponding plane in the full order pwa
# to update the allocation of points to planes.
function _increase_order(f::FunctionEvaluations, pwa_red, pwa, used)
    s = [z - evaluate(pwa_red, x) for (x, z) ∈ f]
    imax = 0
    smax = -1
    for i ∈ 1:length(f)
        if s[i] > smax && !(i in used)
            imax = i
            smax = s[i]
        end
    end
    push!(pwa_red.planes, pwa.planes[imax])
    U = _update_allocation(f, pwa_red)
    return U, vcat(used, imax)
end

# Calculate the approximation error between input data points and the piecewise
# affine approximation in different metrics
function _approx_error(f::FunctionEvaluations, pwa::PWAFunc, penalty = :l1)
    err = 0.0
    for (x, z) ∈ f
        v = evaluate(pwa, x)
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

# Main algorithm
function _progressive_pwa(
    f::FunctionEvaluations,
    optimizer,
    δᵗᵒˡ = 1e-3,
    penalty = :max,
)

    # Find the full order convex approximation, will error if the points are
    # not convex
    pwa = _full_order_pwa(f, optimizer)

    # Start with an approximation with one plane and allocate all points to that segment
    p = 1
    U = [(l, 1) for l ∈ 1:length(f)]

    # Increase the number of planes until the required tolerance is met
    pwa_red = _allocation_improvement(f, p, U, optimizer)
    used = []
    while _approx_error(f, pwa_red, penalty) > δᵗᵒˡ && p < length(f.points)
        U, used = _increase_order(f, pwa_red, pwa, used)
        p += 1
        pwa_red = _allocation_improvement(f, p, U, optimizer)
    end
    @info "Fitting finished, error = $(round(_approx_error(f, pwa_red, penalty); digits = 3)), p = $p"
    return pwa_red
end

# Map algorithm structure to correct parameters
function _progressive_pwa(f::FunctionEvaluations, options::ProgressiveFitting)
    return _progressive_pwa(
        f,
        options.optimizer,
        options.tolerance,
        options.pen,
    )
end
