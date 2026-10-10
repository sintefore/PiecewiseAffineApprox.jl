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

"""
    enforce_curvature(f::FunctionEvaluations, curvature::Curvature, optimizer, metric = :l1)

Create a slightly perturbed version of the function values
to ensure that the data points can be interpolated by a convex/concave
piecewise affine function.

Enforcing the curvature is performed using a linear optimization problem
that adjust the function values and tries to minimize the total deviation.
The total deviation can be measured in different metrics specified
by the `metric` parameter: `:l1` (sum of absolute changes), `:l2` (sum of
squared changes, solved as a quadratic program) or `:max` (largest change).
This function can be useful as a pre-step for running the full-order and progressive
fitting heuristics that require data points that are convex/concave in this sense.
"""
function enforce_curvature(
    f::FunctionEvaluations,
    c::Convex,
    optimizer,
    metric = :l1,
)
    return convexify(f, optimizer, metric)
end

function enforce_curvature(
    f::FunctionEvaluations,
    c::Concave,
    optimizer,
    metric = :l1,
)
    fneg = FunctionEvaluations(f.points, -f.values)
    fconvex = convexify(fneg, optimizer, metric, "concave")
    return FunctionEvaluations(f.points, -fconvex.values)
end

function convexify(
    f::FunctionEvaluations{D},
    optimizer,
    metric = :l1,
    curv_string = "convex",
) where {D}
    if metric ∉ (:l1, :l2, :rms, :max)
        throw(
            ArgumentError(
                "Unsupported metric :$metric for enforcing curvature, " *
                "use :l1, :l2 or :max",
            ),
        )
    end
    K = length(f)

    m = Model(optimizer)
    set_silent(m)
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
    if metric == :l2 || metric == :rms
        obj = sum(s⁺[l]^2 + s⁻[l]^2 for l ∈ 1:K)
    elseif metric == :max
        @variable(m, smax)
        obj = smax
        for l ∈ 1:K
            @constraint(m, smax ≥ s⁺[l] + s⁻[l])
        end
    else
        obj = sum(s⁺[l] + s⁻[l] for l ∈ 1:K)
    end
    @objective(m, Min, obj)

    optimize!(m)

    if !is_solved_and_feasible(m; allow_local = true)
        error(
            "Unable to enforce $(curv_string) curvature: " *
            "$(termination_status(m)) ($(raw_status(m)))",
        )
    end

    # Changes at the level of the solver tolerance are ignored
    δ = value.(s⁺) .- value.(s⁻)
    scale = max(maximum(abs, f.values), 1.0)
    adjusted = abs.(δ) .> 1e-7 * scale
    if !any(adjusted)
        @debug "Data points are $(curv_string)"
        return f
    end
    @debug "Data are not $(curv_string), $(count(adjusted)) point(s) adjusted, " *
           "largest change = $(round(maximum(abs, δ); sigdigits = 3))"
    values = [f.values[l] + δ[l] for l ∈ 1:K]
    return FunctionEvaluations(f.points, values)
end

# Create a full order convex approximation for the given data points,
# i.e. using a separate plane for each data point.
function _full_order_pwa(f::FunctionEvaluations{D}, optimizer, metric) where {D}
    K = length(f)

    m = Model(optimizer)
    set_silent(m)
    @variable(m, a[1:K, 1:D])
    @variable(m, b[1:K])
    @variable(m, s[1:K, 1:K] ≥ 0)
    for (l, (x, z)) ∈ enumerate(f)
        for k ∈ 1:K
            @constraint(m, z - s[l, k] == dot(a[k, :], x) + b[k])
        end
        @constraint(m, z ≤ dot(a[l, :], x) + b[l])
    end
    obj = AffExpr()
    if metric == :l2 || metric == :rms
        obj = sum(s .^ 2)
    elseif metric == :max
        @variable(m, smax)
        obj = smax
        @constraint(m, s .≤ smax)
    elseif metric == :l1
        obj = sum(s)
    elseif metric == :none
        # No regularization added to m
    else
        error("Unrecognized/unsupported metric type $(metric)")
    end
    @objective(m, Min, obj)
    optimize!(m)

    if termination_status(m) ∈ (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
        error(
            "No full order convex approximation exists, since the data points " *
            "are not convex (for example due to noise). Use " *
            "`enforce_curvature(data, Convex(), optimizer)` to adjust the " *
            "function values to be convex first.",
        )
    elseif !is_solved_and_feasible(m; allow_local = true)
        error(
            "Optimization failed: $(termination_status(m)) " *
            "($(raw_status(m)))",
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
# a metric that is minimized in the optimization.
function _reduced_order_pwa(
    f::FunctionEvaluations{D},
    p::Int,
    U,
    optimizer,
    μ = 1e-4,
) where {D}
    yield() # Avoid blocking of main thread causing hangs in some cases, see
    # https://github.com/sintefore/PiecewiseAffineApprox.jl/issues/11
    K = length(f)
    m = Model(optimizer)
    set_silent(m)
    @variable(m, ared[1:p, 1:D])
    @variable(m, bred[1:p])
    @variable(m, s[1:K] ≥ 0)
    @variable(m, sʳ[1:K, 1:p] ≥ 0)
    for (l, (x, z)) ∈ enumerate(f)
        for k ∈ 1:p
            @constraint(m, z - sʳ[l, k] == dot(ared[k, :], x) + bred[k])
        end
    end
    for (l, k) ∈ U
        z = f.values[l]
        x = f.points[l]
        @constraint(m, z ≤ dot(ared[k, :], x) + bred[k] + s[l])
    end
    @objective(m, Min, sum(s) + μ / p * sum(sʳ))
    optimize!(m)
    if !is_solved_and_feasible(m; allow_local = true)
        error(
            "Optimization failed: $(termination_status(m)) " *
            "($(raw_status(m)))",
        )
    end
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

# Approximation error per data point, compared with the tolerance: the mean
# absolute deviation (:l1), the root mean square deviation (:l2, :rms) or the
# maximum absolute deviation (:max)
function _tolerance_error(f::FunctionEvaluations, pwa::PWAFunc, metric)
    r = [abs(evaluate(pwa, x) - z) for (x, z) ∈ f]
    if metric == :l1
        return mean(r)
    elseif metric == :l2 || metric == :rms
        return sqrt(mean(abs2, r))
    elseif metric == :max
        return maximum(r)
    end
    throw(
        ArgumentError(
            "Unsupported metric :$metric for Progressive, use :l1, :l2 or :max",
        ),
    )
end

# Main algorithm
function _progressive_pwa(
    f::FunctionEvaluations,
    optimizer,
    δᵗᵒˡ,
    metric,
    full_order_metric = :l1,
)
    if metric ∉ (:l1, :l2, :rms, :max)
        throw(
            ArgumentError(
                "Unsupported metric :$metric for Progressive, " *
                "use :l1, :l2 or :max",
            ),
        )
    end

    # Find the full order convex approximation, will error if the points are
    # not convex
    pwa = _full_order_pwa(f, optimizer, full_order_metric)

    # Start with an approximation with one plane and allocate all points to that segment
    p = 1
    U = [(l, 1) for l ∈ 1:length(f)]

    # Increase the number of planes until the required tolerance is met
    pwa_red = _allocation_improvement(f, p, U, optimizer)
    used = []
    while _tolerance_error(f, pwa_red, metric) > δᵗᵒˡ && p < length(f.points)
        U, used = _increase_order(f, pwa_red, pwa, used)
        p += 1
        pwa_red = _allocation_improvement(f, p, U, optimizer)
    end
    @debug "Fitting finished, error = " *
           "$(round(_tolerance_error(f, pwa_red, metric); sigdigits = 3)), p = $p"
    return pwa_red
end

# Map algorithm structure to correct parameters
function _progressive_pwa(f::FunctionEvaluations, options::Progressive)
    return _progressive_pwa(
        f,
        options.optimizer,
        options.tolerance,
        options.metric,
    )
end
