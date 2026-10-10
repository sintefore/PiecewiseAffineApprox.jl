# A general 1D piecewise linear function defined by points 
# in sorted order and their associated function value 
struct PWAFunction1D
    x::Vector{Float64}
    z::Vector{Float64}

    function PWAFunction1D(x::Vector{Float64}, z::Vector{Float64})
        @assert length(x) == length(z)
        @assert issorted(x)
        return new(x, z)
    end
end

function PWAFunction1D(x::Vector, z::Vector)
    return PWAFunction1D(
        convert(Vector{Float64}, x),
        convert(Vector{Float64}, z),
    )
end

function _asconvex(pwa::PWAFunction1D)
    N = length(pwa.x)
    x = pwa.x
    z = pwa.z

    c = [(z[i+1] - z[i]) / (x[i+1] - x[i]) for i ∈ 1:(N-1)]
    d = [z[i] - c[i] * x[i] for i ∈ 1:(N-1)]

    return PWAFunc{Convex,1}(Plane.(c, d))
end

# Create a convex approximation to a general pwa function
# by perturbing function values as little as possible (l1-deviation)
function _convexify1D(pwa::PWAFunction1D, options)
    N = length(pwa.x)
    𝒩 = 1:N
    x = pwa.x
    z = pwa.z

    m = Model()

    δ⁺ = @variable(m, [𝒩], lower_bound = 0)
    δ⁻ = @variable(m, [𝒩], lower_bound = 0)

    @objective(m, Min, sum(δ⁺[i] + δ⁻[i] for i ∈ 𝒩))
    for i ∈ 2:(N-1)
        @constraint(
            m,
            (z[i] + δ⁺[i] - δ⁻[i] - z[i-1] - δ⁺[i-1] + δ⁻[i-1]) /
            (x[i] - x[i-1]) <=
            (z[i+1] + δ⁺[i+1] - δ⁻[i+1] - z[i] - δ⁺[i] + δ⁻[i]) /
            (x[i+1] - x[i])
        )
    end

    set_optimizer(m, options.optimizer)
    optimize!(m)
    if !is_solved_and_feasible(m)
        error(
            "Optimization failed: $(termination_status(m)) " *
            "($(raw_status(m)))",
        )
    end

    dp = value.(δ⁺)
    dn = value.(δ⁻)

    return PWAFunction1D(x, [z[i] + dp[i] - dn[i] for i ∈ 𝒩])
end

_convexify1D(x, z, options) = _convexify1D(PWAFunction1D(x, z), options)

# Create a pwa interpolant to the given points (sorted by x) using exactly
# options.planes segments with end points at data points. The segments are
# chosen to minimize the total deviation in options.metric: the sum over the
# segments for :l1 and :l2, and the maximum over the segments for :max.
function _interpolatepw(x, z, options)
    @assert length(x) == length(z)
    p = _segment_penalties(x, z, options.metric)
    minimax = options.metric == :max
    knots, _ = _shortest_path_segments(p, options.planes; minimax)
    return PWAFunction1D(x[knots], z[knots])
end

# Penalty p[i, j] (i < j) for approximating the points i, ..., j by the line
# segment between points i and j: the sum of absolute deviations (:l1), the sum
# of squared deviations (:l2) or the maximum absolute deviation (:max).
function _segment_penalties(x, z, metric)
    if metric == :l2
        return _segment_penalties_l2(x, z)
    elseif metric ∉ (:l1, :max)
        error("Unrecognized/unsupported metric $(metric)")
    end
    N = length(x)
    p = zeros(N, N)
    for i ∈ 1:N-1, j ∈ i+1:N
        c = (z[j] - z[i]) / (x[j] - x[i])
        acc = 0.0
        for k ∈ i+1:j-1
            r = abs(z[i] + c * (x[k] - x[i]) - z[k])
            acc = metric == :l1 ? acc + r : max(acc, r)
        end
        p[i, j] = acc
    end
    return p
end

# Sum of squared deviations for all segments in O(N²) using prefix sums.
# The data is centered first to reduce cancellation errors.
function _segment_penalties_l2(x, z)
    N = length(x)
    x = x .- mean(x)
    z = z .- mean(z)
    Sx = [0.0; cumsum(x)]
    Sz = [0.0; cumsum(z)]
    Sxx = [0.0; cumsum(x .^ 2)]
    Szz = [0.0; cumsum(z .^ 2)]
    Sxz = [0.0; cumsum(x .* z)]
    p = zeros(N, N)
    for i ∈ 1:N-1, j ∈ i+1:N
        # Deviations at points k = i, ..., j are α + c xₖ - zₖ
        c = (z[j] - z[i]) / (x[j] - x[i])
        α = z[i] - c * x[i]
        n = j - i + 1
        sx, sz = Sx[j+1] - Sx[i], Sz[j+1] - Sz[i]
        sxx, szz, sxz = Sxx[j+1] - Sxx[i], Szz[j+1] - Szz[i], Sxz[j+1] - Sxz[i]
        r² =
            n * α^2 + c^2 * sxx + szz + 2α * c * sx - 2α * sz - 2c * sxz
        p[i, j] = max(r², 0.0)
    end
    return p
end

# Find the path 1 → N through increasing indices with exactly K edges that
# minimizes the sum of edge costs p[i, j], or with minimax = true the maximum
# edge cost. Returns the indices on the path and the minimal cost of reaching
# N using 1, ..., K edges.
#
# For minimax, ties are broken by the sum of edge costs: a second pass
# minimizes the sum using only edges with cost at most the minimax value.
function _shortest_path_segments(p, K; minimax = false)
    if !minimax
        return _dp_segments(p, K, +)
    end
    _, costs = _dp_segments(p, K, max)
    q = [v ≤ costs[K] ? v : Inf for v ∈ p]
    path, _ = _dp_segments(q, K, +)
    return path, costs
end

# Dynamic program in O(KN²) for the path 1 → N with exactly K edges that
# minimizes the edge costs p[i, j] aggregated by op (+ or max). Ties are broken
# towards smaller indices.
function _dp_segments(p, K, op)
    N = size(p, 1)
    @assert 1 ≤ K ≤ N - 1
    # cost[s+1, j]: minimal cost of reaching j using s edges
    cost = fill(Inf, K + 1, N)
    pred = zeros(Int, K + 1, N)
    cost[1, 1] = 0.0
    for s ∈ 1:K, j ∈ s+1:N
        for i ∈ s:j-1
            v = op(cost[s, i], p[i, j])
            if v < cost[s+1, j]
                cost[s+1, j] = v
                pred[s+1, j] = i
            end
        end
    end
    path = [N]
    for s ∈ K:-1:1
        pushfirst!(path, pred[s+1, path[1]])
    end
    return path, cost[2:end, N]
end

function _convex_linearization_ipol(x, z, options)
    @assert length(x) == length(z)
    x, z = _sort_unique(x, z)
    N = length(x)
    if !(1 ≤ options.planes ≤ N - 1)
        throw(
            ArgumentError(
                "Interpol with planes = $(options.planes) requires at least " *
                "planes + 1 distinct points, got $N",
            ),
        )
    end
    return _asconvex(_convexify1D(_interpolatepw(x, z, options), options))
end

# Sort the points by x and merge points with identical x, using the average
# of their function values
function _sort_unique(x, z)
    xs = Float64[]
    zs = Float64[]
    n = Int[]
    for i ∈ sortperm(x)
        if !isempty(xs) && x[i] == xs[end]
            n[end] += 1
            zs[end] += (z[i] - zs[end]) / n[end]
        else
            push!(xs, x[i])
            push!(zs, z[i])
            push!(n, 1)
        end
    end
    return xs, zs
end
