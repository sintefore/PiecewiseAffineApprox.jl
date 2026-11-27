
function _convex_linearization_mb_single(
    X::Matrix,
    z::Vector,
    K,
    lᵐᵃˣ,
    metric,
    optimizer,
    strict,
    rng;
    polish = 5,
    quadratic = false,
)
    𝒫 = _random_partition(rng, X, K)
    pwa = _refine_two_phase(X, z, 𝒫, lᵐᵃˣ, metric, optimizer, strict, polish, false)
    if quadratic
        # Quadratic pieces are refined from the random partition, and from the
        # partition of the fitted planes, and the best of these and the planes
        # (quadratic pieces without curvature) is kept, so that the result is
        # never worse than planes from the same start. Pieces with a large
        # curvature fitted to a few points tend to take over the points of the
        # others, and which start avoids this depends on the data.
        candidates = [
            _as_quadratic(pwa),
            _refine_two_phase(X, z, 𝒫, lᵐᵃˣ, metric, optimizer, strict, polish, true),
            _refine_two_phase(X, z, _update_partition(X, pwa), lᵐᵃˣ, metric, optimizer, strict, polish, true),
        ]
        pwa = argmin(f -> approximation_error(f, X, z, metric), candidates)
    end
    e = approximation_error(pwa, X, z, metric)
    return e => pwa
end

# Refine the partition with the requested fit. For the fits that solve an
# optimization problem, in two phases: the partition is first refined with
# closed form least squares, which is fast and moves the pieces smoothly, and
# then for a few iterations with the requested fit. This is both faster and
# better than the requested fit alone, whose solutions at vertices of the
# feasible region tend to get stuck in limit cycles early.
function _refine_two_phase(X, z, 𝒫, lᵐᵃˣ, metric, optimizer, strict, polish, quadratic)
    if strict == :none && metric ∈ (:l2, :rms)
        return _refine_partition(X, z, 𝒫, lᵐᵃˣ, metric, optimizer, strict; quadratic)
    end
    pwa₁ = _refine_partition(X, z, 𝒫, lᵐᵃˣ, :l2, nothing, :none; quadratic)
    𝒫 = _update_partition(X, pwa₁)
    return _refine_partition(X, z, 𝒫, min(lᵐᵃˣ, polish), metric, optimizer, strict; quadratic)
end

# A convex piecewise affine function as a piecewise quadratic one
function _as_quadratic(pwa::PWAFunc{Convex,D}) where {D}
    return PWQFunc{Convex,D}([Quadratic(zeros(D, D), p.α, p.β) for p ∈ pwa.planes])
end

# Finds a pwa convex approximation for the provided data
# X is a matrix with a column for each data point and z is a vector with the 
# corresponding function values
function _convex_linearization_mb(X::Matrix, z::Vector, options::Cluster)
    @assert(size(X, 2) == length(z))

    Nᵗʳ = options.trials    # Number of trials
    lᵐᵃˣ = options.itlim    # Iteration limit
    K = options.planes
    metric = options.metric
    strict = options.strict
    optimizer = options.optimizer

    K ≥ 1 || throw(ArgumentError("Cluster requires at least one plane"))
    if isnothing(optimizer) && !(strict == :none && metric ∈ (:l2, :rms))
        error(
            "Cluster with metric = :$metric and strict = :$strict requires an optimizer. " *
            "Only metric = :l2 with strict = :none can be used without one.",
        )
    end

    Nᵗʳ ≥ 1 || throw(ArgumentError("Cluster requires at least one trial"))
    _check_option(options, "pieces", options.pieces, (:affine, :quadratic))

    # A single trial; several trials are run through Restart, see _fit
    e, pwa = _convex_linearization_mb_single(
        X,
        z,
        K,
        lᵐᵃˣ,
        metric,
        optimizer,
        strict,
        options.rng;
        quadratic = options.pieces == :quadratic,
    )
    @debug "Cluster trial finished with error $e ($metric)"
    return pwa
end

# Euclidean distance between two points
_dist(x, y) = sqrt(sum((x - y) .^ 2))

# Draw K samples from a normal distribution with the same mean and covariance
# as the columns of X. The covariance is factored by its symmetric square root
# rather than Cholesky so that singular covariances (e.g. collinear points)
# work. The symmetric square root does not depend on the choice of basis in
# eigenspaces with repeated eigenvalues (e.g. the identity covariance of
# whitened data), so the samples are robust to rounding in the data.
function _sample_normal(rng::AbstractRNG, X, K)
    μ = vec(mean(X, dims = 2))
    size(X, 2) < 2 && return repeat(μ, 1, K)
    λ, V = eigen(Symmetric(cov(X, dims = 2)))
    L = V * Diagonal(sqrt.(max.(λ, 0))) * V'
    return μ .+ L * randn(rng, length(μ), K)
end

# Create a random partition of the points into K sets
# by generating K random points and group each data point to the
# closest of these random points.
function _random_partition(rng::AbstractRNG, X, K)
    P = _sample_normal(rng, X, K)

    𝒫 = [Int[] for j ∈ 1:K]
    for (i, x) ∈ enumerate(eachcol(X))
        # Find the nearest point amongst the p's
        jmin = argmin(_dist(x, P[:, j]) for j ∈ 1:K)
        push!(𝒫[jmin], i)
    end
    return 𝒫
end

# Improve a partition by fitting a hyperplane for all points 
# belonging to each subset of the partition and then 
# updating the partition such that each point is associated with the
# plane being active at the point. 
#
# Subsets that do not spread out in all directions spanned by the data (e.g. fewer
# than D+1 points, or points lying in a lower dimensional subspace) do not determine
# a plane. Fitting one anyway can give a plane that dominates all others (e.g. a
# flat plane through a single high point). Such subsets are dropped and their
# points reassigned to the active plane among the remaining ones. To keep the number of planes,
# the partition is then refilled by splitting the subsets with the largest error
# (deviates from Magnani & Boyd, who drop empty subsets and continue with
# fewer planes).
#
# The process terminates when the partition no longer changes, or after a given
# number of iterations, returning the best pwa function found (as measured by
# metric). Only fits to complete partitions are considered, to retain the
# properties given by strict.
function _refine_partition(
    X,
    z,
    𝒫,
    lᵐᵃˣ,
    metric,
    optimizer,
    strict;
    refill = true,
    minvar = 1e-4,
    quadratic = false,
)
    D, N = size(X)
    K = length(𝒫)
    determined = _determined_test(X, minvar; quadratic)
    𝒫 = filter(determined, 𝒫)
    if isempty(𝒫)
        𝒫 = [collect(1:N)]
    elseif sum(length, 𝒫) < N
        # Reassign the points of the dropped subsets to the active plane among
        # least squares fits to the remaining ones, so that the partition is
        # complete. A fit to an incomplete partition would not hold strict at
        # the missing points.
        ls = _fit_planes(X, z, 𝒫, :l2, nothing, :none; quadratic)
        𝒫 = filter!(!isempty, _update_partition(X, ls))
    end
    refill && _refill!(𝒫, X, ones(N), K, determined)
    fit(q) = _fit_planes(X, z, q, metric, optimizer, strict; quadratic)

    pwaᵇᵉˢᵗ = fit(𝒫)
    eᵇᵉˢᵗ = approximation_error(pwaᵇᵉˢᵗ, X, z, metric)
    pwa = pwaᵇᵉˢᵗ
    for it ∈ 1:lᵐᵃˣ
        if it > 1
            pwa = fit(𝒫)
            e = approximation_error(pwa, X, z, metric)
            if e < eᵇᵉˢᵗ
                pwaᵇᵉˢᵗ, eᵇᵉˢᵗ = pwa, e
            end
        end

        # Update the partition with the active planes, reassigning points
        # in undetermined subsets to the remaining planes
        𝒫ⁿᵉʷ = _update_partition(X, pwa)
        active = pwa
        if !all(determined, 𝒫ⁿᵉʷ)
            keep = determined.(𝒫ⁿᵉʷ)
            if any(keep)
                active = _keep_pieces(pwa, keep)
                𝒫ⁿᵉʷ = _update_partition(X, active)
            end
            filter!(p -> !isempty(p), 𝒫ⁿᵉʷ)
        end
        if refill && length(𝒫ⁿᵉʷ) < K
            r = [evaluate(active, x) - z[i] for (i, x) ∈ enumerate(eachcol(X))]
            w = metric ∈ (:l2, :rms) ? r .^ 2 : abs.(r)
            _refill!(𝒫ⁿᵉʷ, X, w, K, determined)
        end
        if 𝒫ⁿᵉʷ == 𝒫
            # Converged, pwa is the fit to this partition
            return pwaᵇᵉˢᵗ
        end
        𝒫 = 𝒫ⁿᵉʷ
    end

    # Iteration limit reached, also consider the fit to the final partition
    pwa = fit(𝒫)
    return approximation_error(pwa, X, z, metric) < eᵇᵉˢᵗ ? pwa : pwaᵇᵉˢᵗ
end

# Fit a plane, or a convex quadratic with diagonal curvature, to each subset of
# the partition 𝒫. For strict = :outer, all pieces are required to lie below
# all data points, so that the function is an outer approximation also when
# the partition has not converged. A quadratic piece is also required to lie
# below the data outside its subset (except for strict = :inner, where this
# may be infeasible), see _least_squares_fit_quadratic.
function _fit_planes(X, z, 𝒫, metric, optimizer, strict; quadratic = false)
    D = size(X, 1)
    if quadratic
        pieces = map(𝒫) do p
            # The points outside the subset, which the piece must not exceed
            out = setdiff(axes(X, 2), p)
            γ, a, b = _local_fit(
                X[:, p], z[p], metric, optimizer, strict, X, z;
                quadratic, outside = (X[:, out], z[out]),
            )
            return Quadratic(Matrix(Diagonal(γ)), a, b)
        end
        return PWQFunc{Convex,D}(pieces)
    end
    pwa = PWAFunc{Convex,D}()
    for p ∈ 𝒫
        a, b = _local_fit(X[:, p], z[p], metric, optimizer, strict, X, z)
        _addplane!(pwa, a, b)
    end
    return pwa
end

# The function with the pieces kept
_keep_pieces(pwa::PWAFunc{C,D}, keep) where {C,D} = PWAFunc{C,D}(pwa.planes[keep])
_keep_pieces(f::PWQFunc{C,D}, keep) where {C,D} = PWQFunc{C,D}(f.pieces[keep])

# Returns a function testing whether a subset of the points (columns of X)
# determines a plane: its sample variance must be at least minvar in all directions
# where the full data set has variation. X is assumed to be normalized, i.e. to
# have unit variance in these directions. A quadratic piece with diagonal
# curvature also needs at least 2r + 1 points, with r the number of directions.
function _determined_test(X, minvar; quadratic = false)
    r = if size(X, 2) < 2
        0
    else
        count(>(0.5), eigvals(Symmetric(cov(X, dims = 2))))
    end
    return function (p)
        length(p) > (quadratic ? 2r : r) || return false
        r == 0 && return true
        λ = eigvals(Symmetric(cov(X[:, p], dims = 2, corrected = false)))
        return count(≥(minvar), λ) ≥ r
    end
end

# Increase the number of subsets in the partition 𝒫 to K by repeatedly splitting
# the subset with the largest total weight w (e.g. approximation error), such that
# all subsets still determine a plane.
function _refill!(𝒫, X, w, K, determined)
    splittable = trues(length(𝒫))
    while length(𝒫) < K
        candidates = [j for j ∈ eachindex(𝒫) if splittable[j]]
        isempty(candidates) && break
        j = argmax(j -> sum(w[𝒫[j]]), candidates)
        p₁, p₂ = _split_principal(X, 𝒫[j])
        if !determined(p₁) || !determined(p₂)
            splittable[j] = false
        else
            𝒫[j] = p₁
            push!(𝒫, p₂)
            push!(splittable, true)
        end
    end
    return 𝒫
end

# Split the points p in two halves at the median along their principal direction
function _split_principal(X, p)
    Xc = X[:, p] .- mean(X[:, p], dims = 2)
    d = eigen(Symmetric(Xc * Xc')).vectors[:, end]
    s = vec(d' * Xc)
    m = median(s)
    return p[s .< m], p[s .>= m]
end

# Finds the hypeplane best fitting the given subset of points 
# using the given metric and possible restrictions on whether
# the plane should be strictly above the points (strict = :inner) or
# below the points Xᵒ with values zᵒ (strict = :outer)
function _local_fit(X̄, z̄, metric, optimizer, strict, Xᵒ = X̄, zᵒ = z̄; quadratic = false, outside = nothing)
    if strict == :none && (metric == :l2 || metric == :rms)
        return quadratic ? _least_squares_fit_quadratic(X̄, z̄, outside) : _least_squares_fit(X̄, z̄)
    end

    M, N = size(X̄)

    # Create an optimization model to find the best a and b such that  ax + b ≈ z
    m = Model()

    if Threads.nthreads() > 1
        # Use threading for individual problems, avoid oversubscribing by limiting threads for each problem
        MOI.set(m, MOI.NumberOfThreads(), 1)
    end

    @variable(m, a[1:M])
    @variable(m, b)
    @variable(m, ẑ[1:N])
    # The diagonal curvature of a quadratic piece, nonnegative for convexity
    γ = quadratic ? @variable(m, [1:M], lower_bound = 0) : zeros(M)
    piece(x) = sum(γ[j] * x[j]^2 + a[j] * x[j] for j ∈ 1:M) + b

    for i ∈ 1:N
        @constraint(m, ẑ[i] == piece(view(X̄, :, i)))
    end

    if strict == :inner
        @constraint(m, ẑ .≥ z̄)
    elseif strict == :outer
        for i ∈ axes(Xᵒ, 2)
            @constraint(m, piece(view(Xᵒ, :, i)) ≤ zᵒ[i])
        end
    end
    if quadratic && strict == :none && !isnothing(outside)
        Xᵒᵘᵗ, zᵒᵘᵗ = outside
        for i ∈ axes(Xᵒᵘᵗ, 2)
            @constraint(m, piece(view(Xᵒᵘᵗ, :, i)) ≤ zᵒᵘᵗ[i])
        end
    end

    obj = AffExpr()
    if metric == :l2 || metric == :rms
        obj = sum((z̄[i] - ẑ[i])^2 for i ∈ 1:N)
    elseif metric == :max
        @variable(m, t)
        obj = t
        for i ∈ 1:N
            @constraint(m, t ≥ (z̄[i] - ẑ[i]))
            @constraint(m, t ≥ (ẑ[i] - z̄[i]))
        end
    elseif metric == :l1
        @variable(m, t[1:N])
        obj = sum(t[i] for i ∈ 1:N)
        for i ∈ 1:N
            @constraint(m, t[i] ≥ (z̄[i] - ẑ[i]))
            @constraint(m, t[i] ≥ (ẑ[i] - z̄[i]))
        end
    else
        error("Unrecognized/unsupported metric type $(metric)")
    end

    # TODO: consider adding regularization term

    @objective(m, Min, obj)

    set_optimizer(m, optimizer)
    set_silent(m)
    optimize!(m)

    if !is_solved_and_feasible(m; allow_local = true)
        error(
            "Optimization failed: $(termination_status(m)) ($(raw_status(m)))",
        )
    end

    ā = JuMP.value.(a)
    b̄ = JuMP.value.(b)
    quadratic && return max.(JuMP.value.(γ), 0.0), ā, b̄
    return ā, b̄
end

# Closed form least squares fit of the plane aᵀx + b to the points (Magnani & Boyd, eq. 9).
# Solved via the centered normal equations, with a small ridge term on a to handle
# subsets that do not span all directions (e.g. fewer than D+1 points). This gives
# zero slope along directions with no variation in the data.
function _least_squares_fit(X̄, z̄)
    x̄ = vec(mean(X̄, dims = 2))
    Xc = X̄ .- x̄
    zc = z̄ .- mean(z̄)
    C = Xc * Xc'
    λ = 1e-8 * tr(C) / size(C, 1)
    a = λ > 0 ? (C + λ * I) \ (Xc * zc) : zeros(size(C, 1))
    b = mean(z̄) - dot(a, x̄)
    return a, b
end

# Least squares fit of a convex quadratic with diagonal curvature to the points
# X̄ with values z̄, z ≈ Σⱼ γⱼxⱼ² + aᵀx + b with γ ≥ 0, which lies below the
# values zᵒ at the points Xᵒ given by outside = (Xᵒ, zᵒ) (or nothing).
#
# The points outside are those of the other pieces. A plane fitted to a convex
# function lies below it away from its points, but a quadratic piece can curve
# up faster than the function, and the maximum of the pieces then exceeds the
# data where other pieces are fitted. With the constraint, the deviation of the
# maximum at a point is at most that of the piece fitted to it. Without enough
# points for the curvature, it is zero.
function _least_squares_fit_quadratic(X̄, z̄, outside = nothing)
    D, N = size(X̄)
    quadratic = N > 2D
    # The values of the coefficients at the points (columns of Y), as rows
    function features(Y)
        Yᵗ = permutedims(Y)
        return quadratic ? [Yᵗ .^ 2 Yᵗ ones(size(Y, 2))] : [Yᵗ ones(size(Y, 2))]
    end
    E = features(X̄)
    n = size(E, 2)
    # Constraints Gθ ≥ h: nonnegative curvature, and below the outside values
    G = quadratic ? [Matrix(1.0I, D, D) zeros(D, D + 1)] : zeros(0, n)
    h = zeros(size(G, 1))
    θ = _lsi(E, z̄, G, h)
    if !isnothing(outside) && !isempty(outside[2])
        # The constraints of the outside points are added as they are violated,
        # the most violated first and a few at a time, as only a few are binding
        # and the solution of the least distance problem is accurate for few
        # constraints
        Xᵒ, zᵒ = outside
        Gᵒ = -features(Xᵒ)
        tol = 1e-9 * max(maximum(abs, z̄), maximum(abs, zᵒ), 1.0)
        added = falses(length(zᵒ))
        for _ ∈ 1:50
            slack = Gᵒ * θ .+ zᵒ
            violated = findall(.!added .& (slack .< -tol))
            isempty(violated) && break
            worst = violated[partialsortperm(slack[violated], 1:min(n, length(violated)))]
            added[worst] .= true
            θ = _lsi(E, z̄, [G; Gᵒ[added, :]], [h; -zᵒ[added]])
        end
        # Any violation that remains, e.g. from rounding, is removed by lowering
        # the piece
        θ[end] += min(minimum(Gᵒ * θ .+ zᵒ), 0.0)
    end
    quadratic || return zeros(D), θ[1:D], θ[D+1]
    return max.(θ[1:D], 0.0), θ[D+1:2D], θ[2D+1]
end

# Least squares with inequality constraints, min ‖Eθ - f‖ subject to Gθ ≥ h,
# by the reduction to a least distance problem that is solved by nonnegative
# least squares (Lawson & Hanson, 1974, ch. 23). With the QR factorization
# E = QR and θ = R⁻¹(y + Qᵀf), the problem is min ‖y‖ subject to
# GR⁻¹y ≥ h - GR⁻¹Qᵀf. A small ridge term makes R invertible also for points
# that do not determine θ. The constraints must be feasible.
function _lsi(E, f, G, h)
    n = size(E, 2)
    λ = 1e-10 * max(sum(abs2, E) / n, 1.0)
    F = qr([E; sqrt(λ) * Matrix(1.0I, n, n)])
    R = F.R
    c = (F.Q' * [f; zeros(n)])[1:n]
    isempty(h) && return R \ c
    Ĝ = G / R
    ĥ = h .- Ĝ * c
    # The least distance problem min ‖y‖ subject to Ĝy ≥ ĥ, from the
    # nonnegative least squares solution u of [Ĝᵀ; ĥᵀ] u ≈ (0, …, 0, 1)
    A = [Ĝ'; ĥ']
    u = _nnls(A, [zeros(n); 1.0])
    r = A * u
    r[end] -= 1.0
    abs(r[end]) > eps() || throw(ArgumentError("Infeasible constraints in the fit of a piece"))
    y = -r[1:n] ./ r[end]
    return R \ (y .+ c)
end

# The piece being active at the point x for a convex function
function _active(pwa::PWAFunc{Convex,D}, x) where {D}
    return argmax(collect(evaluate(p, x) for p ∈ pwa.planes))
end
function _active(f::PWQFunc{Convex,D}, x) where {D}
    return argmax(collect(evaluate(q, x) for q ∈ f.pieces))
end

# Creating an updated partition by associating each data point to
# the hyperplane being active for the data point
function _update_partition(X, pwa)
    𝒫 = [Int[] for _ ∈ 1:_planes(pwa)]
    for (i, x) ∈ enumerate(eachcol(X))
        push!(𝒫[_active(pwa, x)], i)
    end
    return 𝒫
end

function _fit(a::Cluster, input::FunctionEvaluations{D}) where {D}
    if a.trials > 1
        # The trials are restarts of single trials with seeded generators
        single = with(a; trials = 1)
        r = Restart(single; trials = a.trials, rng = a.rng, metric = a.metric)
        return _fit(r, input)
    end
    # Work on whitened points to make the heuristic invariant to affine
    # transformations of the input (Magnani & Boyd, Sect. 1.3)
    normalized, n = _normalize(input; whiten = true, values = false)
    X̃ = _matrix(normalized)
    if a.pieces == :quadratic
        X̃, n = _principal_axes(X̃, n, _matrix(input))
    end
    pwa = _convex_linearization_mb(X̃, input.values, a)
    return _denormalize(pwa, n)
end

# The whitened points X̃ = T (x - μ) rotated to the principal axes of the
# points, with the normalization updated. The covariance stays the identity,
# but a diagonal curvature in these coordinates is along the principal axes,
# which move with the data, so that a quadratic fit is invariant to affine
# transformations of the points (up to ties of the axes). Whitening by itself
# keeps the coordinates as close to the original ones as possible.
function _principal_axes(X̃, n::Normalization, X)
    # T = W S, with S scaling by the ranges, and the axes are those of S (x - μ)
    Xs = (X .- n.μ) ./ (2 .* n.σ)
    D = size(Xs, 1)
    V = Matrix{Float64}(I, D, D)
    if size(Xs, 2) ≥ 2
        C = cov(Xs, dims = 2)
        # Uncorrelated coordinates, e.g. a grid, are principal axes. They
        # are kept, also where variances are equal and the axes not unique.
        correlated = any(abs(C[i, j]) > 1e-8 * sqrt(C[i, i] * C[j, j]) for i ∈ 1:D, j ∈ 1:D if i != j)
        correlated && (V = eigen(Symmetric(C)).vectors)
    end
    return V' * X̃, Normalization(n.μ, V' * n.T, n.ν, n.s, n.σ)
end
