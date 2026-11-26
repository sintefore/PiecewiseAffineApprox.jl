"""
    vertices_by_subsets(pwa::PWAFunc{C,D}[, bounds]; atol=1e-9) where {C,D}

Enumerate the vertices of a piecewise affine function.

For a convex PWA function in D dimensions, this method finds the extreme points of the
epigraph epi(f) = {(x,t) | t ≥ αᵢᵀx + βᵢ for all i}. Without bounds, the algorithm:
1. Considers all combinations of (D+1) hyperplanes
2. Solves the system of linear equations to find intersection points
3. Keeps the points that lie on the envelope of the function, i.e. t = f(x)

The number of combinations grows as K^(D+1) with the number of planes K, so this is
only suited for functions with few planes.

If `bounds = [(lo₁, hi₁), …, (lo_D, hi_D)]` is given, the epigraph is restricted to the box
lo ≤ x ≤ hi. The result then also contains the corners of the box and the points where
the planes cross the edges/faces of the box, so that the epigraph over the box is the
convex hull of the vertices plus the ray in the +t direction. The vertices are in this
case found by adding the planes one at a time, cutting the polyhedron found so far (the
double description method), with an effort that depends on the number of vertices
rather than on the number of combinations of planes.

For a concave PWA function, the method negates the function to convert it to convex,
finds vertices, then negates the function values back.

# Arguments
- `pwa::PWAFunc`: A piecewise affine function
- `bounds`: Optional box `[(lo, hi), ...]` with one tuple for each dimension
- `atol::Float64=1e-9`: Absolute tolerance for the envelope and box checks, and for merging
  duplicate vertices

# Returns
- `Vector{Tuple{Vector{Float64},Float64}}`: A vector of tuples `(x, z)` where `x` is the
  D-dimensional vertex and `z` is the function value at that point. Duplicates are removed.
"""
function vertices_by_subsets(pwa::PWAFunc{Convex,D}; atol = 1e-9) where {D}
    return _enumerate_vertices(pwa, nothing, nothing, atol)
end

function vertices_by_subsets(
    pwa::PWAFunc{Convex,D},
    bounds;
    atol = 1e-9,
) where {D}
    lo, hi = _check_bounds(bounds, D)
    return _enumerate_vertices_box(pwa, lo, hi, atol)
end

function vertices_by_subsets(
    pwa::PWAFunc{Concave,D},
    args...;
    atol = 1e-9,
) where {D}
    verts =
        vertices_by_subsets(PWAFunc{Convex,D}(pwa.planes), args...; atol = atol)
    return [(v[1], -v[2]) for v ∈ verts]
end

function _check_bounds(bounds, D)
    length(bounds) == D || throw(
        ArgumentError(
            "Expected bounds for $D dimension(s), got $(length(bounds))",
        ),
    )
    lo = Float64[first(b) for b ∈ bounds]
    hi = Float64[last(b) for b ∈ bounds]
    all(isfinite, lo) && all(isfinite, hi) ||
        throw(ArgumentError("Bounds must be finite, got $bounds"))
    all(lo .≤ hi) || throw(
        ArgumentError("Lower bounds must not exceed upper bounds, got $bounds"),
    )
    return lo, hi
end

# Vertices of the epigraph of a convex pwa function, optionally restricted to the box [lo, hi].
# A vertex is defined by D+1 active constraints. If k of these are box faces, they fix k
# coordinates of x, and the remaining D+1-k planes determine the free coordinates and t.
function _enumerate_vertices(pwa::PWAFunc{Convex,D}, lo, hi, atol) where {D}
    boxed = !isnothing(lo)
    verts = Tuple{Vector{Float64},Float64}[]
    for k ∈ 0:(boxed ? D : 0), fixed ∈ combinations(1:D, k)
        free = setdiff(1:D, fixed)
        n = length(free)
        for sides ∈ Iterators.product(fill((lo, hi), k)...)
            xfix = [sides[j][fixed[j]] for j ∈ 1:k]
            for S ∈ combinations(pwa.planes, n + 1)
                # Solve α_free⋅x_free - t = -β - α_fixed⋅x_fixed for the planes in S
                M = zeros(n + 1, n + 1)
                rhs = zeros(n + 1)
                for (row, plane) ∈ enumerate(S)
                    M[row, 1:n] = [plane.α[i] for i ∈ free]
                    M[row, n+1] = -1.0
                    rhs[row] =
                        -plane.β -
                        sum(plane.α[fixed[j]] * xfix[j] for j ∈ 1:k; init = 0.0)
                end
                rank(M) == n + 1 || continue
                sol = M \ rhs
                x = zeros(D)
                x[fixed] = xfix
                x[free] = sol[1:n]
                if boxed
                    all(lo .- atol .≤ x .≤ hi .+ atol) || continue
                    x = clamp.(x, lo, hi)
                end
                # The point is a vertex if it lies on the envelope of the function
                fx = evaluate(pwa, x)
                abs(sol[n+1] - fx) ≤ atol * max(1.0, abs(fx)) || continue
                if !any(v -> isapprox(v[1], x; atol = atol), verts)
                    push!(verts, (x, fx))
                end
            end
        end
    end
    return verts
end

# Vertices of the epigraph of a convex pwa function restricted to the box [lo, hi], by the
# double description method. The epigraph is first bounded by t ≤ t̄, with t̄ above the
# function over the box, and starts as the prism box × [t̲, t̄], with t̲ below all planes.
# The half-spaces t ≥ αᵢ⋅x + βᵢ of the planes are then added one at a time: the vertices
# below the plane are removed, and a new vertex is added on each edge between a vertex
# above and a vertex below the plane. Two vertices are joined by an edge if and only if no
# other vertex is active in all of the constraints that are active in both, which also
# holds for degenerate vertices. The vertices of the epigraph are finally those that are
# not on t = t̄.
#
# The constraints are numbered with the planes first, then the lower and upper face of the
# box for each coordinate, and finally t ≥ t̲ and t ≤ t̄.
function _enumerate_vertices_box(pwa::PWAFunc{Convex,D}, lo, hi, atol) where {D}
    planes = pwa.planes
    K = length(planes)
    verts = Tuple{Vector{Float64},Float64}[]
    K == 0 && return verts
    A = [p.α[i] for i ∈ 1:D, p ∈ planes]
    β = [p.β for p ∈ planes]
    # The range of each plane over the box
    pmin = [β[k] + sum(min(A[i, k] * lo[i], A[i, k] * hi[i]) for i ∈ 1:D) for k ∈ 1:K]
    pmax = [β[k] + sum(max(A[i, k] * lo[i], A[i, k] * hi[i]) for i ∈ 1:D) for k ∈ 1:K]
    scale = max(1.0, maximum(abs, pmin), maximum(abs, pmax))
    τ = atol * scale
    t̲, t̄ = minimum(pmin) - scale, maximum(pmax) + scale
    bottom, top = K + 2D + 1, K + 2D + 2

    # The corners of the prism. A coordinate with lo = hi is on both of its faces.
    P = Vector{Float64}[]     # points (x, t)
    active = BitSet[]         # active constraints of each point
    for sides ∈ Iterators.product((lo[i] < hi[i] ? (1, 2) : (0,) for i ∈ 1:D)...)
        x = [sides[i] == 2 ? hi[i] : lo[i] for i ∈ 1:D]
        faces = BitSet()
        for i ∈ 1:D
            sides[i] != 2 && push!(faces, K + 2i - 1)
            sides[i] != 1 && push!(faces, K + 2i)
        end
        push!(P, [x; t̲])
        push!(active, union(faces, bottom))
        push!(P, [x; t̄])
        push!(active, union(faces, top))
    end

    slack = Float64[]
    for k ∈ 1:K
        resize!(slack, length(P))
        for (v, p) ∈ enumerate(P)
            s = p[D+1] - β[k]
            for i ∈ 1:D
                s -= A[i, k] * p[i]
            end
            slack[v] = s
        end
        out = findall(<(-τ), slack)
        for v ∈ eachindex(P)
            abs(slack[v]) ≤ τ && push!(active[v], k)
        end
        isempty(out) && continue
        inside = findall(>(τ), slack)
        newP, newactive = Vector{Float64}[], BitSet[]
        for w ∈ out, u ∈ inside
            common = intersect(active[u], active[w])
            # An edge has at least D active constraints
            length(common) ≥ D || continue
            adjacent = true
            for v ∈ eachindex(P)
                if v != u && v != w && common ⊆ active[v]
                    adjacent = false
                    break
                end
            end
            adjacent || continue
            λ = slack[u] / (slack[u] - slack[w])
            push!(newP, P[u] .+ λ .* (P[w] .- P[u]))
            push!(newactive, push!(common, k))
        end
        deleteat!(P, out)
        deleteat!(active, out)
        append!(P, newP)
        append!(active, newactive)
    end

    for (p, a) ∈ zip(P, active)
        top ∈ a && continue
        x = clamp.(p[1:D], lo, hi)
        if !any(v -> isapprox(v[1], x; atol = atol), verts)
            push!(verts, (x, evaluate(pwa, x)))
        end
    end
    return verts
end
