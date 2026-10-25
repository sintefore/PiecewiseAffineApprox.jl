"""
    vertices_by_subsets(pwa::PWAFunc{C,D}[, bounds]; atol=1e-9) where {C,D}

Enumerate vertices of a piecewise affine function by examining all subsets of (D+1)
active constraints.

For a convex PWA function in D dimensions, this method finds the extreme points of the
epigraph epi(f) = {(x,t) | t ≥ αᵢᵀx + βᵢ for all i}. The algorithm:
1. Considers all combinations of (D+1) hyperplanes
2. Solves the system of linear equations to find intersection points
3. Keeps the points that lie on the envelope of the function, i.e. t = f(x)

If `bounds = [(lo₁, hi₁), …, (lo_D, hi_D)]` is given, the epigraph is restricted to the box
lo ≤ x ≤ hi, and the box faces are included as candidate active constraints. The result
then also contains the corners of the box and the points where the planes cross the
edges/faces of the box, so that the epigraph over the box is the convex hull of the
vertices plus the ray in the +t direction.

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
    return _enumerate_vertices(pwa, lo, hi, atol)
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
