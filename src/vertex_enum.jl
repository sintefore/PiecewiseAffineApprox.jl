"""
    vertices_by_subsets(pwa::PWAFunc{C,D}; atol=1e-9) where {C,D}

Enumerate vertices of a piecewise affine function by examining all subsets of (D+1) planes.

For a convex PWA function in D dimensions, this method finds the extreme points of the
epigraph epi(f) = {(x,t) | t ≥ αᵢᵀx + βᵢ for all i}. The algorithm:
1. Considers all combinations of (D+1) hyperplanes
2. Solves the system of linear equations to find intersection points
3. Checks feasibility: the point must satisfy t ≥ αᵢᵀx + βᵢ for all planes
4. Verifies the point lies on the envelope of the function

For a concave PWA function, the method negates the function to convert it to convex,
finds vertices, then negates the function values back.

# Arguments
- `pwa::PWAFunc`: A piecewise affine function
- `atol::Float64=1e-9`: Absolute tolerance for feasibility and envelope checks

# Returns
- `Vector{Tuple}`: A vector of tuples `(x, z)` where `x` is the D-dimensional vertex
  and `z` is the function value at that point
```
"""
function vertices_by_subsets(pwa::PWAFunc{Convex,D}; atol=1e-9) where {D}

    m = PiecewiseAffineApprox._planes(pwa)
    n = D
    verts = []
    for S in combinations(1:m, n+1)
        M = zeros(n+1, n+1)
        rhs = zeros(n+1)
        for (row, p) in enumerate(S)
            plane = pwa.planes[p]
            M[row, 1:n] = collect(plane.α)
            M[row, n+1] = -1.0
            rhs[row] = -plane.β
        end
        # Solve M * [x; t] = rhs
        if rank(M) == n+1
            sol = M \ rhs
            x = sol[1:n]; t = sol[n+1]
            # feasibility: t >= α x + β for all planes
            if all(t >= dot(p.α, x) + p.β - atol for p in pwa.planes)
                # on the envelope:
                if abs(t - evaluate(pwa, x)) <= atol
                    push!(verts, (x, evaluate(pwa, x)))
                end
            end
        end
    end
    return verts
end

function vertices_by_subsets(pwa::PWAFunc{Concave,D}; atol=1e-9) where {D}
    verts = vertices_by_subsets((PWAFunc{Convex,D}(pwa.planes)))
    return [(v[1], -v[2]) for v in verts]
end
