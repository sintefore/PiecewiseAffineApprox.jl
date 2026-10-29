# Benchmark cases: synthetic functions with a known truth, and small real
# datasets used in the literature on piecewise and convex regression.

using PiecewiseAffineApprox
using DelimitedFiles, Downloads, SHA, StableRNGs, Statistics, LinearAlgebra

"""
A benchmark case: the data to fit, and optionally noise-free values at
independent test points for measuring the error against the true function.

- `planes`: the number of planes for methods with a fixed number of planes
  (the true number when known)
- `truek`: the true number of planes, if the truth is piecewise affine
- `tolerance`: the tolerance for methods that add planes until the
  approximation is within a tolerance (Progressive)
- `convexdata`: whether the data itself is convex (noise-free synthetic data),
  as required by the methods that interpolate the data (FullOrder and
  Progressive with `strict = :outer`)
"""
Base.@kwdef struct Case
    name::String
    group::Symbol
    curvature::Any = Convex()
    data::FunctionEvaluations
    test::Union{Nothing,FunctionEvaluations} = nothing
    planes::Int
    truek::Union{Nothing,Int} = nothing
    tolerance::Float64
    convexdata::Bool = false
    source::String = ""
end

dims(c::Case) = length(first(c.data.points))
npoints(c::Case) = length(c.data.values)

# The range of the true values (or of the data), used to report errors relative
# to the size of the function
function valuerange(c::Case)
    v = isnothing(c.test) ? c.data.values : c.test.values
    return maximum(v) - minimum(v)
end

_evals(X::AbstractMatrix, z) =
    FunctionEvaluations([Tuple(X[:, i]) for i ∈ axes(X, 2)], collect(z))

# --- Synthetic functions -------------------------------------------------------

# Specific fuel consumption of an engine as a function of the load fraction
sfc(l) = 155 + 12 / l + 60 * (l - 0.7)^2

# Random max-affine function of k planes in d dimensions on [-1, 1]^d, as in
# the max-affine regression literature (Ghosh et al., 2019): slopes and
# intercepts are standard normal. The planes are redrawn until each one is the
# maximum on at least a fraction 1/(4k) of the domain, so that k is the true
# number of planes.
function maxaffine(rng, d, k)
    U = 2 .* rand(rng, d, 10_000) .- 1
    for _ ∈ 1:1000
        A, b = randn(rng, k, d), randn(rng, k)
        active = [argmax(A * U[:, i] .+ b) for i ∈ axes(U, 2)]
        if all(count(==(j), active) ≥ size(U, 2) / (4k) for j ∈ 1:k)
            return x -> maximum(A * x .+ b)
        end
    end
    return error("No max-affine function with $k active planes found")
end

logsumexp(x) = (m = maximum(x); m + log(sum(exp.(x .- m))))

# Samples N uniform points in the box [lo, hi]^d, with values f(x) for a vector
# x plus normal noise with standard deviation σ times the range of f, and 2000
# noise-free test points.
function synthetic(name, f, d, N; lo = -1.0, hi = 1.0, σ = 0.0, seed = 1,
    planes, truek = nothing, curvature = Convex())
    rng = StableRNG(seed)
    box(n) = lo .+ (hi - lo) .* rand(rng, d, n)
    X, Xt = box(N), box(2000)
    ft = [f(Xt[:, i]) for i ∈ axes(Xt, 2)]
    r = maximum(ft) - minimum(ft)
    z = [f(X[:, i]) for i ∈ axes(X, 2)] .+ σ * r .* randn(rng, N)
    # Methods that add planes until the approximation is within a tolerance of
    # all points use the expected maximum of N normal deviations (plus a margin)
    # for noisy data, and 1 % of the range for noise-free data
    tolerance = σ > 0 ? (sqrt(2log(2N)) + 0.5) * σ * r : 0.01r
    return Case(; name, group = :synthetic, curvature, data = _evals(X, z),
        test = _evals(Xt, ft), planes, truek, tolerance, convexdata = σ == 0,
        source = "synthetic")
end

function synthetic_cases()
    cases = Case[]
    for σ ∈ (0.0, 0.02)
        tag = σ == 0 ? "exact" : "noise2%"
        push!(cases, synthetic("sfc-1d-$tag", x -> sfc(x[1]), 1, 100; lo = 0.2, hi = 1.0,
            σ, planes = 5))
        push!(cases, synthetic("exp-1d-$tag", x -> exp(x[1]), 1, 100; σ, planes = 5))
        f3 = maxaffine(StableRNG(3), 1, 3)
        push!(cases, synthetic("maxaff-1d-k3-$tag", f3, 1, 100; σ, planes = 3,
            truek = 3))
        f2 = maxaffine(StableRNG(4), 2, 4)
        push!(cases, synthetic("maxaff-2d-k4-$tag", f2, 2, 200; σ, planes = 4,
            truek = 4))
        push!(cases, synthetic("lse-2d-$tag", logsumexp,
            2, 200; σ, planes = 6))
        f3d = maxaffine(StableRNG(5), 3, 5)
        push!(cases, synthetic("maxaff-3d-k5-$tag", f3d, 3, 300; σ,
            planes = 5, truek = 5))
    end
    return cases
end

# Non-convex functions, approximated as a difference of convex functions
function nonconvex_cases()
    cases = Case[]
    for σ ∈ (0.0, 0.02)
        tag = σ == 0 ? "exact" : "noise2%"
        push!(cases, synthetic("sin-1d-$tag", x -> sin(3x[1]), 1, 100; σ, planes = 4,
            curvature = DC()))
        # The difference of two random max-affine functions, with 3 planes in
        # each part
        g, h = maxaffine(StableRNG(7), 1, 3), maxaffine(StableRNG(8), 1, 3)
        push!(cases, synthetic("dc-maxaff-1d-k3-$tag", x -> g(x) - h(x), 1, 100; σ,
            planes = 3, truek = 6, curvature = DC()))
        push!(cases, synthetic("sincos-2d-$tag", x -> sin(3x[1]) * cos(2x[2]), 2, 200;
            σ, planes = 6, curvature = DC()))
    end
    # The data of these cases is never convex, whatever the noise
    return [Case(; (k => getfield(c, k) for k ∈ fieldnames(Case))...,
        group = :nonconvex, convexdata = false) for c ∈ cases]
end

# Larger instances for timing and memory
function scaling_cases()
    cases = Case[]
    for N ∈ (500, 1000, 2000)
        push!(cases, synthetic("sfc-1d-noise2%-N$N", x -> sfc(x[1]), 1, N; lo = 0.2,
            hi = 1.0, σ = 0.02, planes = 5))
    end
    f2 = maxaffine(StableRNG(4), 2, 4)
    for N ∈ (1000, 5000, 20_000)
        push!(cases, synthetic("maxaff-2d-k4-noise2%-N$N", f2, 2, N; σ = 0.02,
            planes = 4, truek = 4))
    end
    f5 = maxaffine(StableRNG(6), 5, 8)
    push!(cases, synthetic("maxaff-5d-k8-noise2%-N5000", f5, 5, 5000;
        σ = 0.02, planes = 8, truek = 8))
    return [Case(; (k => getfield(c, k) for k ∈ fieldnames(Case))...,
        group = :scaling) for c ∈ cases]
end

# --- Real datasets -------------------------------------------------------------

const DATA = joinpath(@__DIR__, "data")
const CACHE = joinpath(@__DIR__, "cache")

_readcsv(path) = readdlm(path, ',', Float64; header = true)[1]

# Bacon & Watts (1971): log height of the stagnant surface layer against log
# flow rate of water. A two-segment change-point model with a steeper second
# segment (concave); Carlin, Gelfand & Smith (1992) report a change point at
# about 0.03 with slopes -0.42 and -1.01.
function bacon_watts()
    A = _readcsv(joinpath(DATA, "bacon_watts.csv"))
    return Case(; name = "bacon-watts", group = :real, curvature = Concave(),
        data = _evals(A[:, 1:1]', A[:, 2]), planes = 2, tolerance = 0.04,
        source = "Bacon & Watts (1971), Biometrika 58(3)")
end

# Stopping distance (ft) against speed (mph) for cars in the 1920s (Ezekiel,
# 1930), the `cars` dataset of R. Roughly quadratic, with repeated speeds.
function cars()
    A = _readcsv(joinpath(DATA, "cars.csv"))
    return Case(; name = "cars", group = :real, data = _evals(A[:, 1:1]',
            A[:, 2]), planes = 3, tolerance = 50.0,
        source = "Ezekiel (1930), R datasets::cars")
end

# Total cost (TOTEX) of 89 Finnish electricity distribution firms against
# energy transmitted, network length and number of customers, a standard
# dataset for convex nonparametric least squares (Kuosmanen, 2012). The data
# is distributed with pyStoNED (GPL-3), so it is downloaded on first use rather
# than included here.
const FINNISH_URL = "https://raw.githubusercontent.com/ds2010/pyStoNED/" *
                    "18493242d787d0a55f4a40c0003cb23185de0b6b/" *
                    "pystoned/data/electricityFirms.csv"
const FINNISH_SHA256 =
    "4dd2494e72e95431d118d91fe9ca97d0793010852248baae3a0df5ce2c9102d0"

function finnish_electricity()
    path = joinpath(CACHE, "electricityFirms.csv")
    if !isfile(path)
        mkpath(CACHE)
        Downloads.download(FINNISH_URL, path)
    end
    if bytes2hex(open(sha256, path)) != FINNISH_SHA256
        rm(path)
        error("Checksum mismatch for $path, removed it")
    end
    # OPEX, CAPEX, TOTEX, Energy, Length, Customers, PerUndGr
    A = _readcsv(path)
    return Case(; name = "finnish-electricity", group = :real,
        data = _evals(A[:, 4:6]', A[:, 3]), planes = 4,
        tolerance = 0.05 * (maximum(A[:, 3]) - minimum(A[:, 3])),
        source = "Kuosmanen (2012), Energy Economics 34(6); via pyStoNED")
end

function real_cases()
    cases = [bacon_watts(), cars()]
    try
        push!(cases, finnish_electricity())
    catch e
        @warn "Skipping the Finnish electricity data" exception = e
    end
    return cases
end
