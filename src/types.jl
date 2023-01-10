# Input types
abstract type Curvature end
struct Concave <: Curvature end
struct Convex <: Curvature end

abstract type Algorithm end
struct Heuristic <: Algorithm end
struct Interpol <: Algorithm end
struct Optimized <: Algorithm end

struct FunctionEvaluations{D}
    points::Vector{<:NTuple{D,Number}}
    values::Vector{<:Number}
end

# Result types
struct Plane{D}
    α::NTuple{D,Number}
    β::Number
end
Plane(a::NTuple{N}, b) where {N} = Plane{N}(a, b)
Plane(a::Vector, b) = Plane(Tuple(a), b)
evaluate(p::Plane, x) = dot(p.α, x) + p.β

struct PWLFunc{C<:Curvature,D}
    planes::Vector{Plane{D}}
end
function PWLFunc(planes::Vector{Plane{D}}, C::Curvature) where {D}
    return PWLFunc{typeof(C),D}(planes)
end
PWLFunc{C,D}() where {C,D} = PWLFunc{C,D}(Vector{Plane{D}}())

function evaluate(pwl::PWLFunc{C,D}, x) where {C>:Concave,D}
    return -evaluate(PWLFunc{Convex,D}(pwl.planes), x)
end
function evaluate(pwl::PWLFunc{C,D}, x) where {C<:Convex,D}
    return maximum(dot(p.α, x) + p.β for p ∈ pwl.planes)
end
active(pwl, x) = argmax(collect(evaluate(p, x) for p ∈ pwl.planes))
planes(pwl) = length(pwl.planes)
addplane!(pwl::PWLFunc{C,D}, p::Plane{D}) where {C,D} = push!(pwl.planes, p)
addplane!(pwl::PWLFunc, α, β) = push!(pwl.planes, Plane(α, β))

# Old types
struct PWLFunction
    x::Vector{Float64}
    z::Vector{Float64}

    function PWLFunction(x::Vector{Float64}, z::Vector{Float64})
        return new(x, z)
    end
end

function PWLFunction(x::Vector, z::Vector)
    @assert length(x) == length(z)
    return PWLFunction(convert(Vector{Float64}, x), convert(Vector{Float64}, z))
end

struct ConvexPWLFunction
    x::Vector{Float64}
    z::Vector{Float64}

    c::Vector{Float64}
    d::Vector{Float64}

    function ConvexPWLFunction(
        x::Vector{Float64},
        z::Vector{Float64},
        c::Vector{Float64},
        d::Vector{Float64},
    )
        return new(x, z, c, d)
    end
end

function ConvexPWLFunction(x::Vector, z::Vector)
    @assert length(x) == length(z)
    return with_outer_rep(
        ConvexPWLFunction(
            convert(Vector{Float64}, x),
            convert(Vector{Float64}, z),
            Vector{Float64}(),
            Vector{Float64}(),
        ),
    )
end

function ConvexPWLFunction(c::Vector, d::Vector, xmin, xmax)
    @assert length(c) == length(d)
    return with_inner_rep(
        ConvexPWLFunction(
            Vector{Float64}(),
            Vector{Float64}(),
            convert(Vector{Float64}, c),
            convert(Vector{Float64}, d),
        ),
        xmin,
        xmax,
    )
end

function Base.print(io::IO, pwl::ConvexPWLFunction)
    Printf.@printf("    x        z\n")
    for i ∈ 1:length(pwl.x)
        Printf.@printf("%8.2f %8.2f\n", pwl.x[i], pwl.z[i])
    end
end

function with_outer_rep(pwl::ConvexPWLFunction)
    𝑁 = length(pwl.x)
    𝑥 = pwl.x
    𝑦 = pwl.z

    𝑐 = [(𝑦[i+1] - 𝑦[i]) / (𝑥[i+1] - 𝑥[i]) for i ∈ 1:(𝑁-1)]
    𝑑 = [𝑦[i] - 𝑐[i] * 𝑥[i] for i ∈ 1:(𝑁-1)]

    return ConvexPWLFunction(𝑥, 𝑦, 𝑐, 𝑑)
end

function with_inner_rep(pwl::ConvexPWLFunction, xmin, xmax)
    seg = Dict(zip(pwl.c, pwl.d))
    sseg = sort(collect(seg))

    K = length(sseg)

    x = Vector{Float64}()
    y = Vector{Float64}()

    push!(x, xmin)
    push!(y, sseg[1].first * xmin + sseg[1].second)
    for k ∈ 1:K-1
        xi =
            (sseg[k+1].second - sseg[k].second) /
            (sseg[k].first - sseg[k+1].first)
        push!(x, xi)
        push!(y, sseg[k].first * xi + sseg[k].second)
    end
    push!(x, xmax)
    push!(y, sseg[K].first * xmax + sseg[K].second)

    return ConvexPWLFunction(x, y, pwl.c, pwl.d)
end

function evaluate(pwl::ConvexPWLFunction, x)
    return maximum(pwl.c[i] * x + pwl.d[i] for i ∈ 1:length(pwl.c))
end

function ConvexPWLFunction(x, fz::Function)
    @assert issorted(x)
    return ConvexPWLFunction(
        convert(Vector{Float64}, x),
        map(t -> convert(Float64, fz(t)), x),
    )
end

struct ConcavePWLFunction
    pwl::ConvexPWLFunction

    function ConcavePWLFunction(pwl::ConvexPWLFunction)
        return new(pwl)
    end
end

function ConcavePWLFunction(x::Vector, z::Vector)
    return ConcavePWLFunction(ConvexPWLFunction(x, -z))
end

function evaluate(pwl::ConcavePWLFunction, x)
    return -evaluate(pwl.pwl, x)
end

struct ConvexPWLFunctionND
    x::Vector{Tuple}
    z::Vector{Float64}

    a::Vector{Tuple}
    b::Vector{Float64}

    function ConvexPWLFunctionND(x, z::Vector{Float64}, a, b::Vector{Float64})
        return new(x, z, a, b)
    end
end

function ConvexPWLFunctionND(x, fz::Function)
    ##TODO: create function to calculate coefficients for given points and function values
    return ConvexPWLFunctionND(x, map(t -> convert(Float64, fz(t)), x), (), ())
end

function ConvexPWLFunctionND(
    a::Vector{Tuple{Float64,Float64}},
    b::Vector{Float64},
)
    @assert length(a) == length(b)
    return ConvexPWLFunctionND(Vector{Tuple}(), Vector{Float64}(), a, b)
end

struct ConcavePWLFunctionND
    pwl::ConvexPWLFunctionND

    function ConcavePWLFunctionND(pwl::ConvexPWLFunctionND)
        return new(pwl)
    end
end

function ConcavePWLFunctionND(x::Vector, z::Vector)
    return ConcavePWLFunctionNd(ConvexPWLFunctionNd(x, -z))
end

function evaluate(pwl::ConcavePWLFunctionND, x)
    return -evaluate(pwl.pwl, x)
end

function Base.print(io::IO, pwl::ConvexPWLFunctionND)
    dim = length(pwl.a[1])
    for i ∈ 1:dim
        Printf.@printf("    c%d", i)
    end
    Printf.@printf("    c%d\n", dim + 1)

    for i ∈ 1:length(pwl.a)
        for j ∈ 1:dim
            Printf.@printf("%8.2f", pwl.a[i][j])
        end

        Printf.@printf("%8.2f\n", pwl.b[i])
    end
end

function evaluate(pwl::ConvexPWLFunctionND, x)
    return maximum(dot(pwl.a[i], x) + pwl.b[i] for i ∈ 1:length(pwl.a))
end
