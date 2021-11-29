# Input types
abstract type Curvature end
struct Concave <: Curvature end
struct Convex <: Curvature end

abstract type Algorithm end
struct Heuristic <: Algorithm end
struct Optimized <: Algorithm end

struct FunctionEvaluations{D}
	points::Vector{NTuple{D}}
	values
end

# Skeleton of key methods (to be moved and renamed)
compute_approximation(input, c::Concave, a) = "flip curvature and call for convex"
compute_approximation(input, c::Convex, a::Optimized) = "optimize here"
# If we want to use dispatch for specializing on dimensions (maybe just do branching and call specialized function)
compute_approximation(input::FunctionEvaluations{D}, c::Convex, a::Heuristic) where D = compute_approximation(input, c, a, Val(D))
compute_approximation(input::FunctionEvaluations{D}, c::Convex, a::Heuristic, ::Val{1}) where D = "specialized on 1D"
compute_approximation(input::FunctionEvaluations{D}, c::Convex, a::Heuristic, dims) where D = "other Ds"

# Alternative:
# struct FunctionEvaluations{K,M,D}
# 	points::Vector{NTuple}
# 	values
# end
# compute_approximation(input::FunctionEvaluations{K,M,D}) where {K<:Concave,M,D} = "Concave - flip everything here and call for convex"
# compute_approximation(input::FunctionEvaluations{K,M,D}) where {K<:Convex,M,D} = "Convex"
# compute_approximation(input::FunctionEvaluations{K,M,D}) where {K<:Convex, M<:Optimized,D} = compute_approximation(input, D)
# compute_approximation(input::FunctionEvaluations{K,M,D},dim) where {K<:Convex,M<:Optimized,D} = "$dim"
# test = FunctionEvaluations{Convex,Optimized,3}([(1,2,3),(3,4,5)],[1,1])

# Result types
struct Plane{D,T<:Number}
    α::NTuple{D,T}
    β::T
end
Plane(a::NTuple{N,T}, b::T) where {N,T} = Plane{N,T}(a,b)
Plane(a::Vector{T}, b::T) where T = Plane(Tuple(a),b)

struct ConvexPWLFunc{D,T<:Number}
    planes::Vector{Plane{D,T}}
    M::T
end
ConvexPWLFunc(planes::Vector, M::T) where T = ConvexPWLFunc{length(planes), T}(planes, M)


# Old types
struct PWLFunction
    x::Vector{Float64}
    z::Vector{Float64}

    function PWLFunction(x::Vector{Float64}, z::Vector{Float64})
        new(x,z)
    end
end

function PWLFunction(x::Vector, z::Vector) 
    @assert length(x) == length(z)
    return PWLFunction(convert(Vector{Float64}, x),convert(Vector{Float64}, z))
end

struct ConvexPWLFunction   
    x::Vector{Float64}
    z::Vector{Float64}

    c::Vector{Float64}
    d::Vector{Float64}

    function ConvexPWLFunction(x::Vector{Float64}, z::Vector{Float64},c::Vector{Float64}, d::Vector{Float64})
        new(x, z, c, d)
    end
end


function ConvexPWLFunction(x::Vector, z::Vector) 
    @assert length(x) == length(z)
    return with_outer_rep(ConvexPWLFunction(convert(Vector{Float64}, x),convert(Vector{Float64}, z),Vector{Float64}(), Vector{Float64}()))
end

function ConvexPWLFunction(c::Vector, d::Vector, xmin, xmax) 
    @assert length(c) == length(d)
    return with_inner_rep(ConvexPWLFunction(Vector{Float64}(), Vector{Float64}(), convert(Vector{Float64}, c), convert(Vector{Float64}, d)), xmin, xmax)
end

function Base.print(io::IO, pwl::ConvexPWLFunction)
    Printf.@printf("    x        z\n")
    for i in 1:length(pwl.x)
        Printf.@printf("%8.2f %8.2f\n", pwl.x[i], pwl.z[i])
    end
end

function with_outer_rep(pwl::ConvexPWLFunction)
    
    𝑁  = length(pwl.x)
    𝑥 = pwl.x 
    𝑦 = pwl.z

    𝑐 = [(𝑦[i+1] - 𝑦[i]) / (𝑥[i+1] - 𝑥[i]) for i ∈ 1:(𝑁-1)]
    𝑑 = [𝑦[i] - 𝑐[i]*𝑥[i] for i ∈ 1:(𝑁-1)]

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
    for k in 1:K-1
        xi = (sseg[k+1].second - sseg[k].second) / (sseg[k].first - sseg[k+1].first)
        push!(x, xi)
        push!(y, sseg[k].first * xi + sseg[k].second)
    end
    push!(x, xmax)
    push!(y, sseg[K].first * xmax + sseg[K].second)

    return ConvexPWLFunction(x, y, pwl.c, pwl.d)
end

function evaluate(pwl::ConvexPWLFunction, x)
    return maximum(pwl.c[i]*x + pwl.d[i] for i=1:length(pwl.c))
end


function ConvexPWLFunction(x, fz::Function)
    @assert issorted(x)
    return ConvexPWLFunction(convert(Vector{Float64}, x), map(t->convert(Float64,fz(t)), x))
end

struct ConcavePWLFunction
   
    pwl::ConvexPWLFunction
   
    function ConcavePWLFunction(pwl::ConvexPWLFunction)
        new(pwl)
    end
end



function ConcavePWLFunction(x::Vector, z::Vector)             
    return ConcavePWLFunction(ConvexPWLFunction(x,-z))
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
        new(x, z, a, b)
    end
end

function ConvexPWLFunctionND(x, fz::Function)
    ##TODO: create function to calculate coefficients for given points and function values
    return ConvexPWLFunctionND(x, map(t->convert(Float64, fz(t)), x), (), ())
end

function ConvexPWLFunctionND(a::Vector{Tuple{Float64, Float64}}, b::Vector{Float64})
    @assert length(a) == length(b)
    return ConvexPWLFunctionND(Vector{Tuple}(), Vector{Float64}(), a, b)
end

struct ConcavePWLFunctionND
   
    pwl::ConvexPWLFunctionND
   
    function ConcavePWLFunctionND(pwl::ConvexPWLFunctionND)
        new(pwl)
    end
end

function ConcavePWLFunctionND(x::Vector, z::Vector)             
    return ConcavePWLFunctionNd(ConvexPWLFunctionNd(x,-z))
end

function evaluate(pwl::ConcavePWLFunctionND, x)        
    return -evaluate(pwl.pwl, x)
end

function Base.print(io::IO, pwl::ConvexPWLFunctionND)
    dim = length(pwl.a[1])
    for i ∈ 1:dim
        Printf.@printf("    c%d",i)
    end 
    Printf.@printf("    c%d\n",dim+1)

    for i ∈ 1:length(pwl.a)        
        for j ∈ 1:dim
            Printf.@printf("%8.2f", pwl.a[i][j])        
        end        
        
        Printf.@printf("%8.2f\n", pwl.b[i])
    end        
end

function evaluate(pwl::ConvexPWLFunctionND, x)
    return maximum(dot(pwl.a[i], x) + pwl.b[i] for i=1:length(pwl.a))
end



