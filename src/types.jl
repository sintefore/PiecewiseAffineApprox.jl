
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

    Mᵇⁱᵍ::Float64

    function ConvexPWLFunction(x::Vector{Float64}, z::Vector{Float64},c::Vector{Float64}, d::Vector{Float64}; Mᵇⁱᵍ::Float64=nothing)        
        new(x, z, c, d,Mᵇⁱᵍ )        
    end
end


function ConvexPWLFunction(x::Vector, z::Vector, Mᵇⁱᵍ::Float64 = nothing) 
    @assert length(x) == length(z)
    return with_outer_rep(ConvexPWLFunction(convert(Vector{Float64}, x),convert(Vector{Float64}, z),Vector{Float64}(), Vector{Float64}(), Mᵇⁱᵍ))
end

function ConvexPWLFunction(c::Vector, d::Vector, xmin, xmax, Mᵇⁱᵍ::Float64 = nothing) 
    @assert length(c) == length(d)
    return with_inner_rep(ConvexPWLFunction(Vector{Float64}(), Vector{Float64}(), convert(Vector{Float64}, c), convert(Vector{Float64}, d), Mᵇⁱᵍ), xmin, xmax)
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

function evaluate(concave::ConcavePWLFunction, x)
    return -evaluate(concave.pwl, x)
end

struct Convex2dPWLFunction
   
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}

    c::Vector{Float64}
    d::Vector{Float64}
    e::Vector{Float64}

    Mᵇⁱᵍ::Float64

    function Convex2dPWLFunction(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64},c::Vector{Float64}, d::Vector{Float64}, e::Vector{Float64}, Mᵇⁱᵍ::Float64=nothing)
        new(x, y, z, c, d, e, Mᵇⁱᵍ)
    end
end

function Convex2dPWLFunction(x, y, fz::Function)
    @assert issorted(x)
    @assert issorted(y)

    ##TODO: create function to calculate coefficients for given points and function values
    return Convex2dPWLFunction(convert(Vector{Float64}, x), convert(Vector{Float64}, y), map(t->convert(Float64,fz(t)), x), Vector{Float64}(), Vector{Float64}(), Vector{Float64}())
end

function Convex2dPWLFunction(c::Vector, d::Vector, e::Vector, Mᵇⁱᵍ::Float64=nothing) 
    @assert length(c) == length(d) == length(e)
    return Convex2dPWLFunction(Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), convert(Vector{Float64}, c), convert(Vector{Float64}, d), convert(Vector{Float64}, e), Mᵇⁱᵍ)
end


function Base.print(io::IO, pwl::Convex2dPWLFunction)
    Printf.@printf("    x        y        z\n")
    for i ∈ 1:length(pwl.x), j ∈ 1:length(pwl.y)
        Printf.@printf("%8.2f %8.2f %8.2f\n", pwl.x[i,j], pwl.y[i,j], pwl.z[i,j])
    end
end

function evaluate(pwl::Convex2dPWLFunction, x,y)
    
    return maximum(pwl.c[i]*x + pwl.d[i]*y + pwl.e[i] for i=1:length(pwl.c))
end




