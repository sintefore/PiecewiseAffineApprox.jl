
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

function plot!(p::Plots.Plot, pwl::ConvexPWLFunction)
    xmin = minimum(pwl.x)
    xmax = maximum(pwl.x)
    zmin = minimum(pwl.z)
    zmax = maximum(pwl.z)
    
    for k in 1:length(pwl.c) 
        Plots.plot!(p, [xmin,xmax], [pwl.d[k] + pwl.c[k] * xmin , pwl.d[k] + pwl.c[k] * xmax], legend=:none, ylims=(zmin-0.5,zmax+0.5), linestyle = :dot)
    end
    Plots.plot!(p, pwl.x, pwl.z, ylims=(zmin-0.5,zmax+0.5), legend=:none)
    return p
end

function plot(pwl::ConvexPWLFunction)
    plot!(plot(), pwl)
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

function plot!(p::Plots.Plot, concave::ConcavePWLFunction)

    pwl = concave.pwl
    
    xmin = minimum(pwl.x)
    xmax = maximum(pwl.x)
    zmin = -maximum(pwl.z)
    zmax = -minimum(pwl.z)
    
    for k in 1:length(pwl.c) 
        Plots.plot!(p, [xmin,xmax], [-pwl.d[k] - pwl.c[k] * xmin , -pwl.d[k] - pwl.c[k] * xmax], ylims=(zmin-0.5,zmax+0.5), linestyle = :dot)
    end
    Plots.plot!(p, pwl.x, -pwl.z, ylims=(zmin-0.5,zmax+0.5))
    return p
end
