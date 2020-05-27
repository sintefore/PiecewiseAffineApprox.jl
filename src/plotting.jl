function plot!(p::Plots.Plot, pwl::ConvexPWLFunction)
    xmin = minimum(pwl.x)
    xmax = maximum(pwl.x)
    zmin = minimum(pwl.z)
    zmax = maximum(pwl.z)
    
    for k in 1:length(pwl.c) 
        Plots.plot!(p, [xmin,xmax], [pwl.d[k] + pwl.c[k] * xmin , pwl.d[k] + pwl.c[k] * xmax], legend=:none, ylims=(zmin-0.5,zmax+0.5), linestyle = :dot)
    end
    Plots.plot!(p, pwl.x, pwl.z, linecolor=:black, ylims=(zmin-0.5,zmax+0.5), legend=:none)
    return p
end

function plot(pwl::ConvexPWLFunction)
    plot!(Plots.plot(), pwl)
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
    Plots.plot!(p, pwl.x, -pwl.z, linecolor=:black, ylims=(zmin-0.5,zmax+0.5))
    return p
end

function plot(pwl::ConcavePWLFunction)
    plot!(Plots.plot(), pwl)
end