using .GLMakie


function plotconvND(pwl::ConcavePWLFunctionND, X, z)        
    plotconvND(pwl.pwl, X, z; sign=-1)
end

function plotconvND(pwl::ConvexPWLFunctionND, X, z; sign=1)    
    xg = unique(X[:,1])
    yg = unique(X[:,2])

    plane(a,b,c) = a.*X[:,1] + b.*X[:,2] + c.*ones(size(X[:,1]))
    a_coeff =[pwl.a[i][j] for i ∈ 1:length(pwl.a), j ∈ 1:length(pwl.a[1])]
    p = [plane(a_coeff[i,1].*sign, a_coeff[i,2].*sign, pwl.b[i].*sign) for i in 1:length(pwl.a)]

    scene = surface(xg, yg, reshape(p[1], length(xg), length(yg)), shading=true, transparency=false)
    for i in 2:length(pwl.a)
        surface!(xg, yg, reshape(p[i], length(xg), length(yg)), shading=true, transparency=false)
    end

    zb = reshape(z, length(xg), length(yg))
    scatter!(xg,yg,zb, markersize=30, color=:blue)    

    display(scene)
    return scene
end

function plotconv2D(x, z, pwl::PWLFunc{Convex,2})
    xmin = minimum(x[1,:])
    xmax = maximum(x[1,:])

    ymin = minimum(x[2,:])
    ymax = maximum(x[2,:])

    x̄ = LinRange(xmin, xmax, 20)
    ȳ = LinRange(ymin, ymax, 20)  
    fig = Figure(resolution=(1000,1000)) 
    ax1 = Axis3(fig[1:2,1:2])
    ax2 = Axis(fig[1,3])
    ax3 = Axis(fig[2,3])
    
    scatter!(ax1, x[1,:], x[2,:], z, color=:red, markersize = 2000)
    for p in pwl.planes
        f = [evaluate(p,[x̄[i],ȳ[j]]) for i in 1:length(x̄), j in 1:length(ȳ)] 
        surface!(ax1, x̄, ȳ, f)
    end
    l1 = approx_error(x, z, pwl, :l1)
    l2 = approx_error(x, z, pwl, :l2)
    lmax = approx_error(x, z, pwl, :max)
    ax1.title = "l1 = $(round(l1, digits=2)), l2 = $(round(l2, digits=2)), max = $(round(lmax, digits=2))"
    z̄ = [abs(evaluate(pwl, x[:,i]) - z[i]) for i in 1:length(z)]
    θ = 20 / maximum(z̄) 
    ax2.title = "Max error = $(round(maximum(z̄),digits=2))"
    scatter!(ax2, x[1,:], x[2,:]; markersize=θ*z̄)

    𝒫 = update_partition(x, pwl)
    for j in 1:length(𝒫)
        x̄ = x[:, 𝒫[j] ]
        scatter!(ax3, x̄[1,:], x̄[2,:], marker=:xcross)
    end

    display(fig)
end

function plotconv2D(input::FunctionEvaluations{2}, pwl::PWLFunc{Convex,2}) 
    x = [p[i] for i in 1:2, p in input.points]
    z = input.values
    plotconv2D(x,z,pwl)
end

