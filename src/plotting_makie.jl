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

