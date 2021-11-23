using Revise, PiecewiseLinearApprox, Xpress, JuMP

PWL = PiecewiseLinearApprox

xg = [i for i in -1:0.5:1]
yg = [j for j in -1:0.5:1]
X = [repeat(xg, inner=[size(yg,1)]) repeat(yg, outer=[size(xg,1)])]

𝒫 = collect(Tuple(X'[:,i]) for i in 1:size(X',2))

z = X[:,1].^2 + X[:,2].^2

f(u1, u2) = u1^2 + u2^2
plane(a,b,c) = a.*X[:,1] + b.*X[:,2] + c.*ones(size(X[:,1]))

# Use Xpress without output as the solver
Opt = optimizer_with_attributes(Xpress.Optimizer)

np = 4
pwl1 = convex_linearization(𝒫, z, Opt; nplanes=np, dimensions=2, strict=false, pen=:l2)
pwl2 = convex_linearization(X, z, Opt; nplanes=np, dimensions=2, strict=false, pen=:l2)

z_hat = PWL.evaluate(pwl1, (0.5, 0.5))

a_coeff =[pwl1.a[i][j] for i ∈ 1:length(pwl1.a), j ∈ 1:length(pwl1.a[1])]
p = [plane(a_coeff[i,1], a_coeff[i,2], pwl1.b[i]) for i in 1:length(pwl1.a)]

Makie.surface(xg, yg, reshape(p[1], length(xg), length(yg)), shading=true, transparency=false)
for i in 2:length(pwl1.a)
    Makie.surface!(xg, yg, reshape(p[i], length(xg), length(yg)), shading=true, transparency=false)
end        

zb = reshape(z, length(xg), length(yg))
scatter!(xg,yg,zb, markersize=30, color=:blue)


