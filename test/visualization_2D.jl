using Revise,GLMakie, PiecewiseLinearApprox, Xpress, JuMP

Opt = Xpress.Optimizer
PWL = PiecewiseLinearApprox

xg = [i for i in -1:0.5:1]
yg = [j for j in -1:0.5:1]
X = [repeat(xg, inner=[size(yg,1)]) repeat(yg, outer=[size(xg,1)])]

𝒫 = collect(Tuple(X'[:,i]) for i in 1:size(X',2))

z = X[:,1].^2 + X[:,2].^2
z_concave = z.*-1

f(u1, u2) = u1^2 + u2^2

nplanes = 4
dimensions = 2
pen = :l2

pwl1 = convex_linearization(𝒫, z, Opt; 
    nplanes, dimensions, strict=:none, pen)
pwl2 = convex_linearization(X, z, Opt; 
    nplanes, dimensions, strict=:above, pen)
pwl3 = concave_linearization(X, z_concave, Opt; 
    nplanes, dimensions, strict=:above, pen)

ϵ = 1e-06
diff_1 = abs(PWL.evaluate(pwl1, (0.5, 0.5)) - f(0.5, 0.5)) / (f(0.5, 0.5) + ϵ)
diff_2 = abs(PWL.evaluate(pwl2, (0.5, 0.5)) - f(0.5, 0.5)) / (f(0.5, 0.5) + ϵ)

sc1 = plotconvND(pwl1, X, z)
sc2 = plotconvND(pwl2, X, z)

display(sc1)
display(sc2)



