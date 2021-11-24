using Revise,GLMakie, PiecewiseLinearApprox, Xpress, JuMP

PWL = PiecewiseLinearApprox

xg = [i for i in -1:0.5:1]
yg = [j for j in -1:0.5:1]
X = [repeat(xg, inner=[size(yg,1)]) repeat(yg, outer=[size(xg,1)])]

𝒫 = collect(Tuple(X'[:,i]) for i in 1:size(X',2))

z = X[:,1].^2 + X[:,2].^2
f(u1, u2) = u1^2 + u2^2

# Use Xpress without output as the solver
Opt = optimizer_with_attributes(Xpress.Optimizer)

np = 4
pwl1 = convex_linearization(𝒫, z, Opt; nplanes=np, dimensions=2, strict=:none, pen=:l2)
pwl2 = convex_linearization(X, z, Opt; nplanes=np, dimensions=2, strict=:above, pen=:l2)

eps = 1e-06
diff_1 = abs(PWL.evaluate(pwl1, (0.5, 0.5)) - f(0.5,0.5))/f(0.5,0.5)+eps
diff_2 = abs(PWL.evaluate(pwl2, (0.5, 0.5)) - f(0.5,0.5))/f(0.5,0.5)+eps

sc1 = plotconvND(pwl1, X, z)
sc2 = plotconvND(pwl2, X, z)

display(sc1)
display(sc2)



