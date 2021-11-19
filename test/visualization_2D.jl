using PiecewiseLinearApprox, Xpress, JuMP, LazyGrids

PWL = PiecewiseLinearApprox

xg = [i for i in -1:0.1:1]
yg = [j for j in -1:0.1:1]
x = [xg, yg]

X = ndgrid(x[1], x[2])
z = X[1].^2 + X[1].*X[2] + X[2].^2 
f(u1, u2) = u1^2 + u1*u2 + u2^2
plane(a,b,c) = a.*X[1] + b.*X[2] + c.*ones(size(X[1]))

# Use Xpress without output as the solver
Opt = optimizer_with_attributes(Xpress.Optimizer)

np = 8
pwl1 = convex_linearization(x, z, Opt; nplanes=np, dimensions=2, pen=:l2)

z_hat = PWL.evaluate(pwl1, 0.5, 0.5)

using Polyhedra
using GLMakie

h = hrep([HalfSpace([pwl1.c[i], pwl1.d[i]], -pwl1.e[i]) for i in 1:length(pwl1.c)])         

p = [plane(pwl1.c[i], pwl1.d[i], pwl1.e[i]) for i in 1:length(pwl1.c)]

Makie.surface(xg, yg, p[1], shading=true, transparency=true)
for i in 2:length(pwl1.c)
    Makie.surface!(xg, yg, p[i], shading=true, transparency=true)
end        

Makie.wireframe!(xg, yg,z)
