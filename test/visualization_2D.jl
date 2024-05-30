using Revise, GLMakie, PiecewiseLinearApprox, Xpress, JuMP

Opt = Xpress.Optimizer
PWA = PiecewiseAffineApprox

xg = [i for i ∈ -1:0.5:1]
yg = [j for j ∈ -1:0.5:1]
X = [repeat(xg, inner = [size(yg, 1)]) repeat(yg, outer = [size(xg, 1)])]

𝒫 = collect(Tuple(X'[:, i]) for i ∈ 1:size(X', 2))

z = X[:, 1] .^ 2 + X[:, 2] .^ 2
z_concave = z .* -1

f(u1, u2) = u1^2 + u2^2

planes = 4
dimensions = 2
pen = :l2

pwa1 = convex_linearization(𝒫, z, Opt; planes, dimensions, strict = :none, pen)
pwa2 = convex_linearization(X, z, Opt; planes, dimensions, strict = :above, pen)
pwa3 = concave_linearization(
    X,
    z_concave,
    Opt;
    planes,
    dimensions,
    strict = :above,
    pen,
)

ϵ = 1e-06
diff_1 = abs(PWA.evaluate(pwa1, (0.5, 0.5)) - f(0.5, 0.5)) / (f(0.5, 0.5) + ϵ)
diff_2 = abs(PWA.evaluate(pwa2, (0.5, 0.5)) - f(0.5, 0.5)) / (f(0.5, 0.5) + ϵ)

sc1 = plotconvND(pwa1, X, z)
sc2 = plotconvND(pwa2, X, z)

display(sc1)
display(sc2)
