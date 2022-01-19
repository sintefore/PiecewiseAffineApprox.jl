# big-M for a 1D function
x = [i for i in -1:0.1:1]
z = x.^2

M₁ = PWL.conv_linear_big_M(x,z) 
M₁⁺ = PWL.linear_big_M(x,z)

pwl₁ = approx(FunctionEvaluations(Tuple.(x), z), Convex(), Optimized(), ;optimizer, planes=5, bigM=:conv_linear_big_M)
pwl₁⁺ = approx(FunctionEvaluations(Tuple.(x), z), Convex(), Optimized(), ;optimizer, planes=5, bigM=:linear_big_M)

@test isapprox(PWL.evaluate(pwl₁, 0.4), 0.16, atol=0.015)
@test isapprox(PWL.evaluate(pwl₁⁺, 0.4), 0.16, atol=0.015) 

@test isless(M₁⁺, M₁)  # new big-M is smaller than previous but gives a better approximation

# big-M for a 2D function
xg = [i for i in -1:0.5:1]
yg = [j for j in -1:0.5:1]

X = [repeat(xg, inner=[size(yg,1)]) repeat(yg, outer=[size(xg,1)])]

z = X[:,1].^2 + X[:,2].^2
np = 5

pwl₂ = approx(FunctionEvaluations(PWL.mat2tuples(X), z), Convex(), Optimized() ;optimizer=quadopt(1.7), planes=np, dimensions=2, strict=:above, pen=:l2, bigM = :conv_linear_big_M_ND) # old estimate for big-M

pwl₂⁺ = approx(FunctionEvaluations(PWL.mat2tuples(X), z), Convex(), Optimized() ;optimizer=quadopt(1.7), planes=np, dimensions=2, strict=:above, pen=:l2, bigM = :linear_big_M)

@test isapprox(PWL.evaluate(pwl₂, (0.5, 0.5)), PWL.evaluate(pwl₂⁺, (0.5, 0.5)), atol=0.01)

M₂⁺ = PWL.linear_big_M(x,z)
M₂ = PWL.conv_linear_big_M_ND(x,z)  # old estimate for big-M

@test isless(M₂, M₂⁺)