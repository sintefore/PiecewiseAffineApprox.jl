# Approximate points
Piecewise affine approximations can be obtained either for a given function 𝒻(x) or for a set of sampled points 𝒳 with their corresponding function values 𝒵 = {f(x): ∀x ∈ 𝒳}. Both methods have their applicability depending on the context and the specific application requirements. For example, the method that approximates a function based on a collection of points may be useful for applications that provide a dataset for off-line construction of functions of interest prior to optimization.

Once the datataset is available, the function FunctionEvaluations(𝒳, 𝒵) can be used to create a structure that contains a tuple with the points in the domain 𝒳 and their corresponding function values 𝒵. The function approx will receive this struct as input to calculate the approximations that best fit the dataset.

## 2D 

The following code creates a 2D full-order piecewise affine approximation for a tuple (x,z) where x is a collection of points with their corresponding function values stores in z. The selected approximation is convex, uses 5 planes (linear segments), and applies l₁-norm regulazation as the error measure.

```julia
using PiecewiseAffineApprox, HiGHS

x = collect(range(-1, 1; length = 10))
z = x .^ 2

pwa1 = approx(
    FunctionEvaluations(Tuple.(x), z),
    Convex(),
    Optimized(optimizer = HiGHS.Optimizer, pen = :l1, planes = 5),
)
```

# 3D
A dataset can also be used as input for 3D piecewise affine approximations. The following code creates a uniformly sampled domain X within (-1,1) and calculates the corresponding function values for a concave function z_concave.

```julia
using PiecewiseAffineApprox, HiGHS

xg = [i for i ∈ -1:0.5:1]
yg = [j for j ∈ -1:0.5:1]

X = [repeat(xg, inner = [size(yg, 1)]) repeat(yg, outer = [size(xg, 1)])]

z = X[:, 1] .^ 2 + X[:, 2] .^ 2
z_concave = z .* -1

np = 17

function mat2tuples(x::Matrix)
    return collect(Tuple(x'[:, i]) for i ∈ 1:size(x', 2))
end

pwa1 = approx(
    FunctionEvaluations(mat2tuples(X), z),
    Convex(),
    Optimized(
        optimizer = HiGHS.Optimizer,
        planes = np,
        strict = :outer,
        pen = :l1,
    ),
)
```
