# Approximate points
For piecewise-affine approximations obtained from a dataset, the function FunctionEvaluations(𝒳, 𝒵) can be used to create a structure that contains a tuple with the points in the domain 𝒳 and their corresponding function values 𝒵. The function approx will receive this struct as input to calculate the approximations that best fit the dataset.

## Unidimensional funcions (1D)

The following code creates a 2D full-order piecewise affine approximation for a tuple (x,z) where x is a collection of points with their corresponding function values stores in z. The selected approximation is convex, uses 5 planes (linear segments), and applies l₁-norm regulazation as the error measure.

```jldoctest
using PiecewiseAffineApprox, JuMP, HiGHS
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)

x = collect(range(-1, 1; length = 10))
z = x .^ 2

pwa = approx(
    FunctionEvaluations(Tuple.(x), z),
    Convex(),
    MILP(;optimizer, metric = :l1, planes = 5),
)

length(pwa.planes)

# output
5
```

You may find it more convenient to use Matrix input for the data points rather than wrapping the data in the FunctionEvaluations, consider the alternative version of the example above:

```jldoctest
using PiecewiseAffineApprox, JuMP, HiGHS
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)

x = collect(range(-1, 1; length = 10))
z = x .^ 2
m = hcat(x, z)' # Collect data in a single Matrix

pwa = approx(
    m,
    Convex(),
    MILP(;optimizer, metric = :l1, planes = 5),
)

length(pwa.planes)
# output
5
```

# Bi-variate functions (3D)
A dataset can also be used as input for 3D piecewise affine approximations. The following code creates a uniformly sampled domain X within (-1,1) and calculates the corresponding function values for a concave function z_concave.

```jldoctest
using PiecewiseAffineApprox, JuMP, HiGHS
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)

xg = [i for i ∈ -1:0.5:1]
yg = [j for j ∈ -1:0.5:1]

X = [repeat(xg, inner = [size(yg, 1)]) repeat(yg, outer = [size(xg, 1)])]

z = X[:, 1] .^ 2 + X[:, 2] .^ 2
z_concave = z .* -1

np = 17

pwa = approx(
    FunctionEvaluations(tuple.(eachcol(X)...), z),
    Convex(),
    MILP(
        ;optimizer,
        planes = np,
        strict = :outer,
        metric = :l1,
    ),
)
length(pwa.planes)

# output
17
```

Similarly, the data can be passed as a Matrix if that is more convenient:

```jldoctest
using PiecewiseAffineApprox, JuMP, HiGHS
optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent()=>true)

xg = [i for i ∈ -1:0.5:1]
yg = [j for j ∈ -1:0.5:1]

X = [repeat(xg, inner = [size(yg, 1)]) repeat(yg, outer = [size(xg, 1)])]

z = X[:, 1] .^ 2 + X[:, 2] .^ 2

m = hcat(X, z)'

np = 17

pwa = approx(
    m,
    Convex(),
    MILP(
        ;optimizer,
        planes = np,
        strict = :outer,
        metric = :l1,
    ),
)
length(pwa.planes)

# output
17
```




