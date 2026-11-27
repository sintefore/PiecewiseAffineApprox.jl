```@meta
CurrentModule = PiecewiseAffineApprox
```

## Algorithms

```@docs
MILP
Cluster
Interpol
Progressive
FullOrder
```

## Other types

```@docs
Convex
Concave
DC
Alternating
Select
Restart
Plane
FunctionEvaluations
PWAFunc
DCFunc
Quadratic
PWQFunc
domain
```

## Formulation methods and directions

```@docs
PWAMethod
PsiFormulation
LambdaFormulation
DeltaFormulation
PhiFormulation
Ψ_Formulation
λ_Formulation
Δ_Formulation
Φ_Formulation
Direction
Graph
Epigraph
Hypograph
```

## Extending

The interface for adding an algorithm, see [Add an algorithm](@ref).

```@docs
ApproximationAlgorithm
_fit
planes
concave_planes
metric
strict
optimizer
tolerance
tolerance_metric
with
supported_dimensions
requires_convex_data
supported_metrics
supported_strict
```

The interface for adding a formulation method, see [Add a formulation](@ref).

```@docs
formulate!
record!
```

## Functions

```@docs
approx
fit
FitResult
enforce_curvature
dc_decomposition
approximation_error
error_bounds
evaluate
add_pwaformulation
PWAFormulation
Formulated
```

## Deprecated

```@docs
pwaffine
```
