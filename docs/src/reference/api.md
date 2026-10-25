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
domain
```

## Formulations

```@docs
PWAFormulation
PsiFormulation
LambdaFormulation
DeltaFormulation
PhiFormulation
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
with
supported_dimensions
requires_convex_data
```

## Functions

```@docs
approx
enforce_curvature
dc_decomposition
approximation_error
evaluate
pwaffine
```