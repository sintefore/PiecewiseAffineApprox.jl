```@meta
CurrentModule = PiecewiseAffineApprox
```

## Algorithms

```@docs
#FullOrderFitting
#Heuristic
#Interpol
#Optimized
#ProgressiveFitting
```

## Other types

```@docs

#Convex     # Missing docstring
#Concave    # Missing docstring
DC
Alternating
Select
Plane
FunctionEvaluations
PWAFunc
DCFunc
```

## Formulations

```@docs
PWAFormulation
PsiFormulation
LambdaFormulation
DeltaFormulation
PhiFormulation
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