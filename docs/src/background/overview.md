# Overview

A convex piecewise affine function is the maximum of a set of planes,
``f(x) = \max_k \, (\alpha_k^\top x + \beta_k)``. Fitting one to data means
choosing the planes, and the hard part is deciding which points each plane
should fit: once that is known, each plane is an easy fit to its points. The
algorithms differ in how they make that choice.

| | [`MILP`](@ref) | [`Cluster`](@ref) | [`Interpol`](@ref) | [`Progressive`](@ref) | [`FullOrder`](@ref) |
|:--|:--|:--|:--|:--|:--|
| **Idea** | Optimal assignment of points to planes | Alternate grouping and fitting, from random starts | Knots among the data points: a shortest path with at most ``K`` arcs | Add planes until the tolerance is met | One plane per data point |
| **Result** | Optimal for the number of planes | Local optimum | Optimal for its knots | Within the tolerance | Interpolates convex data |
| **Variables** | Any | Any | One | Any | Any |
| **Solver** | MILP solver | None for least squares | Optional | LP solver | LP solver |
| **You choose** | `planes` | `planes` | `planes` | `tolerance` | — |
| **Default metric** | `:l1` | `:l2` | `:l1` | `:l1` | `:l1` |
| **`strict`** | `:none`, `:outer`, `:inner` | `:none`, `:outer`, `:inner` | — | `:outer` (default), `:none`, `:inner` | always below |

## One interface for all

Every algorithm only fits convex data with a convex function. The rest is
shared by all algorithms:

- **Concave data** is fitted by negating the values, fitting the convex data,
  and negating the planes.
- **Non-convex data** is split into a convex and a concave part, which are
  fitted separately and then refitted alternately, see
  [Difference of convex](@ref).
- **Functions** are sampled on a grid over a box.
- **Input** is checked: the number of points and variables, and the options
  `metric` and `strict`.
- **Wrappers** combine with any algorithm: [`Select`](@ref) chooses the number
  of planes, [`Restart`](@ref) reruns randomized algorithms and keeps the best,
  and [`Alternating`](@ref) fits the two parts of a difference of convex
  approximation with different algorithms.
- **Diagnostics** from [`fit`](@ref): the error, the run time and details from
  the algorithm.

A new algorithm only implements the convex fit and gets all of this, see
[Add an algorithm](@ref).

## Data in any units

The algorithms work on data scaled to the order of one, and map the result
back, with ties (e.g. a point where two planes meet) broken by fixed rules on
the scaled data. The approximation therefore does not depend on the units of
the inputs and the values. [`Cluster`](@ref) also gives the same result after
rotating or stretching the inputs.

## The methods

- [MILP](@ref): the optimal fit for a given number of planes, by mixed-integer
  programming.
- [Cluster](@ref): a fast heuristic for many points and any number of
  variables.
- [Interpol](@ref): exact and fast for one variable.
- [Progressive and FullOrder](@ref): set the accuracy, not the number of
  planes.
- [Difference of convex](@ref): non-convex functions as convex plus concave.
- [Formulations](@ref): how an approximation is added to a JuMP model.
