PiecewiseAffineApprox release notes
===================================

Version 0.7.0 (unreleased)
--------------------------
Breaking:
- Remove the JSON3 extension and the dependency on StructTypes. Serialization
  is provided by the JSON extension, which reads and writes the same file format
- Deprecate `pwaffine`: add an approximation to a model by calling it, `f(x)`,
  or with `add_pwaformulation`
- `Cluster` uses `metric = :l2` by default
- Remove the `:knots` convexify option of `Interpol`

New features:
- Difference of convex approximations of non-convex functions with `DC()` and
  `DCFunc`, and `dc_decomposition`
- Common fitting interface with `ApproximationAlgorithm`, `fit` and `FitResult`
- `Select` to choose the number of planes by tolerance, AIC/BIC or cross
  validation, `Restart` for randomized algorithms and `Alternating`
- `approximation_error` and `error_bounds`, and one-sided approximations over
  the whole domain with `everywhere = true`
- Formulations `PsiFormulation`, `LambdaFormulation`, `DeltaFormulation` and
  `PhiFormulation` (with Greek aliases), recorded in a `PWAFormulation` that
  can be rolled back, and the `formulate!`/`record!` interface for new ones
- Fitted functions record their domain, which is used as default bounds and
  for plotting without data
- Vertex enumeration and plotting of 2D functions with their vertices
- Piecewise quadratic approximations, the maximum of convex quadratic pieces,
  with `Cluster(pieces = :quadratic)`, returned as a `PWQFunc` of `Quadratic`
  pieces and added to a model as rotated second order cones
- Precompilation workload

Improvements:
- `Interpol` uses exact dynamic programming for knot selection, with options
  to convexify the data and adjust values or breakpoints
- `Cluster` normalizes the data, has a closed-form least squares fit for `:l2`,
  an `rng` option for reproducible results and an optional optimizer
- `MILP` solves on normalized data with tighter big-M values and warm starts
- `Progressive` and `FullOrder` support two-sided and inner error bands, lazy
  constraint generation and results independent of the scale of the data
- Replace the dependency on Distributions with an internal sampler
- Rewritten documentation, with a manual, a page on each method and the
  formulations, and figures

Version 0.6.6 (2026-03-26)
--------------------------
- Support Makie 0.24 (#32)

Version 0.6.5 (2026-03-26)
--------------------------
- Serialization with JSON v1 (#27)

Version 0.6.4 (2025-06-23)
--------------------------
- Fix strictness for concave functions (#26)

Version 0.6.3 (2025-02-18)
--------------------------
- Serialization of approximations with StructTypes and a JSON3 extension
- Fix printing of concave functions
- Support newer Makie versions

Version 0.6.2 (2024-10-22)
--------------------------
- Fix hang in the full order fit of `Progressive`
- Log with `@debug` instead of `@info`

Version 0.6.1 (2024-07-14)
--------------------------
- Doctests for the documentation
- Fix plot of the active segment

Version 0.6.0 (2024-07-14)
--------------------------
- Documentation with Documenter
- Rename the algorithm types and the `pen`/`penalty` options to `metric`
- `enforce_curvature` to make data convex or concave
- Improved printing of function evaluations and approximations

Version 0.5.0 (2024-05-31)
--------------------------
- Using affine rather than linear in naming
- Improve BigM calculation
- Pass arguments with Algorithm structs in favor of kwargs
- Add MIT license

Version 0.4.0 (2024-01-05)
--------------------------
- Rename package to PiecewiseAffineApprox to avoid name collision in General registry

Version 0.3.4 (2022-04-07)
--------------------------
- Fix type of pwlinear parameter for JuMP-variables

Version 0.3.3 (2022-03-10)
--------------------------
- Update dependencies to JuMP 0.23/1.0 (#13)
- Fix most broken tests

Version 0.3.2 (2022-01-28)
--------------------------
- Workaround for bigM with normals containing near-zeros (#10)
- Workaround for solver issue with strict=:none for compressor duty (#10)

Version 0.3.1 (2022-01-24)
--------------------------
- Improved calculation of big-M (#9)

Version 0.3.0 (2021-12-16)
--------------------------
- Add support for 2D functions (#3)
- New and breaking typed interface (#7)
- Add support for approximations based on Magnani/Boyd paper (#8)


Version 0.2.0 (2020-05-27)
--------------------------
- Added interpolation method (#1)
- Added optional plotting (#2)
- Added support for concave functions (#2)
- Breaking interface changess

Version 0.1.0 (2020-05-04)
--------------------------
- Initial release
