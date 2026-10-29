# Benchmark

Compares the approximation methods of PiecewiseAffineApprox on synthetic
functions with a known truth and on small real datasets from the literature on
piecewise and convex regression.

```
julia --project=benchmark -e "import Pkg; Pkg.instantiate()"
julia --project=benchmark -t auto benchmark/run.jl [--full] [--cases=REGEX] [--methods=REGEX]
```

The default run takes a few minutes, most of it spent in the MILP method, which
is stopped at its 30 s time limit (exclude it with `--methods="^(?!MILP)"`).
`--full` adds larger instances for timing and memory. Results are written to
`benchmark/results` as a CSV file and a markdown report, together with the
commit, Julia version and number of threads.

## Metrics

For each case and method:

- `planes`: the number of planes of the approximation
- `data_rms`, `data_max`: root mean square and maximum deviation from the data
- `test_rms`, `test_max`: the same against the noise-free truth at 2000
  independent test points (synthetic cases only)
- `time` (s) and `alloc_mib`: for runs below one second, the best of three
- `reproducible`: for randomized methods, whether a second run with the same
  seed gives the same planes

Errors are relative to the range of the function, so that cases are
comparable.

## Cases

Synthetic, each without noise and with normal noise of 2 % of the range:

| case | dims | points | truth |
|---|---|---|---|
| `sfc-1d` | 1 | 100 | engine specific fuel consumption, 155 + 12/l + 60(l − 0.7)² on [0.2, 1] |
| `exp-1d` | 1 | 100 | exp(x) on [−1, 1] |
| `maxaff-1d-k3` | 1 | 100 | random max of 3 affine functions |
| `maxaff-2d-k4` | 2 | 200 | random max of 4 affine functions |
| `lse-2d` | 2 | 200 | log-sum-exp |
| `maxaff-3d-k5` | 3 | 300 | random max of 5 affine functions |

The random max-affine functions follow the max-affine regression literature
(Ghosh et al., 2019): standard normal slopes and intercepts, with each plane
the maximum on a fair share of [−1, 1]^d. Scaling cases (`--full`): `sfc-1d`
with up to 2000 points, `maxaff-2d-k4` with up to 20 000 points and a
5D max of 8 planes with 5000 points.

Non-convex functions, approximated as a difference of convex functions with
`DC()`, each without noise and with 2 % noise; `planes` is the number of planes
in each part:

| case | dims | points | truth |
|---|---|---|---|
| `sin-1d` | 1 | 100 | sin(3x) on [−1, 1] |
| `dc-maxaff-1d-k3` | 1 | 100 | difference of two random max-affine functions of 3 planes |
| `sincos-2d` | 2 | 200 | sin(3x) cos(2y) |

Real data:

| case | dims | points | description |
|---|---|---|---|
| `bacon-watts` | 1 | 29 | log height of the stagnant surface layer against log flow rate (concave, two segments) |
| `cars` | 1 | 50 | stopping distance against speed, with repeated speeds (convex) |
| `finnish-electricity` | 3 | 89 | total cost against energy, network length and customers (convex) |

Methods that need a tolerance (Progressive) use, for noisy synthetic data, the
expected maximum of the noise, so that a convex function within the tolerance
exists, and 1 % of the range for noise-free data. The methods that interpolate
the data (FullOrder, Progressive with `strict = :outer`) only run on
noise-free data.

## Formulations

```
julia --project=benchmark benchmark/formulations.jl [--units=50,200,1000] [--repeats=2]
```

Compares the formulations of `pwaffine` in a dispatch model with one piecewise
affine cost function per unit: the Ψ, λ, Δ and Φ formulations for convex
costs, which are linear programs, and the Ψ-formulation with big-M selection
against the Δ- and Φ-formulations for non-convex costs, which need binary
variables. Reports build and solve times and the model sizes to
`benchmark/results`; each solve is limited to 120 s.

## Data sources

- `data/bacon_watts.csv`: D. W. Bacon and D. G. Watts (1971), Estimating the
  transition between two intersecting straight lines, *Biometrika* 58(3),
  525–534. Values as listed in the SAS/STAT documentation of PROC MCMC
  (change point example), which also reports the fit of Carlin, Gelfand and
  Smith (1992): change point 0.028, slopes −0.42 and −1.01.
- `data/cars.csv`: M. Ezekiel (1930), *Methods of Correlation Analysis*,
  Wiley; the `cars` dataset of R.
- Finnish electricity distribution firms: T. Kuosmanen (2012), Stochastic
  semi-nonparametric frontier estimation of electricity distribution networks,
  *Energy Economics* 34(6), 2189–2199. Distributed with
  [pyStoNED](https://github.com/ds2010/pyStoNED) under GPL-3, and therefore
  downloaded on first use (to `benchmark/cache`, checked by SHA-256) rather
  than included.
