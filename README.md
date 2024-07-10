# PiecewiseAffineApprox

[![Build Status](https://github.com/sintefore/PiecewiseAffineApprox.jl/workflows/CI/badge.svg?branch=main)](https://github.com/sintefore/PiecewiseAffineApprox.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/sintefore/PiecewiseAffineApprox.jl/branch/main/graph/badge.svg?token=2LXGVU04YS)](https://codecov.io/gh/sintefore/PiecewiseAffineApprox.jl)
[![In Development](https://img.shields.io/badge/docs-dev-blue.svg)](https://sintefore.github.io/PiecewiseAffineApprox.jl/dev/)



Add convex (or concave) piecewise linear approximations of functions or a set of points to optimization models modelled in [JuMP](https://jump.dev/). 

This package provides three main methods to fit a set of points: 

1. creates and solves a MILP to fit a set of points, and adds the resulting linear constraints to the optimization model. This method is partially based on [Toriello & Vielma, 2012](https://doi.org/10.1016/j.ejor.2011.12.030). 
2. uses a heuristic to fit the set of points. This method is based on [Magnani & Boyd, 2009](https://doi.org/10.1007/s11081-008-9045-3).
3. a progressive heuristic to add planes until a certain accuracy is met. This method is based on [Kazda & Li, 2024](https://doi.org/10.1016/j.ejor.2023.07.026)

For non-convex functions, consider using [PiecewiseLinearOpt.jl](https://github.com/joehuchette/PiecewiseLinearOpt.jl).

## Usage

```julia
using JuMP, PiecewiseAffineApprox, HiGHS

m = Model(HiGHS.Optimizer)
@variable(m, x)
# Create a piecewise linear approximation to x^2 on the interval [-1, 1]
pwa = approx(x -> x[1]^2, [(-1, 1)], Convex(), MILP(optimizer = HiGHS.Optimizer, planes=5))
# Add the pwa function to the model
z = pwaffine(m, x, pwa)
# Minimize
@objective(m, Min, z)
# Check approximation/solution at x = 0.5
@constraint(m, x >= 0.5)
optimize!(m)
value(z) # 0.2653
```

To keep dependencies light, PiecewiseAffineApprox does not include plotting by default. If the `Makie` or `Plots` package is loaded
before using the module, some simple plotting routines will be available

The following demonstrates how this can be achieved:

```julia
using PiecewiseAffineApprox, GLMakie, HiGHS

x = LinRange(0, 1, 20)
f(x) = first(x)^2
pwa = approx(f, [(0, 1)], Convex(), MILP(optimizer = HiGHS.Optimizer, planes = 3))
p = plot(x, f.(x), pwa)

using CairoMakie
save("approx.svg", p; backend=CairoMakie)
```
![](docs/approx.svg)

Animation showing the accuracy when adding more cuts:
Approximation of 3D function

 <video loop src="docs/rotation.mp4">  video </video> 

![](docs/approx_3D.png)
