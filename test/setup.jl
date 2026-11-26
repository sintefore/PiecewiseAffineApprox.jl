# Packages and definitions shared by the test files
using HiGHS
using JuMP
using JSON
using LinearAlgebra
using Logging
using PiecewiseAffineApprox
using PiecewiseAffineApprox: Direction, Graph, Epigraph, Hypograph
using Random
using StableRNGs
using Test
rng = StableRNG(123)

const PWA = PiecewiseAffineApprox
const optimizer =
    optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)
