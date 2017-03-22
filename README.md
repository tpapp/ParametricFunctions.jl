# ParametricFunctions

[![Build Status](https://travis-ci.org/tpapp/ParametricFunctions.jl.svg?branch=master)](https://travis-ci.org/tpapp/ParametricFunctions.jl)
[![Coverage Status](https://coveralls.io/repos/tpapp/ParametricFunctions.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/tpapp/ParametricFunctions.jl?branch=master)
[![codecov.io](http://codecov.io/github/tpapp/ParametricFunctions.jl/coverage.svg?branch=master)](http://codecov.io/github/tpapp/ParametricFunctions.jl?branch=master)
[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

A unified interface for parametric function families.

## Introduction

This is a package for handling parametric function families of the form `f(x; θ)`.

A particular `f(⋅,⋅)` is a *function family*. It can be queried about its `domain`, degrees of freedom `degf` (the dimension of `θ ∈ ℝⁿ`), evaluated at `x` and `θ`.

This package is mainly useful for solving functional problems using projection methods, ie given some mapping `F`, solve for `F(f(x;θ)) = 0`. The library should be transparent to automatic differentiation tools such as [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). 

The goal of this package is facilitating the management of function families for a mostly "manual" approach to constructing approximations to functions, and then minimizing then some kind of residual. The function families are built from a set of basic building blocks using univariate families (such as Chebyshev polynomials, splines), domain and value transformations (function families need not be linear), combined into tensored multivariate families, or more sophisticated variations such as complete polynomials or sparse approximations.

There is also limited support for calculating `θ` for a given `y = f(x; θ)`, where `x` are the collocation `points` of the function family, but this is mainly useful for initial guesses for `θ`.

## Plans

Currently I am refining the API for univariate functions, using simple function families. Multivariate constructs are coming later. This is work in progress, the API may change without notice.

## Related packages

If your problem is (mostly) linear, without embedded optimization problems in `F` above, you are probably better off with [ApproxFun.jl](https://github.com/JuliaApproximation/ApproxFun.jl).
