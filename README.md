# Extreme value analysis for errors-in-variables models


[![WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Build status](https://github.com/JuliaExtremes/ErrorsInVariablesExtremes.jl/workflows/CI/badge.svg)](https://github.com/JuliaExtremes/ErrorsInVariablesExtremes.jl/actions)
[![codecov](https://codecov.io/gh/JuliaExtremes/ErrorsInVariablesExtremes.jl/branch/master/graph/badge.svg?token=7UGVMF0ENE)](https://codecov.io/gh/JuliaExtremes/ErrorsInVariablesExtremes.jl)
[![documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaextremes.github.io/ErrorsInVariablesExtremes.jl)
[![documentation latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliaextremes.github.io/ErrorsInVariablesExtremes.jl/dev/)


*ErrorsInVariablesExtremes.jl* provides exhaustive high-performance functions for the analysis of extreme values in errors-in-variables models in Julia. In particular, methods for pseudodata models are implemented, such as:
* Parameter estimation by Bayesian procedure.
* Stationary and non-stationary models.
* Diagnostic plots for assessing model accuracy.
* Model selection with deviance information criterion (DIC). 


## Installation

The following **julia** command will install the package:

```julia
julia> import Pkg
julia> Pkg.add("https://github.com/JuliaExtremes/ErrorsInVariablesExtremes.jl")
```

See the [Package Documentation](https://juliaextremes.github.io/ErrorsInVariablesExtremes.jl/stable/) for details and examples.
