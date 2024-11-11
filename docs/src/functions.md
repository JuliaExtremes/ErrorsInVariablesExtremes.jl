# Functions

## Methods for Pseudodata

```@docs
convert(::Type{DataFrame}, ::Vector{Pseudodata})
ensemblemean(::Vector{Pseudodata})
```

## Parameter estimation

```@docs
fitbayes
```

## Methods for fitted models

```@docs
convert(::Type{MaximumLikelihoodEVA}, ::PseudoMaximaEVA, ::Int)
dic
```

## Diagnostic plots

```@docs
diagnosticplots
histplot
probplot
qqplot
returnlevelplot
```

## Types

```@autodocs
Modules = [ErrorsInVariablesExtremes]
Private = false
Order = [:type]
```
