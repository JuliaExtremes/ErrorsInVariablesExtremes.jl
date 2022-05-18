# Stationary

```@setup chaudiere2
using ErrorsInVariablesExtremes, Extremes, Dates, DataFrames, Distributions, Gadfly, Cairo, Fontconfig
```

## Load the data

Loading the pseudo-observations into an object of type [`Pseudoensemble`](@ref):
```@example chaudiere2
filename = "../../../test/data/A2020_Analyse_Historique_QMA_SLSO00003.nc"
pensemble = ErrorsInVariablesExtremes.load_discharge_distribution(filename)
```

The vector of years can be extracted:
```@example chaudiere2
years = pensemble.value[1].year;
```

## MCMC

```@example chaudiere2
res = gevfitbayes(pensemble)
```

## Display of annual maximum estimates



```@example chaudiere2

# Compute the estimates
ŷ = vec(mean( res.maxima[:, :, 1].value, dims=1))

# Compute the 95% credible interval
ymin = [quantile(vec(res.maxima[:, j, 1].value), .025) for j in 1:length(ŷ)]
ymax= [quantile(vec(res.maxima[:, j, 1].value), .975) for j in 1:length(ŷ)]

df = DataFrame(Year = years, Discharge=ŷ, ymin = ymin, ymax=ymax)

# plot the ponctual estimates and the intervals
set_default_plot_size(12cm, 8cm)
plot(df, x=:Year, y=:Discharge, Geom.line, Geom.point,
    ymin=ymin, ymax=ymax, Geom.ribbon
)
```

## Fit of the GEV law adjusted to the estimates of the annual maxima

```@example chaudiere2
set_default_plot_size(21cm ,16cm)
diagnosticplots(res.parameters)
```