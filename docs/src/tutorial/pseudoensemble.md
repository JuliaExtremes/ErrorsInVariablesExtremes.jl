# Pseudoensemble

The [`Pseudoensemble`](@ref) data structure is illustrated using the unobserved riverflows at an ungauged section of the Chaudière River in Québec from 1960 to 2020.

```@setup chaudiere
using ErrorsInVariablesExtremes, Extremes, Dates, DataFrames, Distributions, Gadfly
```

## Load the data

Loading the pseudo-observations into an object of type [`Pseudoensemble`](@ref):
```@example chaudiere
filename = "SLSO00003"
pensemble = ErrorsInVariablesExtremes.load_discharge_distribution(filename)
```

The vector of years can be extracted:
```@example chaudiere
years = pensemble[1].year;
```

## Convert to DataFrame

Converting `Pseudoensemble` to `DataFrame`:
```@example chaudiere
df = convert(DataFrame, pensemble)
first(df,5)
```

## Extracting the statistics
```@example chaudiere
# Interval confidence level
α = 0.5

# Median
df[:, :y] = median.(df.Distribution)

# Lower bound of the confidence interval
df[:, :ymin] = quantile.(df.Distribution, α/2)

# Upper bound of the confidence interval
df[:, :ymax] = quantile.(df.Distribution, 1-α/2)

first(df, 5)
```


## Display the discharge of one specific configuration


```@example chaudiere
config = "LN24HA"
df2 = filter(row -> row.Configuration ==config, df)

set_default_plot_size(12cm, 8cm)
plot(df2, x=:Year, y=:y, Geom.line, Geom.point,
    ymin=:ymin, ymax=:ymax, Geom.ribbon,
    Guide.ylabel("Discharge [m³/s]"),
    Guide.title(config))
```

## Display the discharge of all configurations

```@example chaudiere
set_default_plot_size(12cm, 8cm)
fig = plot(df, x=:Year, y=:y, color=:Configuration, Geom.line,
    ymin=:ymin, ymax=:ymax, Geom.ribbon,
    Guide.ylabel("Discharge [m³/s]"))
```