# Unobserved riverflows in Québec, Canada

The functionalities of *ErrorsInVariablesExtremes.jl* are illustrated by producing an Extreme Value Analysis of unobserved riverflows of a section of the Chaudière River in Québec, Canada, from 1960 to 2020.

```@setup chaudiere
using ErrorsInVariablesExtremes, Extremes, DataFrames, CSV, Distributions, Mamba, Gadfly, Cairo, Fontconfig
```

---
## Load the data

Loading the pseudo-observations into a vector of [`Pseudodata`](@ref):
```@example chaudiere
filename = "../../../test/data/A2020_Analyse_Historique_QMA_SLSO00003.nc"
pdata = load_discharge_distribution(filename)
```

The vector of years can be extracted:
```@example chaudiere
years = pdata[1].year;
```

### Convert to DataFrame

Converting the vector of `Pseudodata` to `DataFrame`:
```@example chaudiere
df = convert(DataFrame, pdata)
first(df,5)
```

### Compute the median and the 1-α confidence interval for each of the configurations

```@example chaudiere
α = 0.5

df[:, :y] = median.(df.Distribution)
df[:, :ymin] = quantile.(df.Distribution, α/2)
df[:, :ymax] = quantile.(df.Distribution, 1-α/2)

first(df, 5)
```

### Display the discharge of all configurations

```@example chaudiere
set_default_plot_size(12cm, 8cm)
plot(df, x=:Year, y=:y, color=:Configuration, Geom.line,
    ymin=:ymin, ymax=:ymax, Geom.ribbon,
    Guide.ylabel("Discharge [m³/s]"))
```


### Display the discharge of one specific configuration

```@example chaudiere
config = "MG24HQ"

df2 = filter(row -> row.Configuration ==config, df)

set_default_plot_size(12cm, 8cm)
plot(df2, x=:Year, y=:y, Geom.line, Geom.point,
    ymin=:ymin, ymax=:ymax, Geom.ribbon,
    Guide.ylabel("Discharge [m³/s]"))
```

---
## Stationary extreme value model

### Statistical model definition

```@example chaudiere
model1 = PseudoMaximaModel(pdata, prior=[Flat(), Flat(), LocationScale(-.5, 1, Beta(6,9))])
```

### Sampling of the posterior distribution by MCMC

```@example chaudiere
fm1 = fitbayes(model1)
```

### Display of the annual maximum estimates

```@example chaudiere
# Compute the estimates
ŷ = vec(mean( fm1.maxima[:, :, 1].value, dims=1))

# Compute the 95% credible interval
ymin = [quantile(vec(fm1.maxima[:, j, 1].value), .025) for j in 1:length(ŷ)]
ymax= [quantile(vec(fm1.maxima[:, j, 1].value), .975) for j in 1:length(ŷ)]

df = DataFrame(Year = years, Discharge=ŷ, ymin = ymin, ymax=ymax)

# Plot the ponctual estimates and the intervals
set_default_plot_size(12cm, 8cm)
plot(df, x=:Year, y=:Discharge, Geom.line, Geom.point,
    ymin=ymin, ymax=ymax, Geom.ribbon
)
```

### Fit of the GEV law adjusted to the estimates of the annual maxima

#### General fit with residuals

```@example chaudiere
set_default_plot_size(21cm ,16cm)
diagnosticplots(fm1, step = 20)
```

#### Fit for a single MCMC iteration

```@example chaudiere
set_default_plot_size(21cm ,16cm)
ErrorsInVariablesExtremes.diagnosticplots_iter(fm1, 200)
```

### Quantile estimation

```@example chaudiere
# Return period
T = 20

r = vec(quantile(fm1, 1-1/T))

set_default_plot_size(12cm, 8cm)
plot(x=r, Geom.density,
    Guide.xlabel(string(T,"-year return level [m³/s]")), Guide.ylabel("Density"))
```

---
## Non stationary extreme value model

The location parameter of the GEV law is a function of the GHG concentration in the atmosphere.

### Loading the covariable

```@example chaudiere
co2data = CSV.read("../../../test/data/RCPdata.csv", DataFrame);
filter!(row->row.Year in pdata[1].year, co2data)

co2data[:, :RCP85std] = (co2data.RCP85 .- mean(co2data.RCP85)) ./ std(co2data.RCP85)

co2data[:, :t] = (co2data.Year .- mean(co2data.Year)) ./ std(co2data.Year)

locationcov= Extremes.buildVariables(co2data, [:RCP85std])
```

### Statistical model definition

```@example chaudiere
model2 = PseudoMaximaModel(pdata, locationcov = locationcov, 
    prior=[Flat(), Flat(), Flat(), LocationScale(-.5, 1, Beta(6,9))])
```

### Sampling of the posterior distribution by MCMC

```@example chaudiere
fm2 = fitbayes(model2)
```

### Display of the annual maximum estimates

```@example chaudiere
# Compute the estimates
ŷ = vec(mean( fm2.maxima[:, :, 1].value, dims=1))

# Compute the 95% credible interval
ymin = [quantile(vec(fm2.maxima[:, j, 1].value), .025) for j in 1:length(ŷ)]
ymax= [quantile(vec(fm2.maxima[:, j, 1].value), .975) for j in 1:length(ŷ)]

df = DataFrame(Year = years, Discharge=ŷ, ymin = ymin, ymax=ymax)

# Plot the ponctual estimates and the intervals
set_default_plot_size(12cm, 8cm)
plot(df, x=:Year, y=:Discharge, Geom.line, Geom.point,
    ymin=ymin, ymax=ymax, Geom.ribbon
)
```

### Fit of the GEV law adjusted to the estimates of the annual maxima

#### General fit with residuals

```@example chaudiere
set_default_plot_size(21cm ,16cm)
diagnosticplots(fm2, step = 20)
```

#### Fit for a single MCMC iteration

```@example chaudiere
set_default_plot_size(21cm ,16cm)
ErrorsInVariablesExtremes.diagnosticplots_iter(fm2, 200)
```

### Quantile estimation

```@example chaudiere
# Return period
T = 20

r = quantile(fm2, 1-1/T)

df = DataFrame(Iteration = 1:size(r,1))

colnames = string.(pdata[1].year)

for j in 1:size(r,2)
   
    df[:, colnames[j]] = r[:,j]
    
end

longdf = stack(df, Not(:Iteration))
rename!(longdf, :variable => :Year )
longdf[!,:Year] = parse.(Int,longdf.Year)


df_quantile = combine(groupby(longdf, :Year),
                :value => (x -> quantile(x, .025)) => :ymin,
                :value => mean => :ȳ,
                :value => (x -> quantile(x, .975)) => :ymax)

set_default_plot_size(12cm, 8cm)
plot(df_quantile, x=:Year, y=:ȳ, Geom.line,
    ymin=:ymin, ymax=:ymax, Geom.ribbon,
    Guide.ylabel(string(T,"-year return level [m³/s]")))
```

---
## Model selection

```@example chaudiere
res = dic.([fm1, fm2])
argmin(res)
```