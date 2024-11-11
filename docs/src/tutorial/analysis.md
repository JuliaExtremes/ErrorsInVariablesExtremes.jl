# Bayesian analysis of the estimated discharges

```@setup chaudiere2
using ErrorsInVariablesExtremes, Extremes, Dates, DataFrames, Distributions, Gadfly, Cairo, Fontconfig
```

## Load the data

Loading the pseudo-observations into an object of type [`Pseudoensemble`](@ref):
```@example chaudiere2
filename = "SLSO00003"
pdata = ErrorsInVariablesExtremes.load_discharge_distribution(filename)
```
Converting in DataFrame:
```@example chaudiere2
df = convert(DataFrame, pdata)
first(df,5)
```

## Combination of the hydrological configurations for the log-discharges

See Eq. (11) and (12) of the paper.

```@example chaudiere2
df[:,:η] = first.(params.(df.Distribution))
df[:,:ζ] = last.(params.(df.Distribution))

m = Float64[]
s = Float64[]

for iYear in unique(df.Year)

    df3 = filter(row -> row.Year==iYear, df)
    mᵢ = sum(df3.η ./ df3.ζ.^2) / sum( 1 ./df3.ζ.^2)
    sᵢ = sqrt(1/sum( 1 ./df3.ζ.^2))
        
    push!(m, mᵢ)
    push!(s, sᵢ)
        
end

pdata[1].year

datadist = Pseudodata("Combined",pdata[1].year ,Normal.(m, s))
```

## Stationary model

### Specifying the stationary extreme value Bayesian model

The stationary model used the combined estimates of discharges and a semi-informative prior.


```@example chaudiere2
model1 = PseudoMaximaModel([datadist], prior=[Flat(), Flat(), LocationScale(-.5, 1, Beta(6,9))])
```

### MCMC

```@example chaudiere2
fm1 = fitbayes(model1, niter=20000, warmup=10001, thin=10, δₒ = .5)
```

### Traces of the parameters

Trace of the GEV parameters:
```@example chaudiere2
fig1 = plot(y=fm1.parameters.value[:,1], Geom.line,
    Guide.xlabel("Iteration"), Guide.ylabel("μ"));

fig2 = plot(y=exp.(fm1.parameters.value[:,2]), Geom.line,
    Guide.xlabel("Iteration"), Guide.ylabel("σ"));

fig3 = plot(y=fm1.parameters.value[:,3], Geom.line,
    Guide.xlabel("Iteration"), Guide.ylabel("ξ"));

set_default_plot_size(24cm, 12cm)
gridstack([fig1 fig2;
    fig3 Gadfly.plot()])
```

## Non-stationary model

### Specifying the non-stationary extreme value Bayesian model

The non-stationary model used the combined estimates of discharges and a semi-informative prior. The GEV location parameter is function of the year.


```@example chaudiere2
t = collect(0:59)

covariate = Variable("t", t );

model2 = PseudoMaximaModel([datadist], locationcov = [covariate], 
    prior=[Flat(), Flat(), Flat(),  LocationScale(-.5, 1, Beta(6,9))])
```

### MCMC

```@example chaudiere2
fm2 = fitbayes(model2, niter=20000, warmup=10001, thin=10, δₒ = .5);
```

### Traces of the parameters

Trace of the GEV parameters:
```@example chaudiere2
fig1 = plot(y=fm2.parameters.value[:,1], Geom.line,
    Guide.xlabel("Iteration"), Guide.ylabel("μ<sub>0</sub>"));

fig2 = plot(y=fm2.parameters.value[:,2], Geom.line,
    Guide.xlabel("Iteration"), Guide.ylabel("μ<sub>1</sub>"));

fig3 = plot(y=exp.(fm2.parameters.value[:,3]), Geom.line,
    Guide.xlabel("Iteration"), Guide.ylabel("σ"));

fig4 = plot(y=fm2.parameters.value[:,4], Geom.line,
    Guide.xlabel("Iteration"), Guide.ylabel("ξ"));

set_default_plot_size(24cm, 12cm)
gridstack([fig1 fig2;
    fig3 fig4])
```