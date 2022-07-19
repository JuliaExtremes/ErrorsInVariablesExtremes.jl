# Tutorial

This tutorial shows the functionalities of *ErrorsInVariablesExtremes.jl*. They are illustrated by producing an Extreme Value Analysis of unobserved riverflows of a section of the Chaudière River in Québec, Canada.

Before executing this tutorial, make sure to have installed the following packages:
- *ErrorsInVariablesExtremes.jl* 
- *Extremes.jl* (for more Extreme Value Analysis tools)
- *DataFrames.jl* (for using the DataFrame type)
- *CSV.jl* (to load CSV files)
- *Distributions.jl* (for using probability distribution objects)
- *Mamba.jl* (for more MCMC tools)
- *Gadfly.jl* (for plotting)


and import them using the following command:
 ```@repl
using ErrorsInVariablesExtremes, Extremes, DataFrames, CSV, Distributions, Mamba, Gadfly
```

