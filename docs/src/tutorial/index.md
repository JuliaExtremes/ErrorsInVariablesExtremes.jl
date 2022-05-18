# Tutorial

This tutorial shows the functionalities of *ErrorsInVariablesExtremes.jl*. They are illustrated by producing an Extreme Value Analysis of unobserved riverflows at ungauged section of the Chaudière River.

Before executing this tutorial, make sure to have installed the following packages:
- *ErrorsInVariablesExtremes.jl* (of course)
- *Extremes.jl* (for more Extreme Value Analysis tools)
- *DataFrames.jl* (for using the DataFrame type)
- *Distributions.jl* (for using probability distribution objects)
- *Gadfly.jl* (for plotting)
- *CSV.jl* (to load CSV files)

and import them using the following command:
 ```@repl
using ErrorsInVariablesExtremes, Extremes, Dates, DataFrames, Distributions, Gadfly, CSV
```

