# Utils

"""
    load_discharge_distribution(filename::String)

Load the discharge distribution from the NetCDF file `filename`.

#### Details

For each configuration of the hydrological model, the function loads the 21 quantiles and transforms them into a Log-normal distribution. The data are encapsulated in the `Pseudo set` structure.
"""

function load_discharge_distribution(filename::String)

    
    att = ncgetatt(filename, "time", "units")
    s = split(att)
    ini = Date(s[3])
    
    t = year.(ini .+ Day.(ncread(filename, "time")))
    
    
    
    id = ncread(filename, "scenario_id")[62:67,:]
    
    S = size(id, 2)
    
    config = String[]
    
    for i in 1:S
       push!(config, String(id[:,i])) 
    end
    
    data = ncread(filename, "Dis")
    
    q2 = quantile(Normal(), 0.75)
    q1 = quantile(Normal(), 0.25)
    x2 = log.(data[16, :, :])
    x1 = log.(data[6, :, :])

    ζ = (x2 .- x1) ./ (q2 - q1)
    η = (x1 .* q2 .- x2 .* q1) ./ (q2 - q1)
       
    
    pdata = Pseudodata[]
    for i in 1:S
        push!(pdata, Pseudodata(config[i], t, LogNormal.(η[:,i], ζ[:,i])))
    end
    
    pensemble = Pseudoensemble("Ensemble", pdata)
    
        
    return pensemble

end