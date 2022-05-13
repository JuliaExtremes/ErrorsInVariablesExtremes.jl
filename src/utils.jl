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

"""
    quantile2gaussian(x::Vector{<:Real}, p::Vector{<:Real})

Determine the Normal distribution with `Φ(x[1]) = p[1]` and `Φ(x[2]) = p[2]`.
"""
function quantile2gaussian(x::Vector{<:Real}, p::Vector{<:Real})
   
    @assert length(x) == length(p) "The size of the quantiles vector must match the size of the probability vector."
    @assert length(unique(p)) == 2 "Two unique quantiles are required to determine the Normal distribution."
    
    σ = (x[2] - x[1]) / (quantile(Normal(), p[2]) - quantile(Normal(), p[1]))
    
    μ = ( x[1]*quantile(Normal(),p[2]) - x[2]*quantile(Normal(),p[1]) ) / ( quantile(Normal(), p[2]) - quantile(Normal(), p[1]) )
    
    return Normal(μ, σ)
    
end

"""
    update_stepsize(δ::Real, accrate::Real)

Update of the random walk step size for the Metropolis-Hastings algorithm

#### Details

#TODO : See Rosenthal ...
"""
function update_stepsize(δ::Real, accrate::Real)
    Δδ = 0.01 * (2 * (accrate > 0.44) - 1)
    return δ * exp(Δδ)
end