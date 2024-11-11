# Utils

# """
#     load_discharge_distribution(filename::String)

# Load the discharge distribution from the NetCDF file `filename`.

# #### Details

# For each configuration of the hydrological model, the function loads the 21 quantiles and transforms them into a Log-normal distribution. The data are encapsulated in the `Pseudo set` structure.
# """
# function load_discharge_distribution(filename::String)

#     att = ncgetatt(filename, "time", "units")
#     s = split(att)
#     ini = Date(s[3])
    
#     t = year.(ini .+ Day.(ncread(filename, "time")))
#     n = length(t)
    
#     id = ncread(filename, "scenario_id")[62:67,:]
    
#     S = size(id, 2)
    
#     config = String[]
    
#     for i in 1:S
#        push!(config, String(id[:,i])) 
#     end
    
#     data = ncread(filename, "Dis")

#     x₁ = data[11, :, :] #Assuming this is the index of the median
#     x₂ = data[16, :, :] #Assuming this is the index of the 75th quantile
    
#     # Convert the quantile to Log-normal distribution
#     μ = log.(x₁)
#     σ = (log.(x₂) .- μ)/ quantile(Normal(), .75)
#     pd = LogNormal.(μ, σ)
    
#      pdata = Pseudodata[]
#     for i in 1:S
#         push!(pdata, Pseudodata(config[i], t, pd[:,i]))
#     end

#     return pdata

# end

"""
    load_discharge_distribution(filename::String)

Load the discharge distribution from the NetCDF file `filename`.

#### Details

For each configuration of the hydrological model, the function loads the 21 quantiles and transforms them into a Log-normal distribution. The data are encapsulated in the `Pseudo set` structure.

Some datasets available:
 - `SLSO00003`: Estimated discharges at the outlet of the Chaudière watershed (QC).
 - `SLSO00025`: Estimated discharges of the Chaudière River at Saint-Lambert-de-Lauzon (QC).

#### Examples
```julia-repl
julia> ErrorsInVariablesExtremes.load_discharge_distribution("SLSO00003")
```

"""
function load_discharge_distribution(name::String)

    filename = joinpath(dirname(@__FILE__), "..", "data", string(name, ".nc"))

    if !isfile(filename)

        error("There is no dataset with the name '$name'")

    else

        att = ncgetatt(filename, "time", "units")
        s = split(att)
        ini = Date(s[3])
        
        t = year.(ini .+ Day.(ncread(filename, "time")))
        n = length(t)
        
        id = ncread(filename, "scenario_id")[62:67,:]
        
        S = size(id, 2)
        
        config = String[]
        
        for i in 1:S
        push!(config, String(id[:,i])) 
        end
        
        data = ncread(filename, "Dis")

        x₁ = data[11, :, :] #Assuming this is the index of the median
        x₂ = data[16, :, :] #Assuming this is the index of the 75th quantile
        
        # Convert the quantile to Log-normal distribution
        μ = log.(x₁)
        σ = (log.(x₂) .- μ)/ quantile(Normal(), .75)
        pd = LogNormal.(μ, σ)
        
        pdata = Pseudodata[]
        for i in 1:S
            push!(pdata, Pseudodata(config[i], t, pd[:,i]))
        end

    end

    return pdata

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

"""
    validateprior(p::Int, prior::Vector{<:Distribution})

Validate that the `prior` distribution vector dimension is of dimension `p`.
"""
function validateprior(prior::Vector{<:Distribution}, p::Int)

    if length(prior) != p
        error("The prior dimension should match the parameter dimension.")
    end
end
