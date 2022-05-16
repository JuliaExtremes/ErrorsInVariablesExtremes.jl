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

    pensemble = Pseudoensemble("Ensemble", pdata)
    
    return pensemble

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
    get_DIC(fm::PseudoMaximaEVA)

Compute the Deviance Information Criterion (DIC).

#### Details

#TODO : See Spiegelhalter et al. (2002) ...

"""
function get_DIC(fm::PseudoMaximaEVA)
    
    y = fm.maxima.value[:, :, 1]
    ȳ = vec(mean(y, dims = 1))

    v = fm.parameters.sim.value[:,:,1]
    V = Extremes.slicematrix(v, dims=2)
    V̄ = vec(mean(V, dims = 1))

    D = Extremes.getdistribution.(fm.parameters.model, V̄)
    d = Extremes.unslicematrix(D, dims=2)

    DIC₁ = sum(logpdf(fm.pseudodata, ȳ)) + sum(logpdf.(d, ȳ))
    
    DIC₂ = 0.0
    for s in 1:length(fm.pseudodata.value)
        DIC₂ += sum(mean(logpdf.(Ref(fm.pseudodata), y, Ref(s)), dims=1))
    end
    DIC₂ += sum(mean(logpdf.(Extremes.getdistribution(fm.parameters), y), dims=1))

    return -2 * DIC₂ + DIC₁
end