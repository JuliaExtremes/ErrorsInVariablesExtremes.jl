struct PseudoMaximaModel
    datadistribution::Vector{Pseudodata}
    location::paramfun
    logscale::paramfun
    shape::paramfun
    prior::Vector{<:ContinuousUnivariateDistribution}
end

"""
    PseudoMaximaModel(data::Vector{Pseudodata};
        locationcov::Vector{Variable} = Vector{Variable}(),
        logscalecov::Vector{Variable} = Vector{Variable}(),
        shapecov::Vector{Variable} = Vector{Variable}())

Creates a PseudoMaximaModel structure.
"""
function PseudoMaximaModel(data::Vector{Pseudodata};
    locationcov::Vector{<:DataItem} = Vector{Variable}(),
    logscalecov::Vector{<:DataItem} = Vector{Variable}(),
    shapecov::Vector{<:DataItem} = Vector{Variable}(),
    prior::Vector{<:ContinuousUnivariateDistribution} = ContinuousUnivariateDistribution[])
    
    n = length(data[1].value)
    y = Vector{Float64}(undef, n)
    
    emptymodel = BlockMaxima(Variable("y", y),
        locationcov = locationcov,
        logscalecov = logscalecov,
        shapecov = shapecov)
    
    p = Extremes.nparameter(emptymodel)
    
    if isempty(prior)
        prior = Vector{Distribution}(undef, p)
        prior[1:p] .= Flat()
    else
        validateprior(prior, p)
    end
    
    emptymodel.location
    
    return PseudoMaximaModel(data, emptymodel.location, emptymodel.logscale, emptymodel.shape, prior)

end

