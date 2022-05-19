"""
    PseudoMaximaEVA(pseudodata::Pseudoensemble, model::EVA, sim::Mamba.Chains)

Construct a PseudoMaximaEVA type.

#### Details

Encapsulates the probability distributions for each of the unobserved data for an ensemble of pseudodata, the extreme value model and a sample from the parameters posterior distribution.

# TODO: Verify if all member has the same year vector.

"""
struct PseudoMaximaEVA
    model::PseudoMaximaModel
    maxima::Mamba.Chains
    parameters::Mamba.Chains
end


function Base.show(io::IO, obj::PseudoMaximaEVA)

    showPseudoMaximaEVA(io, obj)

end

function showPseudoMaximaEVA(io::IO, obj::PseudoMaximaEVA; prefix::String = "")

    println(io, prefix, "PseudoMaximaEVA")
    println(io, prefix, "   model: ", typeof(obj.model))
    println(io::IO, "   maxima:")
    Extremes.showChain(io::IO, obj.maxima, prefix="\t\t")
    println(io::IO, "   parameters:")
    Extremes.showChain(io::IO, obj.parameters; prefix="\t\t")
end


# Methods for PseudoMaximaEVA

"""
    dic(fm::PseudoMaximaEVA)

Compute the Deviance Information Criterion (DIC) described by Gelman et al. (2013) for the PseudoMaximaEVA model `fm`.

#### Details

Reference:
Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A. & Rubin, D.B. (2013). *Bayesian Data Analysis (3rd ed.)*. Chapman and Hall/CRC. https://doi.org/10.1201/b16018
"""
function dic(fm::PseudoMaximaEVA)
   
    Davg = mean(logpdf(fm))
    
    ŷ, θ̂ = ErrorsInVariablesExtremes.findposteriormode(fm)
    
    D = logpdf(fm.model, ŷ, θ̂)
    
    return 2*Davg - D
        
end


"""
    findposteriormode(fm::PseudoMaximaEVA)

Find the point (ŷ, θ̂) among the MCMC iterations of `fm` that maximized the log density.
"""
function findposteriormode(fm::PseudoMaximaEVA)
    
    nsim = size(fm.maxima.value,1)
    
    ll = Vector{Float64}(undef, nsim)
    
    for k in 1:nsim
        y = vec(fm.maxima.value[k,:,1])
        θ = vec(fm.parameters.value[k,:,1])
        
        ll[k] = logpdf(fm.model, y, θ)
    end
    
    ind = argmax(ll)
    
    ŷ = vec(fm.maxima.value[ind,:,1])
    θ̂ = vec(fm.parameters.value[ind,:,1])
    
    return ŷ, θ̂
    
end


"""
    function logpdf(fm::PseudoMaximaEVA)

Compute the log density of `fm` for all MCMC iterations.
"""
function logpdf(fm::PseudoMaximaEVA)
    
      # Number of MCMC simulations
    nsim = size(fm.maxima.value, 1)
    
    # Preallocating the vector
    ll = Vector{Float64}(undef, nsim)

    for k in 1:nsim
    
        y = vec(fm.maxima.value[k, :, 1])
        θ̂ = vec(fm.parameters.value[k,:,1])
        
        ll[k] = logpdf(fm.model, y, θ̂)
    
    end
    
    return ll
    
end
