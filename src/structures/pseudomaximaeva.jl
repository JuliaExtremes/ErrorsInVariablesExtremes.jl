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


# function Base.show(io::IO, obj::PseudoMaximaEVA)

#     showPseudoMaximaEVA(io, obj)

# end

# function showPseudoMaximaEVA(io::IO, obj::PseudoMaximaEVA; prefix::String = "")

#     show(io::IO, obj.pseudodata, prefix ="  ")
#     println(io::IO, "")
#     println(io::IO, "parameters:")
#     Extremes.showfittedEVA(io::IO, obj.parameters; prefix="  ")
#     println(io::IO, "")
#     println(io::IO, "maxima:")
#     Extremes.showChain(io::IO, obj.maxima, prefix="  ")

# end

"""
    dic(fm::PseudoMaximaEVA)

Compute the Deviance Information Criterion (DIC) described by Gelman et al. (2013) for the PseudoMaximaEVA model `fm`.

#### Details

Reference:
Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A. & Rubin, D.B. (2013). *Bayesian Data Analysis (3rd ed.)*. Chapman and Hall/CRC. https://doi.org/10.1201/b16018
"""
# function dic(fm::PseudoMaximaEVA)
   
#     Davg = mean(loglike(fm))
    
#     θ̂ = vec(mean(fm.parameters.sim.value, dims=1))
#     ŷ = vec(mean(fm.maxima.value, dims=1))
    
#     D = loglike(fm, ŷ, θ̂) 
    
#     return 2*Davg - D
        
# end


"""
    function loglike(fm::PseudoMaximaEVA)

Compute the loglikelihood of the ErrorsInVariables extreme value model `fm` for all MCMC iterations.
"""
# function loglike(fm::PseudoMaximaEVA)
    
#     # Number of MCMC simulations
#     nsim = size(fm.parameters.sim.value, 1)
    
#     # Preallocating the vector
#     ll = Vector{Float64}(undef, nsim)
    
#     for k in 1:nsim
    
#         y = vec(fm.maxima.value[k, :, 1])
#         θ̂ = vec(fm.parameters.sim.value[k,:,1])
        
#         ll[k] = loglike(fm, y, θ̂)
    
#     end
    
#     return ll
    
# end

"""
    function loglike(fm::PseudoMaximaEVA, y::Vector{<:Real}, θ::AbstractVector{<:Real})

Compute the loglikelihood of the ErrorsInVariables extreme value model `fm` given the maxima `y` and the GEV parameters `θ`.
"""
# function loglike(fm::PseudoMaximaEVA, y::Vector{<:Real}, θ::AbstractVector{<:Real})
    
#     # Reconstruct the BlockMaxima structure with the data y
#     model = BlockMaxima(Variable("y", y),
#         locationcov = fm.parameters.model.location.covariate,
#         logscalecov = fm.parameters.model.logscale.covariate,
#         shapecov = fm.parameters.model.shape.covariate)
    
#     # Evaluate the loglikelihood knowing the maxima
#     ℓ₁ = Extremes.loglike(model, θ)
    
#     # Evaluate the loglikelihood of the maxima
#     ℓ₂ = sum(logpdf(fm.pseudodata, y))
    
#     return ℓ₁ + ℓ₂
    
# end