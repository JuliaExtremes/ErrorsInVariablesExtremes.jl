"""
    PseudoMaximaEVA(pseudodata::Pseudoensemble, model::EVA, sim::Mamba.Chains)

Construct a PseudoMaximaEVA type.

#### Details

Encapsulates the probability distributions for each of the unobserved data for an ensemble of pseudodata, the extreme value model and a sample from the parameters posterior distribution.

# TODO: Verify if all member has the same year vector.

"""
struct PseudoMaximaEVA
    pseudodata::Pseudoensemble
    parameters::fittedEVA
    maxima::Mamba.Chains
end


function Base.show(io::IO, obj::PseudoMaximaEVA)

    showPseudoMaximaEVA(io, obj)

end

function showPseudoMaximaEVA(io::IO, obj::PseudoMaximaEVA; prefix::String = "")

    show(io::IO, obj.pseudodata, prefix ="  ")
    println(io::IO, "")
    println(io::IO, "parameters:")
    Extremes.showfittedEVA(io::IO, obj.parameters; prefix="  ")
    println(io::IO, "")
    println(io::IO, "maxima:")
    Extremes.showChain(io::IO, obj.maxima, prefix="  ")

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