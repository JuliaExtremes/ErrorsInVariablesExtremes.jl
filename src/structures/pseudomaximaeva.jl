"""
    PseudoMaximaEVA(pseudodata::Pseudoensemble, model::EVA, sim::Mamba.Chains)

Construct a PseudoMaximaEVA type.

#### Details

Encapsulates the probability distributions for each of the unobserved data for an ensemble of pseudodata, the extreme value model and a sample from the parameters posterior distribution.

# TODO: Verify if all member has the same year vector.

"""
struct PseudoMaximaEVA
    pseudodata::Pseudoensemble
    model::EVA
    sim::Mamba.Chains
end


function Base.show(io::IO, obj::PseudoMaximaEVA)

    showPseudoMaximaEVA(io, obj)

end

function showPseudoMaximaEVA(io::IO, obj::PseudoMaximaEVA; prefix::String = "")

    show(io::IO, obj.pseudodata)
    println(io::IO, "")
    show(io::IO, obj.model)
    println(io::IO, "")
    Extremes.showChain(io::IO, obj.sim)

end