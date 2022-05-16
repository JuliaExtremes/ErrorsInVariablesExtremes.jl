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