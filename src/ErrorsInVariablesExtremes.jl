module ErrorsInVariablesExtremes

using CSV, DataFrames, Dates, NetCDF
using Distributions, Extremes, LinearAlgebra, Mamba, Random, Statistics
using ProgressMeter

import Base.convert
import Distributions.pdf, Distributions.logpdf
import Extremes.gevfitbayes, Extremes.loglike

include("structures.jl")
include("utils.jl")

export 

    # Pseudodata type
    Pseudodata, Pseudoensemble, PseudoMaximaEVA,

    # Other functions
    convert,
    ensemblemean,
    gevfitbayes,
    load_discharge_distribution,
    logpdf,
    loglike,
    pdf,
    get_DIC

end # module
