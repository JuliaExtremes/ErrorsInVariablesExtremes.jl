module ErrorsInVariablesExtremes

using CSV, DataFrames, Dates, NetCDF
using Distributions, Extremes, LinearAlgebra, Mamba, Random, Statistics
using ProgressMeter

import Base.convert
import Distributions.pdf, Distributions.logpdf
import Extremes.gevfitbayes

include("structures.jl")
include("utils.jl")
include("parameterestimation.jl")

export 

    # Pseudodata type
    Pseudodata, Pseudoensemble, PseudoMaximaEVA,

    # Other functions
    convert,
    ensemblemean,
    gevfitbayes,
    load_discharge_distribution,
    logpdf,
    pdf,
    get_DIC

end # module
