module ErrorsInVariablesExtremes

using CSV, DataFrames, Dates, NetCDF
using Distributions, Extremes, LinearAlgebra, Mamba, Random, Statistics
using Gadfly
using ProgressMeter

using Test

import Base.convert
import Distributions.pdf, Distributions.logpdf

include("utils.jl")
include("structures.jl")
include("parameterestimation.jl")

export 

    # Pseudodata type
    Pseudodata,

    # Pseudoensemble type
    Pseudoensemble,

    # Other functions
    convert,
    ensemblemean,
    load_discharge_distribution,
    logpdf,
    pdf

end # module
