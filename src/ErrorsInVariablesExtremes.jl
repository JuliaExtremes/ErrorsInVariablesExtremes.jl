module ErrorsInVariablesExtremes

using CSV, DataFrames, Dates, NetCDF
using Distributions, Extremes, LinearAlgebra, Mamba, Random, Statistics
using Gadfly
using ProgressMeter

import Base.convert
import Distributions.pdf, Distributions.logpdf
import Extremes.fitbayes, Extremes.loglike, Extremes.paramfun
import Extremes.diagnosticplots, Extremes.histplot, Extremes.probplot, Extremes.qqplot, Extremes.returnlevelplot
import Mamba.dic

include("structures.jl")
include("utils.jl")

export 

    # Pseudodata type
    Pseudodata, PseudoMaximaEVA, PseudoMaximaModel,

    # Other functions
    convert,
    dic,
    ensemblemean,
    fitbayes,
    load_discharge_distribution,
    logpdf,
    loglike,
    pdf

    # Diagnostic plots
    diagnosticplots, 
    histplot, 
    probplot,
    qqplot,
    returnlevelplot

end # module
