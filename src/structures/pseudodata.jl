# Structures

"""
    Pseudodata(name::String, value::Vector{<:UnivariateDistribution})

Construct a Pseudodata type.

#### Details

Encapsulates the probability distributions for each of the unobserved data. The data are not directly observed but their distributions are known. The distributions can be different for each of the data.
"""
struct Pseudodata
    "Description"
    name::String
    "Year"
    year::Vector{<:Int}
    "Distribution modelisation of pseudodata"
    value::Vector{<:UnivariateDistribution}
end


# Methods for Pseudodata

"""
    logpdf(pdata::Pseudodata, y::Vector{<:Real})

Compute the log density of each of the potential data `y` according to the distributions in `pdata`.
"""
function logpdf(pdata::Pseudodata, y::Vector{<:Real})
    
    pd = pdata.value
    
    return logpdf.(pd, y)
    
end

"""
    pdf(pdata::Pseudodata, y::Vector{<:Real})

Compute the density of each of the potential data `y` according to the distributions in `pdata`.
"""
function pdf(pdata::Pseudodata, y::Vector{<:Real})
       
    pd = pdata.value
    
    return pdf.(pd, y)
    
end   