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

"""
    showpseudodata(io::IO, obj::Pseudodata; prefix::String = "")

Displays a Pseudodata with the prefix `prefix` before every line.
"""
function showpseudodata(io::IO, obj::Pseudodata; prefix::String = "")

    println(io, prefix, "Pseudodata:")
    println(io, prefix, "  name: ", obj.name)
    println(io, prefix, "  year: ", typeof(obj.year), "[", length(obj.year), "]")
    println(io, prefix, "  value:", typeof(obj.value), "[", length(obj.value), "]")

end

"""
    Base.show(io::IO, obj::Pseudodata)

Override of the show function for the objects of type Pseudodata.
"""
function Base.show(io::IO, obj::Pseudodata)

    showpseudodata(io, obj)

end