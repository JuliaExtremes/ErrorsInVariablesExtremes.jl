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


# Methods for Pseudodata

"""
    convert(::Type{DataFrame}, pdata::Vector{Pseudodata})

Convert the vector of Pseudodata type to a DataFrame
"""
function convert(::Type{DataFrame}, pdata::Vector{Pseudodata})
    
    S = length(pdata)
    n = length(pdata[1].value)
    
    df = DataFrame(Year = Int64[], Configuration = String[], Distribution = ContinuousUnivariateDistribution[])
    
    for i in 1:S
        for j in 1:n
            push!(df, [pdata[i].year[j], pdata[i].name, pdata[i].value[j]])
        end
    end
    
    return df
        
end

"""
    ensemblemean(pdata::Vector{Pseudodata})

Compute the ensemble mean for each year.
"""
function ensemblemean(pdata::Vector{Pseudodata})
    
    S = length(pdata)

    pd = [pdata[i].value for i in 1:S]
    
    m = vec(mean(mean.(reduce(hcat,pd)), dims=2))
    
end


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
    logpdf(pensemble::Pseudoensemble, y::Vector{<:Real})

Compute the logpdf of each of the potential data `y` according to the distributions in `pensemble`.

#### Details

Independance is assumed between the members of the `pseudoensemble.value`, i.e. the sum of the logpdf of each member is taken.
"""
function logpdf(datadistribution::Vector{Pseudodata}, y::Vector{<:Real})
    
    ll = zeros(length(y))
    
    for pdata in datadistribution
     
        ll += logpdf(pdata, y)
        
    end
    
    return ll
    
end

function logpdf(datadistribution::Vector{Pseudodata}, y::Real, j::Int)
    
    ll = 0.0
    
    for pdata in datadistribution
     
        ll += logpdf(pdata.value[j], y)
        
    end
    
    return ll
    
end

