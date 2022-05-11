# Structures

"""
    Pseudoensemble(name::String, value::Vector{Pseudodata})

Construct a Pseudoensemble type.

#### Details

Encapsulates the probability distributions for each of the unobserved data for an ensemble of pseudodata. For the moment, independance of each Pseudodata elements is assumed.

# TODO: Verify if all member has the same year vector.

"""
struct Pseudoensemble
    "Description"
    name::String
    "Distribution modelisation of pseudodata ensemble"
    value::Vector{Pseudodata}
end


# Methods for Pseudoensemble

"""
    convert(::Type{DataFrame}, pensemble::Pseudoensemble)

Convert the Pseudoensemble type to a DataFrame
"""
function convert(::Type{DataFrame}, pensemble::Pseudoensemble)
    
    S = length(pensemble.value)
    n = length(pensemble.value[1].year)
    
    df = DataFrame(Year = Int64[], Configuration = String[], Distribution = ContinuousUnivariateDistribution[])
    
    for i in 1:S
        pdata = pensemble.value[i]
        for j in 1:n
            
            push!(df, [pdata.year[j], pdata.name, pdata.value[j]])
            
        end
    end
    
    return df
        
end


"""
    ensemblemean(pdata::Pseudoensemble)

Compute the ensemble mean for each year.
"""
function ensemblemean(pdata::Pseudoensemble)
    
S = length(pdata.value)

    pd = [pdata.value[i].value for i in 1:S]
    
    m = vec(mean(mean.(reduce(hcat,pd)), dims=2))
    
end

"""
    logpdf(pensemble::Pseudoensemble, y::Vector{<:Real})

Compute the logpdf of each of the potential data `y` according to the distributions in `pensemble`.

#### Details

Independance is assumed between the members of the `pseudoensemble.value`, i.e. the sum of the logpdf of each member is taken.
"""
function logpdf(pensemble::Pseudoensemble, y::Vector{<:Real})
    
    ll = zeros(length(y))
    
    for pd in pensemble.value
     
        ll += logpdf(pd, y)
        
    end
    
    return ll
    
end

function logpdf(pensemble::Pseudoensemble, y::Real, j::Int)
    
    ll = 0.0
    
    for pd in pensemble.value
     
        ll += logpdf(pd.value[j], y)
        
    end
    
    return ll
    
end


"""
    pdf(pensemble::Pseudoensemble, y::Vector{<:Real})

Compute the pdf of each of the potential data `y` according to the distributions in `pensemble`.

#### Details

Independance is assumed between the members of the `pseudoensemble.value`, i.e. the product of the pdf of each member is taken.
"""
function pdf(pensemble::Pseudoensemble, y::Vector{<:Real})
       
    return exp.(logpdf(pensemble, y))
    
end

