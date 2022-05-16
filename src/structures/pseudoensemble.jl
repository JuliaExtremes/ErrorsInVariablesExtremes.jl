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


function gevfitbayes(pdata::Pseudoensemble;
    locationcov::Vector{<:DataItem} = Vector{Variable}(),
    logscalecov::Vector{<:DataItem} = Vector{Variable}(),
    shapecov::Vector{<:DataItem} = Vector{Variable}(),
    prior::Vector{ContinuousUnivariateDistribution} = Vector{ContinuousUnivariateDistribution}(),
    δₒ::Real=0,
    δ::Vector{<:Real}=Float64[],
    warmup::Int=10000,
    thin::Int=10,
    niter::Int=20000,
    adapt::Symbol=:warmup)

    n = length(pdata.value[1].value)

    y₀ = Variable("pdata", ensemblemean(pdata))

    data_layer = BlockMaxima(y₀,
        locationcov = locationcov, 
        logscalecov = logscalecov,
        shapecov = shapecov)
        
    fm = Extremes.fit(data_layer)

    if isempty(δ) | (δₒ == 0) 
        pvar = parametervar(fm)
        
        if isempty(δ)
            δ = sqrt.(diag(pvar))
        end
        if δₒ == 0
            δₒ = sqrt(pvar[1])
        end
    end

    Y = Array{Float64}(undef, n, niter)
    Y[:,1] = y₀.value

    params = Array{Float64}(undef, Extremes.nparameter(data_layer), niter)
    params[:,1] = fm.θ̂


    acc_y = falses(n, niter)
    acc = falses(Extremes.nparameter(data_layer), niter)

    if isempty(prior)
        prior = Vector{ContinuousUnivariateDistribution}(undef, Extremes.nparameter(data_layer))

        parindex = Extremes.paramindex(data_layer)
        prior[parindex[:μ]] .= Flat()
        prior[parindex[:ϕ]] .= Flat()
        prior[parindex[:ξ]] .= LocationScale(-.5, 1, Beta(6,9))
    end


    @showprogress for iter in 2:niter
        

        if Extremes.nparameter(data_layer) == 3
            pd = repeat(Extremes.getdistribution(data_layer, params[:, iter-1] ), n)
        else
            pd = Extremes.getdistribution(data_layer, params[:, iter-1])
        end
        
        Y[:, iter] = Y[:, iter-1]
        
        for j in 1:n
            
            ỹ = Y[j, iter] + δₒ*randn()
            
            ll = logpdf(pdata, ỹ, j) - 
                logpdf(pdata, Y[j,iter], j) +
                logpdf(pd[j], ỹ) -
                logpdf(pd[j], Y[j, iter])
            
            if ll > log(rand())
                Y[j, iter] = ỹ
                acc_y[j, iter] = true
            end
            
        end
    
    
        data_layer = BlockMaxima(Variable("y", Y[:, iter]),
            locationcov = locationcov, 
            logscalecov = logscalecov,
            shapecov = shapecov)
        
        params[:, iter] = params[:, iter-1]
        
        for k in 1:Extremes.nparameter(data_layer)
            
            θ̃ = params[:,iter]
            θ̃[k] = params[k, iter] + δ[k]*randn()
            
            ll = Extremes.loglike(data_layer, θ̃) + logpdf(prior[k], θ̃[k]) -
                Extremes.loglike(data_layer, params[:,iter]) - logpdf(prior[k], params[k,iter])
            
            if ll > log(rand())
                params[k, iter] = θ̃[k]
                acc[k, iter] = true
            end
            
        end
                
        # Updating the stepsize
        if iter % 50 == 0
            if !(adapt == :none)
                if (iter <= warmup) | (adapt==:all)
                    accrate = mean(acc_y[:, iter-50+1:iter])
                    δₒ = update_stepsize(δₒ, accrate)
                    
                    accrate = vec(mean(acc[:, iter-50+1:iter], dims=2))
                    δ = update_stepsize.(δ, accrate)
                end
            end
        end
        
    end


    println("GEV parameters acceptance rate")
    println("    ", vec(mean(acc[:, warmup:niter], dims=2)))
    println("")

    println("Pseuso-maxima mean acceptance rate")
    println("    ", mean(acc_y[:, warmup:niter]))
    println("")

    #Extracting output


    paramnames = String[]

    if length(parindex[:μ]) >1
        append!(paramnames, ["μ[$j]" for j in 0:length(parindex[:μ])-1])
    else
        push!(paramnames, "μ")
    end

    if length(parindex[:ϕ]) >1
        append!(paramnames, ["ϕ[$j]" for j in 0:length(parindex[:ϕ])-1])
    else
        push!(paramnames, "ϕ")
    end

    if length(parindex[:ξ]) >1
        append!(paramnames, ["ξ[$j]" for j in 0:length(parindex[:ξ])-1])
    else
        push!(paramnames, "ξ")
    end

    C = Mamba.Chains(collect(params'), names=paramnames)

    ŷ = vec(mean(Y[:, warmup:thin:niter], dims=2))

    data_layer = BlockMaxima(Variable("ŷ", ŷ),
            locationcov = locationcov, 
            logscalecov = logscalecov,
            shapecov = shapecov)

    fm = BayesianEVA(data_layer, C[warmup+1:thin:niter, :, :])

    M_names = ["Y[$j]" for j = 1:n]

    CM = Mamba.Chains(collect(Y'), names=M_names)

    return PseudoMaximaEVA(pdata, fm, CM[warmup:thin:niter, :,:])

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

"""
    showpseudoensemble(io::IO, obj::Pseudoensemble; prefix::String = "")

Displays a Pseudoensemble with the prefix `prefix` before every line.
"""
function showpseudoensemble(io::IO, obj::Pseudoensemble; prefix::String = "")

    println(io, prefix, "Pseudoensemble:")
    println(io, prefix, "  name: ", obj.name)
    println(io, prefix, "  value: ", typeof(obj.value), "[", length(obj.value), "]")
end

"""
    Base.show(io::IO, obj::Pseudoensemble)

Override of the show function for the objects of type Pseudoensemble.
"""
function Base.show(io::IO, obj::Pseudoensemble; prefix::String = "")

    showpseudoensemble(io, obj, prefix=prefix)

end