struct PseudoMaximaModel
    datadistribution::Vector{Pseudodata}
    location::paramfun
    logscale::paramfun
    shape::paramfun
    prior::Vector{<:ContinuousUnivariateDistribution}
end

"""
    PseudoMaximaModel(data::Vector{Pseudodata};
        locationcov::Vector{Variable} = Vector{Variable}(),
        logscalecov::Vector{Variable} = Vector{Variable}(),
        shapecov::Vector{Variable} = Vector{Variable}())

Creates a PseudoMaximaModel structure.
"""
function PseudoMaximaModel(data::Vector{Pseudodata};
    locationcov::Vector{<:DataItem} = Vector{Variable}(),
    logscalecov::Vector{<:DataItem} = Vector{Variable}(),
    shapecov::Vector{<:DataItem} = Vector{Variable}(),
    prior::Vector{<:ContinuousUnivariateDistribution} = ContinuousUnivariateDistribution[])
    
    n = length(data[1].value)
    y = Vector{Float64}(undef, n)
    
    emptymodel = BlockMaxima(Variable("y", y),
        locationcov = locationcov,
        logscalecov = logscalecov,
        shapecov = shapecov)
    
    p = Extremes.nparameter(emptymodel)
    
    if isempty(prior)
        prior = Vector{Distribution}(undef, p)
        prior[1:p] .= Flat()
    else
        validateprior(prior, p)
    end
    
    emptymodel.location
    
    return PseudoMaximaModel(data, emptymodel.location, emptymodel.logscale, emptymodel.shape, prior)

end


function fitbayes(model::PseudoMaximaModel,
    δₒ::Real=0,
    δ::Vector{<:Real}=Float64[],
    warmup::Int=10000,
    thin::Int=10,
    niter::Int=20000,
    adapt::Symbol=:warmup)

    pdata = model.datadistribution
    
    n = length(pdata[1].value)

    y₀ = Variable("pdata", ensemblemean(pdata))

    data_layer = BlockMaxima(y₀,
        model.location,
        model.logscale,
        model.shape)
        
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
            model.location,
            model.logscale,
            model.shape)
        
        params[:, iter] = params[:, iter-1]
        
        for k in 1:Extremes.nparameter(data_layer)
            
            θ̃ = params[:,iter]
            θ̃[k] = params[k, iter] + δ[k]*randn()
            
            ll = Extremes.loglike(data_layer, θ̃) + logpdf(model.prior[k], θ̃[k]) -
                Extremes.loglike(data_layer, params[:,iter]) - logpdf(model.prior[k], params[k,iter])
            
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

    parindex = Extremes.paramindex(data_layer)

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
    
    maxima_chain = Mamba.Chains(collect(Y'), names=["Y[$j]" for j = 1:n])

    parameters_chain = Mamba.Chains(collect(params'), names=paramnames)
    
    return PseudoMaximaEVA(model, maxima_chain, parameters_chain)

end