"""
    PseudoMaximaEVA(pseudodata::Pseudoensemble, model::EVA, sim::Mamba.Chains)

Construct a PseudoMaximaEVA type.

#### Details

Encapsulates the probability distributions for each of the unobserved data for an ensemble of pseudodata, the extreme value model and a sample from the parameters posterior distribution.

# TODO: Verify if all member has the same year vector.

"""
struct PseudoMaximaEVA
    model::PseudoMaximaModel
    maxima::Mamba.Chains
    parameters::Mamba.Chains
end

Base.Broadcast.broadcastable(obj::PseudoMaximaEVA) = Ref(obj)

function Base.show(io::IO, obj::PseudoMaximaEVA)

    showPseudoMaximaEVA(io, obj)

end

function showPseudoMaximaEVA(io::IO, obj::PseudoMaximaEVA; prefix::String = "")

    println(io, prefix, "PseudoMaximaEVA")
    println(io, prefix, "   model: ", typeof(obj.model))
    println(io::IO, "   maxima:")
    Extremes.showChain(io::IO, obj.maxima, prefix="\t\t")
    println(io::IO, "   parameters:")
    Extremes.showChain(io::IO, obj.parameters; prefix="\t\t")
end


# Methods for PseudoMaximaEVA


"""
    convert(::Type{MaximumLikelihoodEVA}, fm::PseudoMaximaEVA, iter::Int)

Convert a single MCMC iteration of the fm model to a pseudo MaximumLikelihoodEVA for hacking Extremes.jl methods. 
"""
function convert(::Type{MaximumLikelihoodEVA}, fm::PseudoMaximaEVA, iter::Int)
    
        bm = BlockMaxima(Variable("y", fm.maxima.value[iter, :, 1]),
        fm.model.location, fm.model.logscale, fm.model.shape)

        return MaximumLikelihoodEVA(bm, fm.parameters.value[iter, :, 1])
    
end

"""
    dic(fm::PseudoMaximaEVA)

Compute the Deviance Information Criterion (DIC) described by Gelman et al. (2013) for the PseudoMaximaEVA model `fm`.

#### Details

Reference:
Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A. & Rubin, D.B. (2013). *Bayesian Data Analysis (3rd ed.)*. Chapman and Hall/CRC. https://doi.org/10.1201/b16018
"""
function dic(fm::PseudoMaximaEVA)
   
    Davg = mean(logpdf(fm))
    
    ŷ = vec(mean(fm.maxima.value[:,:,1], dims=1))
    θ̂ = vec(mean(fm.parameters.value[:,:,1], dims=1))
    
    D = logpdf(fm.model, ŷ, θ̂)
    
    return 2*Davg - D
        
end


"""
    getdistribution(fm::PseudoMaximaEVA, iter::Int)

Return the underlying extreme value distribution of the MCMC iteration `iter`.
"""
function getdistribution(fm::PseudoMaximaEVA, iter::Int)
    
    @assert iter > 0 "The specified MCMC iteration shoud be greater than 0."
    
    θᵢ = fm.parameters.value[iter,:,1]

    return ErrorsInVariablesExtremes.getdistribution(fm.model, θᵢ)

end

"""
    getdistribution(fm::PseudoMaximaEVA)

Return the underlying extreme value distribution for each MCMC iterations.
"""
function getdistribution(fm::PseudoMaximaEVA)
    
    nsim = size(fm.maxima.value,1)
    
    D = getdistribution.(fm, 1:nsim)

    return Extremes.unslicematrix(D, dims=2)

end


"""
    function logpdf(fm::PseudoMaximaEVA)

Compute the log density of `fm` for all MCMC iterations.
"""
function logpdf(fm::PseudoMaximaEVA)
    
      # Number of MCMC simulations
    nsim = size(fm.maxima.value, 1)
    
    # Preallocating the vector
    ll = Vector{Float64}(undef, nsim)

    for k in 1:nsim
    
        y = vec(fm.maxima.value[k, :, 1])
        θ̂ = vec(fm.parameters.value[k,:,1])
        
        ll[k] = logpdf(fm.model, y, θ̂)
    
    end
    
    return ll
    
end

"""
    thin(fm::PseudoMaximaEVA, step::Int)

Discard all but every `step`th MCMC iteration of `fm.maxima` and `fm.parameters`.
"""
function thin(fm::PseudoMaximaEVA, step::Int)
   
    @assert step > 0 "The thinning step shoud be a strictly positive integer."
    
    res = PseudoMaximaEVA(fm.model, fm.maxima[1:step:end,:,:], fm.parameters[1:step:end, :,:])
    
end


# Diagnostic plots for a given MCMC iteration

function diagnosticplots_iter(pme::PseudoMaximaEVA, iter::Int)

    fbm = convert(MaximumLikelihoodEVA, pme, iter)
    
    return diagnosticplots(fbm)
    
end

function histplot_iter(pme::PseudoMaximaEVA, iter::Int)

    fbm = convert(MaximumLikelihoodEVA, pme, iter)
    
    return histplot(fbm)
    
end

function probplot_iter(pme::PseudoMaximaEVA, iter::Int)

    fbm = convert(MaximumLikelihoodEVA, pme, iter)
    
    return probplot(fbm)
    
end

function qqplot_iter(pme::PseudoMaximaEVA, iter::Int)

    fbm = convert(MaximumLikelihoodEVA, pme, iter)
    
    return qqplot(fbm)
    
end


# Residual diagnostic plots

function diagnosticplots(fm::PseudoMaximaEVA; step::Int=1)
    
    f1 = probplot(fm, step = step)
    f2 = qqplot(fm, step = step)
    f3 = histplot(fm, step = step)
    f4 = plot()
    
    return gridstack([f1 f2; f3 f4])
    
end

function histplot(pme::PseudoMaximaEVA; step::Int=1)
    
    fm = ErrorsInVariablesExtremes.thin(pme, step)
    
    nsim, n, nchains = size(fm.maxima.value)
    
    if nsim*n > 3000
        @warn("Consider thinning the chains to have a more manageable plot. The function `probplot()` with the optional argument `step` could be used.")
    end
    
    Z = Matrix{Float64}(undef, nsim, n)

    for k in 1:nsim
        Z[k,:] = ErrorsInVariablesExtremes.standardize(fm.model, fm.maxima.value[k,:,1], fm.parameters.value[k,:,1])
    end

    z = vec(Z)
    
    nbin = floor(Int64, sqrt(length(z)))
    
    x, p = Extremes.ecdf(z)
    
    xmin, xmax = quantile.(Gumbel(),[1/1000, 1-1/1000])
    xp = range(xmin, xmax, length=1000)

    a = repeat([0.55, 0.85], outer=nbin)
    
    h = layer(x=x, Geom.histogram(bincount=nbin, density=true), alpha=[a;a])
    d = layer(x=xp, y=pdf.(Gumbel(), xp), Geom.line, Theme(default_color="red") )
    
    plot(d,h, 
        Coord.cartesian(xmin = xp[1], xmax = xp[end]),
        Guide.xlabel("Data"), Guide.ylabel("Density"),
        Guide.title("Residual density plot"))
    
end

function probplot(pme::PseudoMaximaEVA; step::Int=1)
    
    fm = ErrorsInVariablesExtremes.thin(pme, step)
    
    nsim, n, nchains = size(fm.maxima.value)
    
    if nsim*n > 3000
        @warn("Consider thinning the chains to have a more manageable plot. The function `probplot()` with the optional argument `step` could be used.")
    end
    
    Z = Matrix{Float64}(undef, nsim, n)

    for k in 1:nsim
        Z[k,:] = ErrorsInVariablesExtremes.standardize(fm.model, fm.maxima.value[k,:,1], fm.parameters.value[k,:,1])
    end

    z = vec(Z)
    
    y, p = Extremes.ecdf(z)

    plot(x=cdf.(Gumbel(),y), y=p, Geom.point, 
        Geom.abline(color="red", style=:dash),
        Guide.xlabel("Model"), Guide.ylabel("Empirical"), Guide.title("Residual probability plot"),
            Theme(discrete_highlight_color=c->nothing))
    
end

function qqplot(pme::PseudoMaximaEVA; step::Int=1)

    fm = ErrorsInVariablesExtremes.thin(pme, step)
        
    nsim, n, nchains = size(fm.maxima.value)
    
    if nsim*n > 3000
        @warn("Consider thinning the chains to have a more manageable plot. The function `probplot()` with the optional argument `step` could be used.")
    end

    Z = Matrix{Float64}(undef, nsim, n)

    for k in 1:nsim
        Z[k,:] = ErrorsInVariablesExtremes.standardize(fm.model, fm.maxima.value[k,:,1], fm.parameters.value[k,:,1])
    end

    z = vec(Z)
    y, p = Extremes.ecdf(z)

    plot(x=quantile.(Gumbel(),p), y=y, Geom.point, 
        Geom.abline(color="red", style=:dash),
        Guide.xlabel("Model"), Guide.ylabel("Empirical"), Guide.title("Residual quantile-quantile plot"),
            Theme(discrete_highlight_color=c->nothing))
    
end