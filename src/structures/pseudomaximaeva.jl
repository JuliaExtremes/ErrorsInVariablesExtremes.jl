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
    dic(fm::PseudoMaximaEVA)

Compute the Deviance Information Criterion (DIC) described by Gelman et al. (2013) for the PseudoMaximaEVA model `fm`.

#### Details

Reference:
Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D.B., Vehtari, A. & Rubin, D.B. (2013). *Bayesian Data Analysis (3rd ed.)*. Chapman and Hall/CRC. https://doi.org/10.1201/b16018
"""
function dic(fm::PseudoMaximaEVA)
   
    Davg = mean(logpdf(fm))
    
    ŷ, θ̂ = ErrorsInVariablesExtremes.findposteriormode(fm)
    
    D = logpdf(fm.model, ŷ, θ̂)
    
    return 2*Davg - D
        
end


"""
    findposteriormode(fm::PseudoMaximaEVA)

Find the point (ŷ, θ̂) among the MCMC iterations of `fm` that maximized the log density.
"""
function findposteriormode(fm::PseudoMaximaEVA)
    
    nsim = size(fm.maxima.value,1)
    
    ll = Vector{Float64}(undef, nsim)
    
    for k in 1:nsim
        y = vec(fm.maxima.value[k,:,1])
        θ = vec(fm.parameters.value[k,:,1])
        
        ll[k] = logpdf(fm.model, y, θ)
    end
    
    ind = argmax(ll)
    
    ŷ = vec(fm.maxima.value[ind,:,1])
    θ̂ = vec(fm.parameters.value[ind,:,1])
    
    return ŷ, θ̂
    
end

"""
    getdistribution(fm::PseudoMaximaEVA)

Return the underlying extreme value distribution for each MCMC iterations.
"""
function getdistribution(fm::PseudoMaximaEVA)

    v = fm.parameters.value[:,:,1]

    V = Extremes.slicematrix(v, dims=2)

    D = ErrorsInVariablesExtremes.getdistribution.(fm.model, V)

    d = Extremes.unslicematrix(D, dims=2)

    return d

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

# Diagnostic plots

"""
    diagnosticplots(fm::PseudoMaximaEVA)

Diagnostic plots
"""
function diagnosticplots(fm::PseudoMaximaEVA)

    f1 = probplot(fm)
    f2 = qqplot(fm)
    f3 = histplot(fm)

    if ErrorsInVariablesExtremes.isstationary(fm.model)
        f4 = returnlevelplot(fm)
    else
        f4 = plot()
    end

    return gridstack([f1 f2; f3 f4])
end

function probplot(fm::PseudoMaximaEVA)
   
    if ErrorsInVariablesExtremes.isstationary(fm.model)
    
        p, y = ErrorsInVariablesExtremes.probplot_data(fm)
        
    else
        
        p, y = ErrorsInVariablesExtremes.probplot_std_data(fm)
        
    end
    
    ȳ = vec(mean(y, dims=1))
    
    plot(x=p, y=ȳ, Geom.point, 
        Geom.abline(color="red", style=:dash),
        Guide.xlabel("Model"), Guide.ylabel("Empirical"),
        Guide.title("Probability Plot"))
    
end

function probplot_data(fm::PseudoMaximaEVA)
   
    nsim, n, nchain  = size(fm.maxima.value)

    p = collect(1:n) / (n+1)
    
    y = Matrix{Float64}(undef, nsim, n)
  
    for k in 1:nsim
     
        θᵢ = fm.parameters.value[k, :, 1]
    
        pd = ErrorsInVariablesExtremes.getdistribution(fm.model, θᵢ)
        
        x = sort(fm.maxima.value[k,:,1])
    
        y[k,:] = cdf.(pd, x)
        
    end
    
    return p, y
    
end


function probplot_std_data(fm::PseudoMaximaEVA)
    
    nsim, n, nchain  = size(fm.maxima.value)

    p = collect(1:n)/ (n+1)

    y = Matrix{Float64}(undef, nsim, n)
   
    for k in 1:nsim
     
        θᵢ = fm.parameters.value[k, :, 1]
    
        pd = ErrorsInVariablesExtremes.getdistribution(fm.model, θᵢ)
        
        x = sort(fm.maxima.value[k,:,1])
        
        z = ErrorsInVariablesExtremes.standardize(fm.model, x, θᵢ)
    
        y[k,:] = cdf.(Gumbel(), z)
        
    end
    
    return p, y
    
end

function qqplot(fm::PseudoMaximaEVA)
    
    if ErrorsInVariablesExtremes.isstationary(fm.model)
    
        x, y = ErrorsInVariablesExtremes.qqplot_data(fm)
        
    else
        
        x, y = ErrorsInVariablesExtremes.qqplot_std_data(fm)
        
    end
    
    x̄ = vec(mean(x, dims=1))
    ȳ = vec(mean(y, dims=1))
    
    pmin = .9*minimum([x̄[1], ȳ[1]])
    pmax = 1.1*maximum([x̄[end], ȳ[end]])

    plot(x=x̄, y=ȳ, Geom.point, Geom.abline(color="red", style=:dash),
        Guide.xlabel("Model"), Guide.ylabel("Empirical"),
        Coord.cartesian(xmin=pmin, ymin=pmin),
        Guide.title("Quantile-Quantile Plot"))
    
end

function qqplot_data(fm::PseudoMaximaEVA)
    
    nsim, n, nchain  = size(fm.maxima.value)

    p = collect(1:n) / (n+1)

    x = Matrix{Float64}(undef, nsim, n)
    y = Matrix{Float64}(undef, nsim, n)
   
    for k in 1:nsim
   
        x[k,:] = sort(fm.maxima.value[k,:,1])
    
        θᵢ = fm.parameters.value[k, :, 1]
    
        pd = ErrorsInVariablesExtremes.getdistribution(fm.model, θᵢ)
    
        y[k,:] = quantile.(pd, p)
        
    end
    
    return x, y
    
end

function qqplot_std_data(fm::PseudoMaximaEVA)
    
    nsim, n, nchain  = size(fm.maxima.value)

    p = collect(1:n)/ (n+1)

    z = Matrix{Float64}(undef, nsim, n)
    y = Matrix{Float64}(undef, nsim, n)
   
    for k in 1:nsim
       
        θᵢ = fm.parameters.value[k, :, 1]
        
        x = sort(fm.maxima.value[k,:,1])
        
        z[k, :] = ErrorsInVariablesExtremes.standardize(fm.model, x, θᵢ)
    
        y[k,:] = quantile.(Gumbel(), p)
        
    end
    
    return z, y
    
end

function histplot(fm::PseudoMaximaEVA)

    if ErrorsInVariablesExtremes.isstationary(fm.model)
    
        x, xp, y = ErrorsInVariablesExtremes.histplot_data(fm)
        
    else
        
        x, xp, y = ErrorsInVariablesExtremes.histplot_std_data(fm)
        
    end
    
    x̄ = vec(mean(x, dims=1))
    ȳ = vec(mean(y, dims=1))
    
    nbin = floor(Int64, sqrt(size(x,2 )))

    h = layer(x=x̄, Geom.histogram(bincount = nbin, density=true))
    d = layer(x=xp, y=ȳ, Geom.line, Theme(default_color="red") )

    plot(d,h, 
        Coord.cartesian(xmin = xp[1], xmax = xp[end]),
        Guide.xlabel("Data"), Guide.ylabel("Density"),
        Guide.title("Density plot"))
    
end

function histplot_data(fm::PseudoMaximaEVA)
   
    nsim, n, nchain  = size(fm.maxima.value)
    
    xmin = quantile(fm.maxima.value[:], 1/1000)
    xmax = quantile(fm.maxima.value[:], 1 - 1/1000)
    xp = range(xmin, xmax, length=1000)
    
    x = Matrix{Float64}(undef, nsim, n)
    y = Matrix{Float64}(undef, nsim, length(xp))
    
    for k in 1:nsim
     
        x[k, :] = sort(fm.maxima.value[k, :, 1]) 
        
        θᵢ = fm.parameters.value[k, :, 1]
    
        pd = ErrorsInVariablesExtremes.getdistribution(fm.model, θᵢ)[]
        
        y[k,:] = pdf.(pd, xp)
        
    end
    
    return x, xp, y
    
end

function histplot_std_data(fm::PseudoMaximaEVA)
   
    nsim, n, nchain  = size(fm.maxima.value)
    
    xmin = quantile(Gumbel(), 1/1000)
    xmax = quantile(Gumbel(), 1 - 1/1000)
    xp = range(xmin, xmax, length=1000)
    
    z = Matrix{Float64}(undef, nsim, n)
    y = Matrix{Float64}(undef, nsim, length(xp))
    
    for k in 1:nsim
     
        θᵢ = fm.parameters.value[k, :, 1]
        
        x = sort(fm.maxima.value[k, :, 1])
        
        z[k,:] = ErrorsInVariablesExtremes.standardize(fm.model, x, θᵢ)
        
        y[k,:] = pdf.(Gumbel(), xp)
        
    end
    
    return z, xp, y
    
end

function returnlevelplot(fm::PseudoMaximaEVA)
    
    @assert ErrorsInVariablesExtremes.isstationary(fm.model) "The model should be stationary."
    
    T, x, y = ErrorsInVariablesExtremes.returnlevelplot_data(fm)

    x̄ = vec(mean(x, dims=1))
    ȳ = vec(mean(y, dims=1))

    l1 = layer(x=T, y=ȳ, Geom.line, Theme(default_color="red", line_style=[:dash]))
    l2 = layer(x=T, y=x̄, Geom.point)


    plot(l1,l2, Scale.x_log10, Guide.xlabel("Return Period"), Guide.ylabel("Return Level"),
                Guide.title("Return Level Plot"), Theme(discrete_highlight_color=c->nothing))
    
end


function returnlevelplot_data(fm::PseudoMaximaEVA)

    nsim, n, nchain  = size(fm.maxima.value)

    p = collect(1:n)/ (n+1)  
    
    T = 1 ./ (1 .- p)
    
    x = Matrix{Float64}(undef, nsim, n)
    y = Matrix{Float64}(undef, nsim, n)
    
    for k in 1:nsim
     
        x[k, :] = sort(fm.maxima.value[k, :, 1]) 
        
        θᵢ = fm.parameters.value[k, :, 1]
    
        pd = ErrorsInVariablesExtremes.getdistribution(fm.model, θᵢ)
        
        y[k,:] = quantile.(pd, p)
        
    end
    
    return T, x, y

end

