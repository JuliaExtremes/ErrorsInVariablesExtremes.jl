
@testset "findposteriormode(::PseudoMaximaEVA)" begin

    Y = fill(100,3)

    η = log.(Y)
    ζ = 1/100

    pdata = Pseudodata("test", collect(1:length(Y)), LogNormal.(η, ζ))

    model = PseudoMaximaModel([pdata], prior=[Flat(), Flat(), Flat()])

    fm = PseudoMaximaEVA(model, 
        Mamba.Chains([Y'; Y' .+ 10 ;  Y' .+ 20]), 
        Mamba.Chains([100 log(10) -.1; 50 log(10) -.1; 150 log(10) -.1]))

    ŷ, θ̂ = ErrorsInVariablesExtremes.findposteriormode(fm)
    
    @test ŷ ≈ Y
    @test θ̂ ≈ [100, log(10), -.1]
    
end

@testset "convert(::Type{MaximumLikelihoodEVA}, fm::PseudoMaximaEVA, iter::Int)" begin
    
    y = [90., 100., 110.]
   
    pdata = Pseudodata("y", collect(0:2), Normal.(y, 1/100))
    
    Y = [y' .- .1 ; y' ; y' .+ .1]
    
    θ = [90 log(10) -.1; 100 log(10) -.1; 110 log(10) -.1]
    
    pmm = PseudoMaximaModel([pdata], prior=[Flat(), Flat(), Flat()])
    
    fm = PseudoMaximaEVA(pmm, 
        Mamba.Chains(Y), 
        Mamba.Chains(θ))
    
    res = convert(MaximumLikelihoodEVA, fm, 2)
    
    @test all(res.model.data.value .≈ Y[2,:])
    @test all(res.θ̂ .≈ θ[2,:])
    
end


@testset "dic(::PseudoMaximaEVA)" begin
   
    # Test dic() with 3 maxima and 3 MCMC iterations
    
    Y = fill(100,3)

    η = log.(Y)
    ζ = 1/100

    pdata = Pseudodata("test", collect(1:length(Y)), LogNormal.(η, ζ))

    model = PseudoMaximaModel([pdata], prior=[Flat(), Flat(), Flat()])

    fm = PseudoMaximaEVA(model, 
            Mamba.Chains([Y'; Y' .+ 10 ;  Y' .+ 20]), 
            Mamba.Chains([100 log(10) -.1; 50 log(10) -.1; 150 log(10) -.1]))

    res = 2*mean(logpdf(fm)) - logpdf(fm.model, Y, [100, log(10), -.1])
    
    @test dic(fm) ≈ res
    
end

@testset "getdistribution(::PseudoMaximaEVA)" begin
    
    y = [90., 100., 110.]
   
    pdata = Pseudodata("y", collect(0:2), Normal.(y, 1/100))
    
    Y = [y' .- .1 ; y' ; y' .+ .1]
    
    @testset "stationary model" begin
    
        θ = [90 log(10) -.1; 100 log(10) -.1; 110 log(10) -.1]

        pmm = PseudoMaximaModel([pdata], prior=[Flat(), Flat(), Flat()])

        fm = PseudoMaximaEVA(pmm, 
            Mamba.Chains(Y), 
            Mamba.Chains(θ))

        @test ErrorsInVariablesExtremes.getdistribution(fm, 2)[] ==
            GeneralizedExtremeValue(θ[2,1], exp(θ[2,2]), θ[2,3])

        @test all(ErrorsInVariablesExtremes.getdistribution(fm) .== 
            GeneralizedExtremeValue.(θ[:,1], exp.(θ[:,2]), θ[:,3]))
        
    end
    
    @testset "nonstationary model" begin
       
        θ = [80 10 log(10) -.1; 90 10 log(10) -.1; 100 10 log(10) -.1]

        pmm = PseudoMaximaModel([pdata], 
            locationcov = [Variable("x", collect(0:2))],
            prior=[Flat(), Flat(), Flat(), Flat()])

        fm = PseudoMaximaEVA(pmm, 
                Mamba.Chains(Y), 
                Mamba.Chains(θ))

        μ = θ[:,1] .+ θ[:,2].*collect(0:2)'
        σ = exp.(θ[:,3])
        ξ = θ[:,4]
        
        @test all(ErrorsInVariablesExtremes.getdistribution(fm, 2) .== 
            GeneralizedExtremeValue.(μ[2,:], σ[2], ξ[2]))
        
        @test all(ErrorsInVariablesExtremes.getdistribution(fm) .== 
            GeneralizedExtremeValue.(μ, σ[2], ξ[2]))
        
    end
    
end

@testset "logpdf(::PseudoMaximaEVA)" begin

    Y = fill(100,3)

    η = log.(Y)
    ζ = 1/100

    pdata = Pseudodata("test", collect(1:length(Y)), LogNormal.(η, ζ))

    model = PseudoMaximaModel([pdata], prior=[Flat(), Flat(), Flat()])

    fm = PseudoMaximaEVA(model, 
            Mamba.Chains([Y'; Y' .+ 10 ;  Y' .+ 20]), 
            Mamba.Chains([100 log(10) -.1; 50 log(10) -.1; 150 log(10) -.1]))
    
    res = logpdf(fm)
    
    @test res[1] ≈ sum(logpdf(pdata, Y) + logpdf.(GeneralizedExtremeValue(100, 10, -.1), Y))
    @test res[2] ≈ sum(logpdf(pdata, Y .+ 10) + logpdf.(GeneralizedExtremeValue(50, 10, -.1), Y .+ 10))
    @test res[3] ≈ sum(logpdf(pdata, Y .+ 20) + logpdf.(GeneralizedExtremeValue(150, 10, -.1), Y .+ 20))
    
end

# Tests on diagnostic plots

@testset "diagnostic plots for stationary model" begin
   
    y = [90., 100., 110.]
    
    n = length(y)
    
    Y = [y' .- 1 ; y' ; y' .+ 1]
    
    p₀ = [.25 .5 .75;
         .25 .5 .75;
         .25 .5 .75]

    Θ = [90 log(10) -.1;
        100 log(10) -.1;
        110 log(10) -.1]

    pdata = Pseudodata("y", collect(1:3), Normal.(y, 1/100))

    pmm = PseudoMaximaModel([pdata], prior=[Flat(), Flat(), Flat()])

    fmm = PseudoMaximaEVA(pmm, Mamba.Chains(Y), Mamba.Chains(Θ))
    
    pd = GeneralizedExtremeValue.(Θ[:,1], exp.(Θ[:,2]), Θ[:,3])
    
    @testset "probplot_data" begin

        p, p̂ = ErrorsInVariablesExtremes.probplot_data(fmm)

        @test all(p .≈ p₀[1,:])
        @test all(p̂ .≈ cdf.(pd, Y))
    end
    
    @testset "qqplot_data" begin

        q, q̂ = ErrorsInVariablesExtremes.qqplot_data(fmm)

        @test all(q .≈ Y)
        @test all(q̂ .≈ quantile.(pd, p₀))
    end
    
    @testset "histplot_data" begin

        d, xp, d̂ = ErrorsInVariablesExtremes.histplot_data(fmm)

        @test all(d .≈ Y) 

        @test all(d̂ .≈ pdf.(pd, collect(xp')))
    end
    
    @testset "returnlevel_data" begin

        T, r, r̂ = ErrorsInVariablesExtremes.returnlevelplot_data(fmm)

        @test all(T .≈ 1 ./ (1 .- p₀[1,:]))

        @test all(r .≈ Y)

        @test all(r̂ .≈ quantile.(pd, p₀ ))
    end

    
end


@testset "thin(::PseudoMaximaEVA, step::Int)" begin
    
    y = [90., 100., 110.]
   
    pdata = Pseudodata("y", collect(0:2), Normal.(y, 1/100))
    
    Y = [y' .- .1 ; y' ; y' .+ .1]
    
    θ = [90 log(10) -.1; 100 log(10) -.1; 110 log(10) -.1]
    
    pmm = PseudoMaximaModel([pdata], prior=[Flat(), Flat(), Flat()])
    
    fm = PseudoMaximaEVA(pmm, 
        Mamba.Chains(Y), 
        Mamba.Chains(θ))
    
    res = ErrorsInVariablesExtremes.thin(fm, 2)
    
    @test all(res.maxima.value .≈ Y[1:2:end, :])
    @test all(res.parameters.value .≈ θ[1:2:end, :])
    
end



@testset "diagnostic plots for nonstationary model" begin
   
    y = [90., 100., 110.]
    
    n = length(y)
    
    Y = [y' .- 1 ; y' ; y' .+ 1]
    
    p₀ = [.25 .5 .75;
         .25 .5 .75;
         .25 .5 .75]

    Θ = [80 5 log(10) -.1;
        90 5 log(10) -.1;
        100 5 log(10) -.1]
    
    nsim = 3

    pdata = Pseudodata("y", collect(1:3), Normal.(y, 1/100))

    pmm = PseudoMaximaModel([pdata],
        locationcov = [Variable("x", collect(0:2))],
        prior=[Flat(), Flat(), Flat(), Flat()])

    fmm = PseudoMaximaEVA(pmm, Mamba.Chains(Y), Mamba.Chains(Θ))
    
    pd = ErrorsInVariablesExtremes.getdistribution(fmm)
    
    Z = Matrix{Float64}(undef, nsim, n)
    for k in 1:nsim
        Z[k,:] = ErrorsInVariablesExtremes.standardize(fmm.model, Y[k,:], Θ[k,:])
    end
    
    @testset "probplot_std_data" begin

        p, p̂ = ErrorsInVariablesExtremes.probplot_std_data(fmm)

        @test all(p .≈ p₀[1,:])
        
        @test all(p̂ .≈ cdf.(Gumbel(), Z))
    end
    
    @testset "qqplot_std_data" begin

        q, q̂ = ErrorsInVariablesExtremes.qqplot_std_data(fmm)

        @test all(q .≈ Z)
        @test all(q̂ .≈ quantile.(Gumbel(), p₀))
    end
    
    @testset "histplot_std_data" begin

        d, xp, d̂ = ErrorsInVariablesExtremes.histplot_std_data(fmm)

        @test all(d .≈ Z) 

        @test all(d̂ .≈ pdf.(Gumbel(), repeat(xp', 3)))
    end
    
end