
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

        μ = θ[:,1]
        σ = exp.(θ[:,2])
        ξ = θ[:,3]

        @test all(ErrorsInVariablesExtremes.getdistribution(fm) .== GeneralizedExtremeValue.(μ, σ, ξ))
        
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

        @test all(ErrorsInVariablesExtremes.getdistribution(fm) .== GeneralizedExtremeValue.(μ, σ, ξ))
        
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