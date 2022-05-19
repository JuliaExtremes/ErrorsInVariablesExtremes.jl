
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