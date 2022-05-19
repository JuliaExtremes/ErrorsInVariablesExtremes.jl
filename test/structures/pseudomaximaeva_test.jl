
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


# @testset "dic(::PseudoMaximaEVA)" begin
   
#     # Test dic() with 2 maxima and two MCMC iterations
    
#     model = BlockMaxima(Variable("y", [.5, 1]))

#     fm = BayesianEVA(model, Mamba.Chains([0 0 -.1; -1 0 -.1], names=["μ", "ϕ", "ξ"]))

#     p₁ = Pseudodata("y₁", [1] ,[Normal(0,1)])
#     p₂ = Pseudodata("y₂", [1] ,[Normal(1,1)])

#     pensemble = Pseudoensemble("test", [p₁, p₂])

#     eiv_model = PseudoMaximaEVA(pensemble, fm, Mamba.Chains([.5 1; 0 .5], names=["Y[1]", "Y[2]"]))

#     ŷ = vec(mean([.5 1; 0 .5], dims=1))
#     θ̂ = vec(mean([0 0 -.1; -1 0 -.1], dims=1))
    
#     @test dic(eiv_model) ≈ (2*mean(loglike(eiv_model)) - loglike(eiv_model, ŷ, θ̂))
    
# end

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