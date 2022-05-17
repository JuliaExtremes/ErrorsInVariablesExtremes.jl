
@testset "loglike(::PseudoMaximaEVA)" begin
   
    # Test for 2 maxima and two MCMC iterations
    
    model = BlockMaxima(Variable("y", [.5, 1]))

    fm = BayesianEVA(model, Mamba.Chains([0 0 -.1; -1 0 -.1], names=["μ", "ϕ", "ξ"]))

    p₁ = Pseudodata("y₁", [1] ,[Normal(0,1)])
    p₂ = Pseudodata("y₂", [1] ,[Normal(1,1)])

    pensemble = Pseudoensemble("test", [p₁, p₂])

    eiv_model = PseudoMaximaEVA(pensemble, fm, Mamba.Chains([.5 1; 0 .5], names=["Y[1]", "Y[2]"]))
    
    res = loglike(eiv_model)
    
    @testset "loglike values" begin
        @test res[1] ≈ sum(logpdf.(Normal(0,1), [.5, 1]) + logpdf.(Normal(1,1), [.5, 1]) + logpdf.(GeneralizedExtremeValue(0,1,-.1), [.5, 1]))
        @test res[2] ≈ sum(logpdf.(Normal(0,1), [0, .5]) + logpdf.(Normal(1,1), [0, .5]) + logpdf.(GeneralizedExtremeValue(-1,1,-.1), [0, .5]))
    end
    
end