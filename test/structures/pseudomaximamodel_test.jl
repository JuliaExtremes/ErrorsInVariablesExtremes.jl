# @testset "fitbayes(::PseudoMaximaModel)" begin
   
#     n = 100
#     t = collect(1:n)

#     Y = rand(GeneralizedExtremeValue(100, 10, -.1), n)

#     η = log.(Y)
#     ζ = 1/1000

#     pdata = Pseudodata("test", t, LogNormal.(η, ζ))
    
#     @testset "One member fit" begin
    
#         model = PseudoMaximaModel([pdata])

#         fm = fitbayes(model)

#         V = Extremes.parametervar(fm.parameters)
#         δ = 3*sqrt.(diag(V))

        # θ̂ = Extremes.findposteriormode(fm.parameters)

        # @test θ̂[1] ≈ 100 atol = δ[1]
        # @test θ̂[2] ≈ log(10) atol = δ[2]
        # @test θ̂[3] ≈ -.1 atol = δ[3]

        # @test fm.parameters.model.data.value ≈ Y rtol=.01
    # end
    
    # @testset "Fit on several members" begin
    
    #     pensemble = Pseudoensemble("test", [pdata, pdata, pdata])
    #     fm = gevfitbayes(pensemble)

    #     V = Extremes.parametervar(fm.parameters)
    #     δ = 3*sqrt.(diag(V))

    #     θ̂ = Extremes.findposteriormode(fm.parameters)

    #     @test θ̂[1] ≈ 100 atol = δ[1]
    #     @test θ̂[2] ≈ log(10) atol = δ[2]
    #     @test θ̂[3] ≈ -.1 atol = δ[3]

    #     @test fm.parameters.model.data.value ≈ Y rtol=.01
    # end
    
# end

@testset "logpdf(::PseudoMaximaModel, y, θ)" begin
    
    @testset "stationary model" begin
        
        n = 100
        t = collect(1:n)

        Y = rand(GeneralizedExtremeValue(100, 10, -.1), n)

        η = log.(Y)
        ζ = 1/1000

        pdata = Pseudodata("test", t, LogNormal.(η, ζ))
        
        @testset "Flat priors" begin
            
            model = PseudoMaximaModel([pdata], prior=[Flat(), Flat(), Flat()])

            @test logpdf(model, Y, [100, log(10), -.1]) ≈ 
                sum(logpdf(pdata, Y) + logpdf.(GeneralizedExtremeValue(100, 10, -.1), Y))
        end
        
        @testset "Custom priors" begin
            
            model = PseudoMaximaModel([pdata], prior=[Flat(), Flat(), Normal()])

            @test logpdf(model, Y, [100, log(10), -.1]) ≈ 
                sum(logpdf(pdata, Y) + logpdf.(GeneralizedExtremeValue(100, 10, -.1), Y)) + logpdf(Normal(), -.1)
        end
        
    end
    
    @testset "nonstationary model" begin
       
        n = 100
        t = collect(1:n)

        μ = 100 .+ t

        Y = rand.(GeneralizedExtremeValue.(μ, 10, -.1))

        η = log.(Y)
        ζ = 1/1000

        pdata = Pseudodata("test", t, LogNormal.(η, ζ))

        model = PseudoMaximaModel([pdata],locationcov = [Variable("t", t)], prior=[Flat(), Flat(), Flat(), Flat()])

        @test logpdf(model, Y, [100, 1, log(10), -.1]) ≈ 
            sum(logpdf(pdata, Y) + logpdf.(GeneralizedExtremeValue.(μ, 10, -.1), Y))
        
    end
    
end

