@testset "fitbayes(::PseudoMaximaModel)" begin
   
    n = 100
    t = collect(1:n)

    Y = rand(GeneralizedExtremeValue(100, 10, -.1), n)
    
    eva_fm = gevfit(Y)

    η = log.(Y)
    ζ = 1/1000

    pdata = Pseudodata("test", t, LogNormal.(η, ζ))
    
    @testset "One member fit" begin
    
        model = PseudoMaximaModel([pdata])

        fm = fitbayes(model)

        V = Extremes.parametervar(eva_fm)
        δ = 3*sqrt.(diag(V))

        ŷ = vec(mean(fm.maxima.value[:,:,1], dims=1))
        θ̂ = vec(mean(fm.parameters.value[:,:,1], dims=1))

        @test θ̂[1] ≈ 100 atol = δ[1]
        @test θ̂[2] ≈ log(10) atol = δ[2]
        @test θ̂[3] ≈ -.1 atol = δ[3]

        @test ŷ ≈ Y rtol=.01
    end
    
    @testset "Fit on several members" begin
    
        model = PseudoMaximaModel([pdata, pdata])

        fm = fitbayes(model)

        V = Extremes.parametervar(eva_fm)
        δ = 3*sqrt.(diag(V))

        ŷ = vec(mean(fm.maxima.value[:,:,1], dims=1))
        θ̂ = vec(mean(fm.parameters.value[:,:,1], dims=1))

        @test θ̂[1] ≈ 100 atol = δ[1]
        @test θ̂[2] ≈ log(10) atol = δ[2]
        @test θ̂[3] ≈ -.1 atol = δ[3]

        @test ŷ ≈ Y rtol=.01
        
    end
    
end

@testset "getdistribution(::PseudoMaximaModel)" begin

    y = [90., 100., 110.]

    pdata = Pseudodata("y", collect(0:2), Normal.(y, 1/100))
    
    @testset "stationary model" begin
       
        model = PseudoMaximaModel([pdata],
        prior=[Flat(), Flat(), Flat()])

        θ = [100., log(10), -.1]

        @test ErrorsInVariablesExtremes.getdistribution(model, θ)[] .== GeneralizedExtremeValue(θ[1], exp(θ[2]), θ[3])
        
    end
    
    @testset "nonstationary model" begin
       
        model = PseudoMaximaModel([pdata], locationcov = [Variable("x", collect(0:2))],
        prior=[Flat(), Flat(), Flat(), Flat()])

        θ = [90., 10., log(10), -.1]

        @test all(ErrorsInVariablesExtremes.getdistribution(model, θ) .== GeneralizedExtremeValue.(θ[1] .+ θ[2].*collect(0:2), exp(θ[3]), θ[4]))
        
    end
    
end

@testset "isstationary(::PseudoMaximaModel)" begin

    @testset "stationary model" begin
    
        pm = PseudoMaximaModel([Pseudodata("y", [1], [Normal()])],
            prior = [Flat(), Flat(), Flat()])
    
        @test ErrorsInVariablesExtremes.isstationary(pm)
    
    end

    @testset "nonstationary model" begin
    
        pm = PseudoMaximaModel([Pseudodata("y", [1], [Normal()])],
            locationcov = [Variable("x", [1])],
            prior = [Flat(), Flat(), Flat(), Flat()])
    
        @test !(ErrorsInVariablesExtremes.isstationary(pm))
    
    end

end

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

@testset "standardize(::PseudoMaximaModel)" begin
   
    @testset "stationary model" begin
       
        pdata = Pseudodata("Empty", Int64[], UnivariateDistribution[])
        model = PseudoMaximaModel([pdata], prior=[Flat(), Flat(), Flat()])

        y = [90., 100., 110.]
        θ = [100., log(10), -.1]

        @test all(ErrorsInVariablesExtremes.standardize(model, y, θ) .≈ Extremes.standardize.(y, θ[1], exp(θ[2]), θ[3]))

    end
    
    @testset "nonstationary model" begin
       
        y = [90., 100., 110.]

        pdata = Pseudodata("Empty", collect(0:2), Normal.(y, 1/100))
        model = PseudoMaximaModel([pdata], locationcov = [Variable("x", collect(0:2))],
            prior=[Flat(), Flat(), Flat(), Flat()])

        θ = [90., 10., log(10), -.1]

        @test all(ErrorsInVariablesExtremes.standardize(model, y, θ) .≈ Extremes.standardize.(y, θ[1] .+ θ[2].*collect(0:2), exp(θ[3]), θ[4]))
        
    end
    
end
