# Test on Structures

@testset "Pseudoensemble constructor" begin
    pdata = Pseudodata("Normal", collect(1:10), Normal.(zeros(10), 1))
    pensemble = Pseudoensemble("Ensemble", [pdata for i=1:5])
    
    @test pensemble.name == "Ensemble"
    @test pensemble.value[1].name == "Normal"
    
    @test pensemble.value[1] == pdata
    
end


# Tests on Pseudoensemble Methods

@testset "Pseudoensemble probability calculations" begin
    pdata = Pseudodata("Normal", collect(1:10), Normal.(zeros(10), 1))
    pensemble = Pseudoensemble("Ensemble", [pdata for i=1:5])
    
    @test_throws MethodError logpdf(pensemble, 0)
    @test_throws MethodError pdf(pensemble, 0)
    
    @test ensemblemean(pensemble) ≈ zeros(10)
    @test logpdf(pensemble, 0, 1) ≈ 5*logpdf(Normal(), 0)
    @test all( logpdf(pensemble, zeros(10)) .≈ 5*logpdf(Normal(), 0) )
    @test all( pdf(pensemble, zeros(10)) .≈ (pdf(Normal(), 0))^5 )
    
end

@testset "Pseudoensemble gevfitbayes" begin
   
    n = 100
    t = collect(1:n)

    Y = rand(GeneralizedExtremeValue(100, 10, -.1), n)

    η = log.(Y)
    ζ = 1/1000

    pdata = Pseudodata("test", t, LogNormal.(η, ζ))
    
    @testset "One member fit" begin
    
        pensemble = Pseudoensemble("test", [pdata])
        fm = gevfitbayes(pensemble)

        V = Extremes.parametervar(fm.parameters)
        δ = 3*sqrt.(diag(V))

        θ̂ = Extremes.findposteriormode(fm.parameters)

        @test θ̂[1] ≈ 100 atol = δ[1]
        @test θ̂[2] ≈ log(10) atol = δ[2]
        @test θ̂[3] ≈ -.1 atol = δ[3]

        @test fm.parameters.model.data.value ≈ Y rtol=.01
    end
    
    @testset "Fit on several members" begin
    
        pensemble = Pseudoensemble("test", [pdata, pdata, pdata])
        fm = gevfitbayes(pensemble)

        V = Extremes.parametervar(fm.parameters)
        δ = 3*sqrt.(diag(V))

        θ̂ = Extremes.findposteriormode(fm.parameters)

        @test θ̂[1] ≈ 100 atol = δ[1]
        @test θ̂[2] ≈ log(10) atol = δ[2]
        @test θ̂[3] ≈ -.1 atol = δ[3]

        @test fm.parameters.model.data.value ≈ Y rtol=.01
    end
    
end