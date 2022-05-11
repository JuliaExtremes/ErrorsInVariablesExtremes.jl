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