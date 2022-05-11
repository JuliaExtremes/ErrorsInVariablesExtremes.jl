# Test on Structures

@testset "Pseudodata constructor" begin
    pdata = Pseudodata("Normal", collect(1:10), Normal.(1:10, 1))
    @test pdata.name == "Normal"
    @test pdata.year == collect(1:10)
    @test pdata.value == Normal.(1:10, 1)
end


# Tests on Pseudodata Methods

@testset "Pseudodata probability calculations" begin
    pdata = Pseudodata("Normal", collect(1:10), Normal.(zeros(10), 1))
    
    @test_throws MethodError logpdf(pdata, 0)
    @test_throws MethodError pdf(pdata, 0)
    
    @test all( pdf(pdata, zeros(10)) .≈ pdf(Normal(), 0) )
    @test all( logpdf(pdata, zeros(10)) .≈ logpdf(Normal(), 0) )
    
end