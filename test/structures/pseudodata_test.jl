# Test on Structures

@testset "Pseudodata constructor" begin
    pdata = Pseudodata("Normal", collect(1:10), Normal.(1:10, 1))
    @test pdata.name == "Normal"
    @test pdata.year == collect(1:10)
    @test pdata.value == Normal.(1:10, 1)
end


# Tests on Pseudodata Methods

@testset "ensemblemean(::Vector{Pseudodata})" begin
   
    pdata1 = Pseudodata("member 1", [1, 2], [Normal(0,1), Normal(1,2)])
    pdata2 = Pseudodata("member 2", [1, 2], [Normal(0,1), Normal(1,2)])
    pdata = [pdata1, pdata2]
    
    m = ensemblemean(pdata)
    
    @test m[1] ≈ 0
    @test m[2] ≈ 1
    
end

@testset "Pseudodata probability calculations" begin
    pdata = Pseudodata("Normal", collect(1:10), Normal.(zeros(10), 1))
    
    @test_throws MethodError logpdf(pdata, 0)
    @test_throws MethodError pdf(pdata, 0)
    
    @test all( pdf(pdata, zeros(10)) .≈ pdf(Normal(), 0) )
    @test all( logpdf(pdata, zeros(10)) .≈ logpdf(Normal(), 0) )
    
end

@testset "logpdf(::Vector{Pseudodata})" begin
   
    pdata1 = Pseudodata("member 1", [1, 2], [Normal(0,1), Normal(1,2)])
    pdata2 = Pseudodata("member 2", [1, 2], [Normal(0,1), Normal(1,2)])
    pdata = [pdata1, pdata2]
    
    @test logpdf(pdata, 0, 1) ≈ logpdf(Normal(0,1), 0) + logpdf(Normal(0,1), 0)
    @test logpdf(pdata, 0, 2) ≈ logpdf(Normal(1,2), 0) + logpdf(Normal(1,2), 0)
    
    ll = logpdf(pdata, [0, 0])
    @test ll[1] ≈ logpdf(Normal(0,1), 0) + logpdf(Normal(0,1), 0)
    @test ll[2] ≈ logpdf(Normal(1,2), 0) + logpdf(Normal(1,2), 0)
    
end

