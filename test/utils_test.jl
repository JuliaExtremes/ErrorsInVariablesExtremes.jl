
@testset "load_discharge_distribution" begin
        
    name = "SLSO00003"
    pdata = load_discharge_distribution(name)
    @test typeof(pdata) == Vector{Pseudodata}

    name = "unknown_file"

    @test_throws "There is no dataset with the name '$name'" load_discharge_distribution(name)
    
end

@testset "validateprior" begin
        
    prior = [Normal(), Normal()]

    @test_throws ErrorException ErrorsInVariablesExtremes.validateprior(prior, 1)
    
end