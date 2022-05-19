
@testset "load_discharge_distribution" begin
        
    filename = "data/A2020_Analyse_Historique_QMA_SLSO00003.nc"
    pdata = load_discharge_distribution(filename)
    
    data = ncread(filename, "Dis")

    p = collect(.05:.05:.95)
    x = data[2:end-1,:,:] # I don't know what the first and last values are.
    
    for i in size(data, 3)
        for j in size(data, 2)

            pd = pdata[i].value[j]

            q = quantile.(pd, p)
            
            @test x[:,j,i] ≈ q rtol=1e-6
            
        end
    end
    
end

@testset "validateprior" begin
        
    prior = [Normal(), Normal()]

    @test_throws ErrorException ErrorsInVariablesExtremes.validateprior(prior, 1)
    
end