
@testset "quantile2gaussian" begin
    
    @test_throws AssertionError ErrorsInVariablesExtremes.quantile2gaussian([0, 0], [.5, .5])
    @test_throws AssertionError ErrorsInVariablesExtremes.quantile2gaussian([0, 0], [.5, .5, .5])
    @test_throws AssertionError ErrorsInVariablesExtremes.quantile2gaussian([0], [.5])
    
    pd = Normal(1,2)
    
    p = [.25, .75]
    x = quantile.(pd, p)
    
    @test ErrorsInVariablesExtremes.quantile2gaussian(x,p) == pd
    
end