@testset "structures.jl" begin
    include(joinpath("structures", "pseudodata_test.jl"))
    include(joinpath("structures", "pseudoensemble_test.jl"))
    include(joinpath("structures", "pseudomaximaeva_test.jl"))
end