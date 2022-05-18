using DataFrames, Dates, NetCDF
using Distributions, ErrorsInVariablesExtremes, Extremes
using Test
using LinearAlgebra, Random
using Mamba
using Statistics

# Set the seed for reproductible test results
Random.seed!(12)

@testset "ErrorsInVariablesExtremes.jl" begin
    include("structures_test.jl")
    include("utils_test.jl")
end;
