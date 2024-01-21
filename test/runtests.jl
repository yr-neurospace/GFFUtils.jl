using GFFUtils
using Test
using DataFrames

@testset "All tests" begin
    include("gff_io_tests.jl")
    include("gff_utils_tests.jl")
end