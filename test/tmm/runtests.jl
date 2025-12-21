# TMM module tests

using Test
using PhoXonic

@testset "TMM" begin
    include("test_layers.jl")
    include("test_matrices.jl")
end
