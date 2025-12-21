# TMM module tests

using Test
using PhoXonic

@testset "TMM" begin
    include("test_layers.jl")
    include("test_matrices.jl")
    include("test_solver.jl")
    include("test_bands.jl")
    include("test_oblique.jl")
    include("test_phononic.jl")
end
