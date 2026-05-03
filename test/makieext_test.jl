using Test
using PhoXonic
using CairoMakie
using StaticArrays

@testset "PhoXonicMakieExt smoke" begin
    bs = PhoXonic.BandStructure{2}(
        [SVector(0.0, 0.0), SVector(0.5, 0.0)],
        [0.0, 1.0],
        Float64[
            0.10 0.30
            0.20 0.40
        ],
        [(1, "Γ"), (2, "X")],
    )
    fig = plot_bands(bs)
    @test fig isa Makie.Figure
end
