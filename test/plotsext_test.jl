using Test
using PhoXonic
using Plots
using StaticArrays

@testset "PhoXonicPlotsExt smoke" begin
    bs = PhoXonic.BandStructure{2}(
        [SVector(0.0, 0.0), SVector(0.5, 0.0)],
        [0.0, 1.0],
        Float64[
            0.10 0.30
            0.20 0.40
        ],
        [(1, "Γ"), (2, "X")],
    )
    p = plot_bands(bs)
    @test p isa Plots.Plot
end
