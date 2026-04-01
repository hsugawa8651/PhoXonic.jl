#=
RecipesBase extension for PhoXonic.jl

This extension is automatically loaded when RecipesBase is available.
Provides a lightweight @recipe for BandStructure plotting.
=#

module PhoXonicRecipesBaseExt

using PhoXonic
using RecipesBase

import PhoXonic: band_plot_data

# Recipe for plot(bs::BandStructure) — works with any Plots.jl backend
@recipe function f(bs::BandStructure)
    data = band_plot_data(bs)

    xlabel --> "Wave vector"
    ylabel --> "Frequency"
    legend --> false
    grid --> true

    for b in 1:Base.size(data.frequencies, 2)
        @series begin
            seriestype := :path
            linewidth --> 2
            data.distances, data.frequencies[:, b]
        end
    end
end

end # module
