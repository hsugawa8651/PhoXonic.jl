#=
Plots.jl extension for PhoXonic.jl

This extension is automatically loaded when both PhoXonic and Plots are used.
=#

module PhoXonicPlotsExt

using PhoXonic
using Plots

import PhoXonic: plot_bands, plot_bands!, band_plot_data

"""
    plot_bands(bs::BandStructure; kwargs...)

Plot a band structure diagram.

# Keyword Arguments
- `color`: Line/marker color (default: :blue)
- `linewidth`: Line width (default: 2)
- `scatter`: Use scatter plot instead of lines (default: false).
  Recommended for 3D band structures where Γ point anomalies may occur.
- `markersize`: Marker size for scatter plot (default: 3)
- `markershape`: Marker shape for scatter plot (default: :circle)
- `title`: Plot title (default: "Band Structure")
- `ylabel`: Y-axis label (default: "Frequency")
- `xlabel`: X-axis label (default: "Wave vector")
- `show_gaps`: Highlight band gaps (default: false)
- `gap_color`: Color for gap highlighting (default: :yellow)
- `gap_alpha`: Transparency for gap highlighting (default: 0.2)
- `gap_threshold`: Minimum gap size to highlight (default: 0.0)
- `normalize`: Normalization factor for frequency (default: 1.0)
- `size`: Plot size (default: (600, 400))

# Example
```julia
# Line plot (default, good for 2D)
plot_bands(bands; color=:blue)

# Scatter plot (recommended for 3D)
plot_bands(bands; scatter=true, markersize=2, color=:blue)
```
"""
function PhoXonic.plot_bands(
    bs::BandStructure;
    color=:blue,
    linewidth::Real=2,
    scatter::Bool=false,
    markersize::Real=3,
    markershape::Symbol=:circle,
    title::String="Band Structure",
    ylabel::String="Frequency",
    xlabel::String="Wave vector",
    show_gaps::Bool=false,
    gap_color=:yellow,
    gap_alpha::Real=0.2,
    gap_threshold::Real=0.0,
    normalize::Real=1.0,
    size::Tuple{Int,Int}=(600, 400),
    kwargs...,
)
    data = band_plot_data(bs; normalize=normalize)

    # Create plot
    p = plot(;
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        legend=false,
        grid=true,
        size=size,
        kwargs...,
    )

    # Plot each band
    nbands = Base.size(data.frequencies, 2)
    for b in 1:nbands
        if scatter
            Plots.scatter!(
                p,
                data.distances,
                data.frequencies[:, b];
                markersize=markersize,
                markershape=markershape,
                color=color,
                markerstrokewidth=0,
            )
        else
            plot!(
                p, data.distances, data.frequencies[:, b]; linewidth=linewidth, color=color
            )
        end
    end

    # Add high-symmetry point markers
    if !isempty(data.label_positions)
        vline!(p, data.label_positions; color=:gray, linestyle=:dash, alpha=0.5)
        xticks!(p, data.label_positions, data.label_names)
    end

    # Highlight band gaps
    if show_gaps
        gaps = find_all_gaps(bs; threshold=gap_threshold)
        for g in gaps
            hspan!(
                p,
                [g.max_lower * normalize, g.min_upper * normalize];
                alpha=gap_alpha,
                color=gap_color,
                label="",
            )
        end
    end

    return p
end

"""
    plot_bands!(p, bs::BandStructure; kwargs...)

Add band structure to an existing plot.

# Keyword Arguments
- `color`: Line/marker color (default: :blue)
- `linewidth`: Line width (default: 2)
- `linestyle`: Line style (default: :solid)
- `scatter`: Use scatter plot instead of lines (default: false)
- `markersize`: Marker size for scatter plot (default: 3)
- `markershape`: Marker shape for scatter plot (default: :circle)
- `normalize`: Normalization factor for frequency (default: 1.0)
"""
function PhoXonic.plot_bands!(
    p,
    bs::BandStructure;
    color=:blue,
    linewidth::Real=2,
    linestyle=:solid,
    scatter::Bool=false,
    markersize::Real=3,
    markershape::Symbol=:circle,
    normalize::Real=1.0,
    kwargs...,
)
    data = band_plot_data(bs; normalize=normalize)

    nbands = Base.size(data.frequencies, 2)
    for b in 1:nbands
        if scatter
            Plots.scatter!(
                p,
                data.distances,
                data.frequencies[:, b];
                markersize=markersize,
                markershape=markershape,
                color=color,
                markerstrokewidth=0,
                kwargs...,
            )
        else
            plot!(
                p,
                data.distances,
                data.frequencies[:, b];
                linewidth=linewidth,
                color=color,
                linestyle=linestyle,
                kwargs...,
            )
        end
    end

    return p
end

# Plots.jl recipe for BandStructure
@recipe function f(bs::BandStructure)
    data = band_plot_data(bs)

    # Set default attributes
    xlabel --> "Wave vector"
    ylabel --> "Frequency"
    legend --> false
    grid --> true

    # Plot each band as a separate series
    for b in 1:Base.size(data.frequencies, 2)
        @series begin
            seriestype := :path
            linewidth --> 2
            data.distances, data.frequencies[:, b]
        end
    end
end

end # module
