# Last-Modified: 2025-12-12T22:52:30+09:00

#=
Plotting utilities for PhoXonic.jl

This module provides convenient plotting functions for band structures.
Requires Plots.jl to be loaded separately.

Usage:
    using PhoXonic
    using Plots

    bands = compute_bands(solver, kpath)
    plot_bands(bands)
=#

"""
    plot_bands(bs::BandStructure; kwargs...)

Plot a band structure diagram.

# Arguments
- `bs`: BandStructure object from compute_bands()

# Keyword Arguments
- `color`: Line/marker color (default: :blue)
- `linewidth`: Line width (default: 2)
- `scatter`: Use scatter plot instead of lines (default: false).
  Recommended for 3D band structures where Γ point anomalies may occur.
- `markersize`: Marker size for scatter plot (default: 3)
- `markershape`: Marker shape for scatter plot (default: :circle)
- `title`: Plot title (default: "Band Structure")
- `ylabel`: Y-axis label (default: "Frequency")
- `show_gaps`: Highlight band gaps (default: false)
- `gap_color`: Color for gap highlighting (default: :yellow)
- `gap_alpha`: Transparency for gap highlighting (default: 0.2)
- `normalize`: Normalization factor for frequency (default: 1.0)

# Returns
A Plots.jl plot object.

# Example
```julia
using PhoXonic, Plots
solver = Solver(TEWave(), geo, (64, 64))
bands = compute_bands(solver, kpath)
p = plot_bands(bands; title="TE Bands", show_gaps=true)
savefig(p, "bands.png")

# For 3D, use scatter to avoid line artifacts at Γ point
p = plot_bands(bands; scatter=true, markersize=2)
```
"""
function plot_bands end

"""
    plot_bands!(p, bs::BandStructure; kwargs...)

Add band structure to an existing plot.
"""
function plot_bands! end

# Note: Actual implementations are in ext/PhoXonicPlotsExt.jl
# These stubs only exist to provide docstrings and allow pre-declaration

"""
    band_plot_data(bs::BandStructure; normalize=1.0)

Extract data for plotting a band structure.

Returns a NamedTuple with:
- `distances`: K-path distances
- `frequencies`: Normalized frequency matrix
- `label_positions`: Positions of high-symmetry points
- `label_names`: Names of high-symmetry points
"""
function band_plot_data(bs::BandStructure; normalize::Real=1.0)
    dists = bs.distances
    freqs = bs.frequencies * normalize

    label_positions = Float64[]
    label_names = String[]

    for (idx, name) in bs.labels
        push!(label_positions, dists[idx])
        push!(label_names, name)
    end

    (
        distances=dists,
        frequencies=freqs,
        label_positions=label_positions,
        label_names=label_names,
    )
end
