#=
Plots.jl extension for PhoXonic.jl

This extension is automatically loaded when both PhoXonic and Plots are used.
=#

module PhoXonicPlotsExt

using PhoXonic
using Plots

import PhoXonic: plot_bands, plot_bands!, band_plot_data
import PhoXonic: plot_field, plot_field!, plot_epsilon
import PhoXonic: reconstruct_field, get_epsilon_field, fix_phase
import PhoXonic: Dim1, Dim2, Dim3, Solver, WaveType

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

# =============================================================================
# Field Visualization Functions
# =============================================================================

# -----------------------------------------------------------------------------
# plot_field - 1D
# -----------------------------------------------------------------------------
"""
    plot_field(solver::Solver{Dim1}, eigenvector; kwargs...)

Plot a 1D field distribution.

# Keyword Arguments
- `quantity::Symbol = :real` - Quantity to plot: `:real`, `:imag`, `:abs`, `:phase`
- `fix_phase_flag::Bool = true` - Apply phase normalization
- `title::String = "Field"` - Plot title
- `xlabel::String = "x/a"` - X-axis label
- `ylabel::String = "Field"` - Y-axis label
"""
function PhoXonic.plot_field(
    solver::Solver{Dim1},
    eigenvector::AbstractVector;
    quantity::Symbol=:real,
    fix_phase_flag::Bool=true,
    title::String="Field",
    xlabel::String="x/a",
    ylabel::String="Field",
    kwargs...
)
    field = reconstruct_field(solver, eigenvector)
    if fix_phase_flag
        field = fix_phase(field)
    end

    data = _extract_quantity(field, quantity)
    N = length(data)
    x = range(0, 1, length=N)

    Plots.plot(x, data;
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        legend=false,
        kwargs...
    )
end

# -----------------------------------------------------------------------------
# plot_field - 2D
# -----------------------------------------------------------------------------
"""
    plot_field(solver::Solver{Dim2}, eigenvector; kwargs...)

Plot a 2D field distribution as a heatmap.

# Keyword Arguments
- `component::Symbol = :auto` - Component to plot for vector fields
- `quantity::Symbol = :real` - Quantity to plot: `:real`, `:imag`, `:abs`, `:phase`
- `fix_phase_flag::Bool = true` - Apply phase normalization
- `colormap::Symbol = :RdBu` - Colormap
- `title::String = "Field"` - Plot title
- `xlabel::String = "x/a"` - X-axis label
- `ylabel::String = "y/a"` - Y-axis label
"""
function PhoXonic.plot_field(
    solver::Solver{Dim2,W},
    eigenvector::AbstractVector;
    component::Symbol=:auto,
    quantity::Symbol=:real,
    fix_phase_flag::Bool=true,
    colormap::Symbol=:RdBu,
    title::String="Field",
    xlabel::String="x/a",
    ylabel::String="y/a",
    kwargs...
) where {W}
    field = reconstruct_field(solver, eigenvector)

    # Handle tuple (vector field) or matrix (scalar field)
    if field isa Tuple
        comp_idx = _get_component_index(W, component)
        data_complex = field[comp_idx]
    else
        data_complex = field
    end

    if fix_phase_flag
        data_complex = fix_phase(data_complex)
    end

    data = _extract_quantity(data_complex, quantity)

    # Transpose for correct orientation (x horizontal, y vertical)
    Plots.heatmap(data';
        c=colormap,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        aspect_ratio=:equal,
        kwargs...
    )
end

# -----------------------------------------------------------------------------
# plot_field! - Add field to existing plot
# -----------------------------------------------------------------------------
"""
    plot_field!(p, solver, eigenvector; kwargs...)

Add a field to an existing plot.
"""
function PhoXonic.plot_field!(
    p,
    solver::Solver{Dim1},
    eigenvector::AbstractVector;
    quantity::Symbol=:real,
    fix_phase_flag::Bool=true,
    kwargs...
)
    field = reconstruct_field(solver, eigenvector)
    if fix_phase_flag
        field = fix_phase(field)
    end

    data = _extract_quantity(field, quantity)
    N = length(data)
    x = range(0, 1, length=N)

    plot!(p, x, data; kwargs...)
    return p
end

function PhoXonic.plot_field!(
    p,
    solver::Solver{Dim2,W},
    eigenvector::AbstractVector;
    component::Symbol=:auto,
    quantity::Symbol=:real,
    fix_phase_flag::Bool=true,
    colormap::Symbol=:RdBu,
    kwargs...
) where {W}
    field = reconstruct_field(solver, eigenvector)

    if field isa Tuple
        comp_idx = _get_component_index(W, component)
        data_complex = field[comp_idx]
    else
        data_complex = field
    end

    if fix_phase_flag
        data_complex = fix_phase(data_complex)
    end

    data = _extract_quantity(data_complex, quantity)

    heatmap!(p, data'; c=colormap, aspect_ratio=:equal, kwargs...)
    return p
end

# -----------------------------------------------------------------------------
# plot_epsilon - 1D
# -----------------------------------------------------------------------------
"""
    plot_epsilon(solver::Solver{Dim1}; kwargs...)

Plot 1D permittivity distribution.
"""
function PhoXonic.plot_epsilon(
    solver::Solver{Dim1};
    title::String="Permittivity",
    xlabel::String="x/a",
    ylabel::String="ε",
    kwargs...
)
    eps = get_epsilon_field(solver)
    N = length(eps)
    x = range(0, 1, length=N)

    Plots.plot(x, eps;
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        fill=0,
        alpha=0.3,
        legend=false,
        kwargs...
    )
end

# -----------------------------------------------------------------------------
# plot_epsilon - 2D
# -----------------------------------------------------------------------------
"""
    plot_epsilon(solver::Solver{Dim2}; kwargs...)

Plot 2D permittivity distribution as a heatmap.
"""
function PhoXonic.plot_epsilon(
    solver::Solver{Dim2};
    colormap::Symbol=:grays,
    title::String="Permittivity",
    xlabel::String="x/a",
    ylabel::String="y/a",
    kwargs...
)
    eps = get_epsilon_field(solver)

    Plots.heatmap(eps';
        c=colormap,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        aspect_ratio=:equal,
        kwargs...
    )
end

# -----------------------------------------------------------------------------
# Helper functions for field visualization
# -----------------------------------------------------------------------------
function _extract_quantity(field::AbstractArray{<:Complex}, quantity::Symbol)
    if quantity == :real
        real.(field)
    elseif quantity == :imag
        imag.(field)
    elseif quantity == :abs
        abs.(field)
    elseif quantity == :phase
        angle.(field)
    else
        error("Unknown quantity: $quantity. Use :real, :imag, :abs, or :phase")
    end
end

function _get_component_index(::Type{<:WaveType}, component::Symbol)
    # Default: first component
    component == :auto && return 1
    component in (:x, :first, :1) && return 1
    component in (:y, :second, :2) && return 2
    component in (:z, :third, :3) && return 3
    error("Unknown component: $component")
end

end # module
