#=
Makie.jl extension for PhoXonic.jl

This extension is automatically loaded when both PhoXonic and Makie are used.
Works with CairoMakie, GLMakie, or WGLMakie.
=#

module PhoXonicMakieExt

using PhoXonic
using Makie

import PhoXonic: plot_bands, plot_bands!, band_plot_data
import PhoXonic: plot_field, plot_field!, plot_epsilon, surface_field
import PhoXonic: reconstruct_field, get_epsilon_field, get_material_field, fix_phase
import PhoXonic: Dim1, Dim2, Dim3, Solver, WaveType, BandStructure

# =============================================================================
# plot_bands
# =============================================================================

"""
    plot_bands(bs::BandStructure; kwargs...)

Plot a band structure diagram using Makie.

# Keyword Arguments
- `color`: Line/marker color (default: :blue)
- `linewidth`: Line width (default: 2)
- `scatter`: Use scatter plot instead of lines (default: false)
- `markersize`: Marker size for scatter plot (default: 6)
- `title`: Plot title (default: "Band Structure")
- `ylabel`: Y-axis label (default: "Frequency")
- `xlabel`: X-axis label (default: "Wave vector")
- `show_gaps`: Highlight band gaps (default: false)
- `gap_color`: Color for gap highlighting (default: :yellow)
- `gap_alpha`: Transparency for gap highlighting (default: 0.2)
- `size`: Figure size (default: (600, 400))

# Returns
A Makie `Figure` object.
"""
function PhoXonic.plot_bands(
    bs::BandStructure;
    color=:blue,
    linewidth::Real=2,
    scatter::Bool=false,
    markersize::Real=6,
    title::String="Band Structure",
    ylabel::String="Frequency",
    xlabel::String="Wave vector",
    show_gaps::Bool=false,
    gap_color=:yellow,
    gap_alpha::Real=0.2,
    gap_threshold::Real=0.0,
    normalize::Real=1.0,
    size::Tuple{Int,Int}=(600, 400),
)
    data = band_plot_data(bs; normalize=normalize)

    fig = Figure(; size=size)
    ax = Axis(
        fig[1, 1];
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        xticks=(data.label_positions, data.label_names),
    )

    # Plot each band
    nbands = Base.size(data.frequencies, 2)
    for b in 1:nbands
        if scatter
            Makie.scatter!(
                ax,
                data.distances,
                data.frequencies[:, b];
                markersize=markersize,
                color=color,
            )
        else
            lines!(
                ax, data.distances, data.frequencies[:, b]; linewidth=linewidth, color=color
            )
        end
    end

    # Add high-symmetry point markers
    if !isempty(data.label_positions)
        vlines!(ax, data.label_positions; color=:gray, linestyle=:dash, alpha=0.5)
    end

    # Highlight band gaps
    if show_gaps
        gaps = find_all_gaps(bs; threshold=gap_threshold)
        for g in gaps
            hspan!(
                ax,
                g.max_lower * normalize,
                g.min_upper * normalize;
                color=(gap_color, gap_alpha),
            )
        end
    end

    return fig
end

"""
    plot_bands!(ax, bs::BandStructure; kwargs...)

Add band structure to an existing Makie Axis.
"""
function PhoXonic.plot_bands!(
    ax,
    bs::BandStructure;
    color=:blue,
    linewidth::Real=2,
    linestyle=:solid,
    scatter::Bool=false,
    markersize::Real=6,
    normalize::Real=1.0,
)
    data = band_plot_data(bs; normalize=normalize)

    nbands = Base.size(data.frequencies, 2)
    for b in 1:nbands
        if scatter
            Makie.scatter!(
                ax,
                data.distances,
                data.frequencies[:, b];
                markersize=markersize,
                color=color,
            )
        else
            lines!(
                ax,
                data.distances,
                data.frequencies[:, b];
                linewidth=linewidth,
                color=color,
                linestyle=linestyle,
            )
        end
    end

    return ax
end

# =============================================================================
# plot_field - 1D
# =============================================================================

"""
    plot_field(solver::Solver{Dim1}, eigenvector; kwargs...)

Plot a 1D field distribution using Makie.

# Keyword Arguments
- `quantity::Symbol = :real` - Quantity to plot: `:real`, `:imag`, `:abs`, `:phase`
- `fix_phase_flag::Bool = true` - Apply phase normalization
- `title::String = "Field"` - Plot title
- `xlabel::String = "x/a"` - X-axis label
- `ylabel::String = "Field"` - Y-axis label
- `size::Tuple{Int,Int} = (600, 400)` - Figure size

# Returns
A Makie `Figure` object.
"""
function PhoXonic.plot_field(
    solver::Solver{Dim1},
    eigenvector::AbstractVector;
    quantity::Symbol=:real,
    fix_phase_flag::Bool=true,
    title::String="Field",
    xlabel::String="x/a",
    ylabel::String="Field",
    linewidth::Real=2,
    color=:blue,
    size::Tuple{Int,Int}=(600, 400),
)
    field = reconstruct_field(solver, eigenvector)
    if fix_phase_flag
        field = fix_phase(field)
    end

    data = _extract_quantity(field, quantity)
    N = length(data)
    x = range(0, 1; length=N)

    fig = Figure(; size=size)
    ax = Axis(fig[1, 1]; xlabel=xlabel, ylabel=ylabel, title=title)
    lines!(ax, x, data; linewidth=linewidth, color=color)

    return fig
end

# =============================================================================
# plot_field - 2D
# =============================================================================

"""
    plot_field(solver::Solver{Dim2}, eigenvector; kwargs...)

Plot a 2D field distribution as a heatmap using Makie.

# Keyword Arguments
- `component::Symbol = :auto` - Component to plot for vector fields
- `quantity::Symbol = :real` - Quantity to plot: `:real`, `:imag`, `:abs`, `:phase`
- `fix_phase_flag::Bool = true` - Apply phase normalization
- `colormap = :RdBu` - Colormap
- `title::String = "Field"` - Plot title
- `xlabel::String = "x/a"` - X-axis label
- `ylabel::String = "y/a"` - Y-axis label
- `size::Tuple{Int,Int} = (650, 500)` - Figure size

# Returns
A Makie `Figure` object.
"""
function PhoXonic.plot_field(
    solver::Solver{Dim2,W},
    eigenvector::AbstractVector;
    component::Symbol=:auto,
    quantity::Symbol=:real,
    fix_phase_flag::Bool=true,
    colormap=:RdBu,
    title::String="Field",
    xlabel::String="x/a",
    ylabel::String="y/a",
    size::Tuple{Int,Int}=(650, 500),
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

    # Symmetric colorrange for real/imag
    if quantity in (:real, :imag)
        maxval = maximum(abs.(data))
        colorrange = (-maxval, maxval)
    else
        colorrange = extrema(data)
    end

    fig = Figure(; size=size)
    ax = Axis(fig[1, 1]; xlabel=xlabel, ylabel=ylabel, title=title, aspect=DataAspect())
    hm = heatmap!(ax, data'; colormap=colormap, colorrange=colorrange)
    Colorbar(fig[1, 2], hm)

    return fig
end

# =============================================================================
# plot_field! - Add to existing axis
# =============================================================================

"""
    plot_field!(ax, solver::Solver{Dim1}, eigenvector; kwargs...)

Add a 1D field to an existing Makie Axis.
"""
function PhoXonic.plot_field!(
    ax,
    solver::Solver{Dim1},
    eigenvector::AbstractVector;
    quantity::Symbol=:real,
    fix_phase_flag::Bool=true,
    linewidth::Real=2,
    color=:blue,
    label=nothing,
)
    field = reconstruct_field(solver, eigenvector)
    if fix_phase_flag
        field = fix_phase(field)
    end

    data = _extract_quantity(field, quantity)
    N = length(data)
    x = range(0, 1; length=N)

    lines!(ax, x, data; linewidth=linewidth, color=color, label=label)
    return ax
end

"""
    plot_field!(ax, solver::Solver{Dim2}, eigenvector; kwargs...)

Add a 2D field heatmap to an existing Makie Axis.
"""
function PhoXonic.plot_field!(
    ax,
    solver::Solver{Dim2,W},
    eigenvector::AbstractVector;
    component::Symbol=:auto,
    quantity::Symbol=:real,
    fix_phase_flag::Bool=true,
    colormap=:RdBu,
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

    if quantity in (:real, :imag)
        maxval = maximum(abs.(data))
        colorrange = (-maxval, maxval)
    else
        colorrange = extrema(data)
    end

    heatmap!(ax, data'; colormap=colormap, colorrange=colorrange)
    return ax
end

# =============================================================================
# plot_epsilon
# =============================================================================

"""
    plot_epsilon(solver::Solver{Dim1}; kwargs...)

Plot 1D permittivity distribution using Makie.
"""
function PhoXonic.plot_epsilon(
    solver::Solver{Dim1};
    title::String="Permittivity",
    xlabel::String="x/a",
    ylabel::String="ε",
    size::Tuple{Int,Int}=(600, 400),
)
    eps = get_epsilon_field(solver)
    N = length(eps)
    x = range(0, 1; length=N)

    fig = Figure(; size=size)
    ax = Axis(fig[1, 1]; xlabel=xlabel, ylabel=ylabel, title=title)
    band!(ax, x, zeros(N), eps; color=(:blue, 0.3))
    lines!(ax, x, eps; color=:blue, linewidth=2)

    return fig
end

"""
    plot_epsilon(solver::Solver{Dim2}; kwargs...)

Plot 2D permittivity distribution as a heatmap using Makie.
"""
function PhoXonic.plot_epsilon(
    solver::Solver{Dim2};
    colormap=:grays,
    title::String="Permittivity",
    xlabel::String="x/a",
    ylabel::String="y/a",
    size::Tuple{Int,Int}=(650, 500),
)
    eps = get_epsilon_field(solver)

    fig = Figure(; size=size)
    ax = Axis(fig[1, 1]; xlabel=xlabel, ylabel=ylabel, title=title, aspect=DataAspect())
    hm = heatmap!(ax, eps'; colormap=colormap)
    Colorbar(fig[1, 2], hm)

    return fig
end

# =============================================================================
# surface_field - 2D field as 3D surface (Part 2)
# =============================================================================

"""
    surface_field(solver::Solver{Dim2}, eigenvector; kwargs...)

Plot a 2D field as a 3D surface (amplitude as height).

This is a Makie-specific function that provides enhanced visualization
by displaying field amplitude as the z-coordinate.

# Keyword Arguments
- `component::Symbol = :auto` - Component to plot for vector fields
- `quantity::Symbol = :real` - Quantity to plot: `:real`, `:imag`, `:abs`
- `fix_phase_flag::Bool = true` - Apply phase normalization
- `colormap = :RdBu` - Colormap
- `title::String = "Field"` - Plot title
- `size::Tuple{Int,Int} = (700, 600)` - Figure size

# Returns
A Makie `Figure` object with a 3D surface plot.

# Example
```julia
using GLMakie  # or CairoMakie
using PhoXonic

solver = Solver(TMWave(), geo, (64, 64); cutoff=7)
freqs, vecs = solve_at_k_with_vectors(solver, [0.0, 0.0], DenseMethod())
fig = surface_field(solver, vecs[:, 1])
```
"""
function PhoXonic.surface_field(
    solver::Solver{Dim2,W},
    eigenvector::AbstractVector;
    component::Symbol=:auto,
    quantity::Symbol=:real,
    fix_phase_flag::Bool=true,
    colormap=:RdBu,
    title::String="Field",
    xlabel::String="x/a",
    ylabel::String="y/a",
    zlabel::String="Field",
    size::Tuple{Int,Int}=(700, 600),
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

    Nx, Ny = Base.size(data)
    x = range(0, 1; length=Nx)
    y = range(0, 1; length=Ny)

    fig = Figure(; size=size)
    ax = Axis3(fig[1, 1]; xlabel=xlabel, ylabel=ylabel, zlabel=zlabel, title=title)
    surface!(ax, x, y, data; colormap=colormap)

    return fig
end

# =============================================================================
# plot_field - 3D cross-section (Part 2)
# =============================================================================

"""
    plot_field(solver::Solver{Dim3}, eigenvector; plane=:xy, slice=0.5, ...)

Plot a cross-section of a 3D field as a 2D heatmap.

# Keyword Arguments
- `plane::Symbol = :xy` - Cross-section plane: `:xy`, `:xz`, or `:yz`
- `slice::Real = 0.5` - Position of the slice (0 to 1, normalized)
- `component::Symbol = :auto` - Component for vector fields
- `quantity::Symbol = :real` - Quantity to plot
- `fix_phase_flag::Bool = true` - Apply phase normalization
- `colormap = :RdBu` - Colormap
- `size::Tuple{Int,Int} = (650, 500)` - Figure size

# Returns
A Makie `Figure` object.
"""
function PhoXonic.plot_field(
    solver::Solver{Dim3,W},
    eigenvector::AbstractVector;
    plane::Symbol=:xy,
    slice::Real=0.5,
    component::Symbol=:auto,
    quantity::Symbol=:real,
    fix_phase_flag::Bool=true,
    colormap=:RdBu,
    title::String="Field",
    size::Tuple{Int,Int}=(650, 500),
) where {W}
    field = reconstruct_field(solver, eigenvector)

    # Handle tuple (vector field) or array (scalar field)
    if field isa Tuple
        comp_idx = _get_component_index(W, component)
        data_complex = field[comp_idx]
    else
        data_complex = field
    end

    if fix_phase_flag
        data_complex = fix_phase(data_complex)
    end

    data_3d = _extract_quantity(data_complex, quantity)

    # Extract cross-section
    Nx, Ny, Nz = Base.size(data_3d)
    if plane == :xy
        iz = clamp(round(Int, slice * Nz), 1, Nz)
        data_2d = data_3d[:, :, iz]
        xlabel, ylabel = "x/a", "y/a"
        slice_label = "z = $(round(slice, digits=2))"
    elseif plane == :xz
        iy = clamp(round(Int, slice * Ny), 1, Ny)
        data_2d = data_3d[:, iy, :]
        xlabel, ylabel = "x/a", "z/a"
        slice_label = "y = $(round(slice, digits=2))"
    elseif plane == :yz
        ix = clamp(round(Int, slice * Nx), 1, Nx)
        data_2d = data_3d[ix, :, :]
        xlabel, ylabel = "y/a", "z/a"
        slice_label = "x = $(round(slice, digits=2))"
    else
        error("plane must be :xy, :xz, or :yz")
    end

    # Symmetric colorrange for real/imag
    if quantity in (:real, :imag)
        maxval = maximum(abs.(data_2d))
        colorrange = (-maxval, maxval)
    else
        colorrange = extrema(data_2d)
    end

    full_title = isempty(title) ? slice_label : "$title ($slice_label)"

    fig = Figure(; size=size)
    ax = Axis(
        fig[1, 1]; xlabel=xlabel, ylabel=ylabel, title=full_title, aspect=DataAspect()
    )
    hm = heatmap!(ax, data_2d'; colormap=colormap, colorrange=colorrange)
    Colorbar(fig[1, 2], hm)

    return fig
end

# =============================================================================
# Helper functions
# =============================================================================

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
    component == :auto && return 1
    component in (:x, :first, :1) && return 1
    component in (:y, :second, :2) && return 2
    component in (:z, :third, :3) && return 3
    return error("Unknown component: $component")
end

end # module
