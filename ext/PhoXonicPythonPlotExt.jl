# SPDX-License-Identifier: MIT
# Copyright (C) 2026 Hiroharu Sugawara
# Part of PhoXonic.jl - PythonPlot extension
module PhoXonicPythonPlotExt

using PhoXonic
import PhoXonic: BandStructure, band_plot_data, savefig_publication
import PhoXonic: plot_on_axis!, figure_publication
using PythonPlot: PythonPlot

const MM_PER_INCH = 25.4

"""
Plot a single BandStructure on the given matplotlib Axes.
"""
function _plot_bs_on_ax!(
    ax,
    bs::BandStructure;
    color="black",
    linewidth=1.5,
    linestyle="-",
    show_gaps::Bool=false,
    normalize::Real=1.0,
    ylabel::AbstractString="Frequency",
    title::AbstractString="",
)
    data = band_plot_data(bs; normalize)
    nbands = size(data.frequencies, 2)

    for i in 1:nbands
        ax.plot(data.distances, data.frequencies[:, i]; color, linewidth, linestyle)
    end

    # Vertical lines at high-symmetry points
    for pos in data.label_positions
        ax.axvline(pos; color="gray", linewidth=0.5)
    end

    # Optional gap shading between consecutive bands
    if show_gaps && nbands >= 2
        for b in 1:(nbands - 1)
            max_lower = maximum(data.frequencies[:, b])
            min_upper = minimum(data.frequencies[:, b + 1])
            if min_upper > max_lower
                ax.axhspan(max_lower, min_upper; alpha=0.2, color="yellow")
            end
        end
    end

    ax.set_xticks(data.label_positions)
    ax.set_xticklabels(data.label_names)
    ax.set_ylabel(ylabel)
    if !isempty(title)
        ax.set_title(title)
    end
    ax.set_xlim(data.distances[1], data.distances[end])
end

# ── L1: plot_on_axis! ──

"""
    plot_on_axis!(ax, bs::BandStructure; kwargs...) -> ax

Draw a [`BandStructure`](@ref) onto `ax` and return `ax`. Bands are drawn as
lines with vertical guides and tick labels at the high-symmetry points.
Requires `using PythonPlot`; not exported (call as `PhoXonic.plot_on_axis!`).

# Arguments

- `ax`: a matplotlib axis (PythonCall `Py`) to draw onto
- `bs::BandStructure`: the band structure to plot

# Keywords

- `color="black"`: line color
- `linewidth=1.5`: line width
- `linestyle="-"`: line style
- `show_gaps::Bool=false`: shade band gaps between consecutive bands with `axhspan`
- `normalize::Real=1.0`: frequency normalization factor (forwarded to `band_plot_data`)
- `xlabel::AbstractString=""`: x axis label; empty leaves the axis labelled only by
    the high-symmetry-point xticks
- `ylabel::AbstractString="Frequency"`: y axis label
- `title::AbstractString=""`: title; empty leaves the axes untitled
"""
function plot_on_axis!(ax, bs::BandStructure; xlabel::AbstractString="", kwargs...)
    _plot_bs_on_ax!(ax, bs; kwargs...)
    isempty(xlabel) || ax.set_xlabel(xlabel)
    return ax
end

"""
Compute figure size (inches) and axes positions from axis dimensions and layout.
"""
function _layout_axes(
    axis_width_mm,
    axis_height_mm,
    n;
    margin_left_mm=15.0,
    margin_right_mm=3.0,
    margin_bottom_mm=10.0,
    margin_top_mm=8.0,
    hgap_mm=18.0,
    vgap_mm=15.0,
    nrows=1,
    ncols=1,
)
    widths = axis_width_mm isa AbstractVector ? axis_width_mm : fill(axis_width_mm, ncols)
    heights =
        axis_height_mm isa AbstractVector ? axis_height_mm : fill(axis_height_mm, nrows)

    fig_w_mm = margin_left_mm + sum(widths) + hgap_mm * (ncols - 1) + margin_right_mm
    fig_h_mm = margin_bottom_mm + sum(heights) + vgap_mm * (nrows - 1) + margin_top_mm

    fig_w = fig_w_mm / MM_PER_INCH
    fig_h = fig_h_mm / MM_PER_INCH

    positions = Vector{NTuple{4,Float64}}()
    for row in 1:nrows
        for col in 1:ncols
            left =
                (margin_left_mm + sum(widths[1:(col - 1)]) + hgap_mm * (col - 1)) / fig_w_mm
            bottom =
                (margin_bottom_mm + sum(heights[(row + 1):end]) + vgap_mm * (nrows - row)) /
                fig_h_mm
            w = widths[col] / fig_w_mm
            h = heights[row] / fig_h_mm
            push!(positions, (left, bottom, w, h))
        end
    end

    return fig_w, fig_h, positions
end

# ── L2: figure_publication (single BandStructure only) ──

"""
    figure_publication(bs; axis_width_mm=80.0, axis_height_mm=60.0, ylims=nothing, kwargs...) -> (fig, ax)

Create a publication matplotlib figure and a single axis for `bs`, draw it via
[`plot_on_axis!`](@ref), and return `(fig, ax)` so you can tweak it before
saving. The caller owns the figure and must call `PythonPlot.close(fig)`. Single
`BandStructure` only; for a `Vector` use [`savefig_publication`](@ref) or compose
with [`plot_on_axis!`](@ref) and your own `subplots()`. Requires `using
PythonPlot`; not exported (call as `PhoXonic.figure_publication`).

# Arguments

- `bs::BandStructure`: the band structure to plot

# Keywords

- `axis_width_mm=80.0`: plotting area width in mm
- `axis_height_mm=60.0`: plotting area height in mm
- `ylims=nothing`: pass a tuple to override the y axis limits
- `kwargs...`: forwarded to [`plot_on_axis!`](@ref)
"""
function figure_publication(
    bs::BandStructure; axis_width_mm=80.0, axis_height_mm=60.0, ylims=nothing, kwargs...
)
    fig_w, fig_h, positions = _layout_axes(axis_width_mm, axis_height_mm, 1)
    fig = PythonPlot.figure(; figsize=(fig_w, fig_h))
    try
        ax = fig.add_axes(collect(positions[1]))
        plot_on_axis!(ax, bs; kwargs...)
        if ylims !== nothing
            ax.set_ylim(ylims...)
        end
        return fig, ax
    catch
        PythonPlot.close(fig)
        rethrow()
    end
end

# ── L3: savefig_publication (single BandStructure → L2 delegation) ──

function PhoXonic.savefig_publication(bs::BandStructure, path::AbstractString; kwargs...)
    fig, _ = figure_publication(bs; kwargs...)
    try
        fig.savefig(path)
    finally
        PythonPlot.close(fig)
    end
    return path
end

# Default colors/styles for overlay
const _OVERLAY_STYLES = [
    (color="black", linewidth=1.5, linestyle="-"),
    (color="blue", linewidth=1.2, linestyle="--"),
    (color="red", linewidth=1.2, linestyle="-."),
    (color="green", linewidth=1.2, linestyle=":"),
    (color="purple", linewidth=1.2, linestyle="--"),
]

# ── L3: savefig_publication (Vector{BandStructure} → L1 direct loop, no L2) ──

function PhoXonic.savefig_publication(
    bss::AbstractVector{<:BandStructure},
    path::AbstractString;
    axis_width_mm=80.0,
    axis_height_mm=60.0,
    layout=(1, length(bss)),
    overlay::Bool=false,
    ylims=nothing,
    title::AbstractString="",
    kwargs...,
)
    if overlay
        fig_w, fig_h, positions = _layout_axes(axis_width_mm, axis_height_mm, 1)
    else
        nrows, ncols = layout
        fig_w, fig_h, positions = _layout_axes(
            axis_width_mm, axis_height_mm, length(bss); nrows, ncols
        )
    end

    fig = PythonPlot.figure(; figsize=(fig_w, fig_h))
    try
        if overlay
            ax = fig.add_axes(collect(positions[1]))
            for (i, bs) in enumerate(bss)
                style = _OVERLAY_STYLES[mod1(i, length(_OVERLAY_STYLES))]
                plot_on_axis!(ax, bs; style..., kwargs...)
            end
            if ylims !== nothing
                ax.set_ylim(ylims...)
            end
            if !isempty(title)
                ax.set_title(title)
            end
        else
            for (i, bs) in enumerate(bss)
                i > length(positions) && break
                ax = fig.add_axes(collect(positions[i]))
                plot_on_axis!(ax, bs; kwargs...)
                if ylims !== nothing
                    ax.set_ylim(ylims...)
                end
            end
        end

        fig.savefig(path)
    finally
        PythonPlot.close(fig)
    end
    return path
end

end # module
