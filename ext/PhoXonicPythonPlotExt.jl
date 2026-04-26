# SPDX-License-Identifier: MIT
# Copyright (C) 2026 Hiroharu Sugawara
# Part of PhoXonic.jl - PythonPlot extension
module PhoXonicPythonPlotExt

using PhoXonic
import PhoXonic: BandStructure, band_plot_data, savefig_publication
import PythonPlot

const CM_PER_INCH = 2.54

"""
Plot a single BandStructure on the given matplotlib Axes.
"""
function _plot_bs_on_ax!(ax, bs::BandStructure;
                         color="black", linewidth=1.5, linestyle="-",
                         show_gaps::Bool=false,
                         normalize::Real=1.0,
                         title::AbstractString="")
    data = band_plot_data(bs; normalize)
    nbands = size(data.frequencies, 2)

    for i in 1:nbands
        ax.plot(data.distances, data.frequencies[:, i];
                color, linewidth, linestyle)
    end

    # Vertical lines at high-symmetry points
    for pos in data.label_positions
        ax.axvline(pos; color="gray", linewidth=0.5)
    end

    # Optional gap shading between consecutive bands
    if show_gaps && nbands >= 2
        for b in 1:nbands-1
            max_lower = maximum(data.frequencies[:, b])
            min_upper = minimum(data.frequencies[:, b+1])
            if min_upper > max_lower
                ax.axhspan(max_lower, min_upper; alpha=0.2, color="yellow")
            end
        end
    end

    ax.set_xticks(data.label_positions)
    ax.set_xticklabels(data.label_names)
    ax.set_ylabel("Frequency")
    if !isempty(title)
        ax.set_title(title)
    end
    ax.set_xlim(data.distances[1], data.distances[end])
end

"""
Compute figure size (inches) and axes positions from axis dimensions and layout.
"""
function _layout_axes(axis_width_cm, axis_height_cm, n;
                      margin_left_cm=1.5, margin_right_cm=0.3,
                      margin_bottom_cm=1.0, margin_top_cm=0.8,
                      hgap_cm=1.8, vgap_cm=1.5,
                      nrows=1, ncols=1)
    widths = axis_width_cm isa AbstractVector ? axis_width_cm : fill(axis_width_cm, ncols)
    heights = axis_height_cm isa AbstractVector ? axis_height_cm : fill(axis_height_cm, nrows)

    fig_w_cm = margin_left_cm + sum(widths) + hgap_cm * (ncols - 1) + margin_right_cm
    fig_h_cm = margin_bottom_cm + sum(heights) + vgap_cm * (nrows - 1) + margin_top_cm

    fig_w = fig_w_cm / CM_PER_INCH
    fig_h = fig_h_cm / CM_PER_INCH

    positions = Vector{NTuple{4,Float64}}()
    for row in 1:nrows
        for col in 1:ncols
            left = (margin_left_cm + sum(widths[1:col-1]) + hgap_cm * (col - 1)) / fig_w_cm
            bottom = (margin_bottom_cm + sum(heights[row+1:end]) + vgap_cm * (nrows - row)) / fig_h_cm
            w = widths[col] / fig_w_cm
            h = heights[row] / fig_h_cm
            push!(positions, (left, bottom, w, h))
        end
    end

    return fig_w, fig_h, positions
end

# ── Single BandStructure ──

function PhoXonic.savefig_publication(
    bs::BandStructure, path::AbstractString;
    axis_width_cm=8.0, axis_height_cm=6.0,
    ylims=nothing,
    kwargs...)

    fig_w, fig_h, positions = _layout_axes(axis_width_cm, axis_height_cm, 1)

    fig = PythonPlot.figure(figsize=(fig_w, fig_h))
    ax = fig.add_axes(collect(positions[1]))

    _plot_bs_on_ax!(ax, bs; kwargs...)
    if ylims !== nothing
        ax.set_ylim(ylims...)
    end

    fig.savefig(path)
    PythonPlot.close(fig)
    return path
end

# Default colors/styles for overlay
const _OVERLAY_STYLES = [
    (color="black",  linewidth=1.5, linestyle="-"),
    (color="blue",   linewidth=1.2, linestyle="--"),
    (color="red",    linewidth=1.2, linestyle="-."),
    (color="green",  linewidth=1.2, linestyle=":"),
    (color="purple", linewidth=1.2, linestyle="--"),
]

# ── Vector{BandStructure} ──

function PhoXonic.savefig_publication(
    bss::AbstractVector{<:BandStructure}, path::AbstractString;
    axis_width_cm=8.0, axis_height_cm=6.0,
    layout=(1, length(bss)),
    overlay::Bool=false,
    ylims=nothing,
    title::AbstractString="",
    kwargs...)

    if overlay
        fig_w, fig_h, positions = _layout_axes(axis_width_cm, axis_height_cm, 1)
        fig = PythonPlot.figure(figsize=(fig_w, fig_h))
        ax = fig.add_axes(collect(positions[1]))

        for (i, bs) in enumerate(bss)
            style = _OVERLAY_STYLES[mod1(i, length(_OVERLAY_STYLES))]
            _plot_bs_on_ax!(ax, bs;
                color=style.color, linewidth=style.linewidth,
                linestyle=style.linestyle)
        end

        if ylims !== nothing
            ax.set_ylim(ylims...)
        end
        if !isempty(title)
            ax.set_title(title)
        end
    else
        nrows, ncols = layout
        n = length(bss)

        fig_w, fig_h, positions = _layout_axes(
            axis_width_cm, axis_height_cm, n;
            nrows, ncols)

        fig = PythonPlot.figure(figsize=(fig_w, fig_h))

        for (i, bs) in enumerate(bss)
            i > length(positions) && break
            ax = fig.add_axes(collect(positions[i]))
            _plot_bs_on_ax!(ax, bs; kwargs...)
            if ylims !== nothing
                ax.set_ylim(ylims...)
            end
        end
    end

    fig.savefig(path)
    PythonPlot.close(fig)
    return path
end

end # module
