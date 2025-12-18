# Last-Modified: 2025-12-15T19:00:00+09:00
# Plot unit cell structures for all examples
# Visualize the geometry of photonic crystals using shape drawing

using PhoXonic
using Plots

println("Generating structure plots...")

# ============================================================================
# Helper function to draw a circle
# ============================================================================
function draw_circle!(p, cx, cy, r; color=:blue, fillalpha=1.0, n=100)
    θ = range(0, 2π, length=n)
    xs = cx .+ r .* cos.(θ)
    ys = cy .+ r .* sin.(θ)
    plot!(p, xs, ys, seriestype=:shape, fillcolor=color, fillalpha=fillalpha,
          linecolor=color, linewidth=0.5, label="")
end

# ============================================================================
# Helper function to plot unit cell with shapes
# ============================================================================
function plot_unit_cell_shapes(lat, bg_eps, inclusions, title_str;
                               xlim=(-0.6, 1.6), ylim=(-0.6, 1.6))
    a1 = lat.vectors[1]
    a2 = lat.vectors[2]

    # Determine colors based on background
    if bg_eps > 5
        bg_color = :goldenrod
        inc_color = :navy
    else
        bg_color = :navy
        inc_color = :goldenrod
    end

    # Create plot with background
    p = plot(
        aspect_ratio=:equal,
        xlabel="x/a",
        ylabel="y/a",
        title=title_str,
        xlim=xlim,
        ylim=ylim,
        size=(550, 500),
        legend=false,
        grid=false,
        background_color_inside=bg_color
    )

    # Draw inclusions at all visible unit cells
    for di in -1:2
        for dj in -1:2
            offset = di * a1 + dj * a2
            for (center, radius, eps) in inclusions
                cx = center[1] + offset[1]
                cy = center[2] + offset[2]
                col = eps > bg_eps ? :goldenrod : :navy
                draw_circle!(p, cx, cy, radius; color=col)
            end
        end
    end

    # Draw unit cell boundary
    corners = [
        [0.0, 0.0],
        collect(a1),
        collect(a1 + a2),
        collect(a2),
        [0.0, 0.0]
    ]
    xs = [c[1] for c in corners]
    ys = [c[2] for c in corners]
    plot!(p, xs, ys, color=:white, linewidth=3, linestyle=:solid, label="")

    return p
end

# ============================================================================
# 1. Triangular lattice rods (tri-rods)
# ============================================================================
println("  1. Triangular rods...")
lat_tri = hexagonal_lattice(1.0)
inclusions_tri_rods = [([0.0, 0.0], 0.2, 12.0)]
p1 = plot_unit_cell_shapes(lat_tri, 1.0, inclusions_tri_rods,
                           "Triangular Rods\n(ε=12, r=0.2a)")

# ============================================================================
# 2. Square lattice rods (sq-rods)
# ============================================================================
println("  2. Square rods...")
lat_sq = square_lattice(1.0)
inclusions_sq_rods = [([0.0, 0.0], 0.2, 11.56)]
p2 = plot_unit_cell_shapes(lat_sq, 1.0, inclusions_sq_rods,
                           "Square Rods\n(ε=11.56, r=0.2a)";
                           xlim=(-0.3, 1.3), ylim=(-0.3, 1.3))

# ============================================================================
# 3. Honeycomb lattice rods (honey-rods)
# ============================================================================
println("  3. Honeycomb rods...")
lat_hex = hexagonal_lattice(1.0)
a1 = lat_hex.vectors[1]
a2 = lat_hex.vectors[2]
pos1 = (1/6) * a1 + (1/6) * a2
pos2 = (-1/6) * a1 + (-1/6) * a2
inclusions_honey = [
    (collect(pos1), 0.14, 12.0),
    (collect(pos2), 0.14, 12.0)
]
p3 = plot_unit_cell_shapes(lat_hex, 1.0, inclusions_honey,
                           "Honeycomb Rods\n(ε=12, r=0.14a)")

# ============================================================================
# 4. Triangular lattice holes (tri-holes)
# ============================================================================
println("  4. Triangular holes...")
inclusions_tri_holes = [([0.0, 0.0], 0.45, 1.0)]
p4 = plot_unit_cell_shapes(lat_tri, 12.0, inclusions_tri_holes,
                           "Triangular Holes\n(ε=12 bg, r=0.45a)")

# ============================================================================
# Save individual plots
# ============================================================================
savefig(p1, joinpath(@__DIR__, "801_structure_tri_rods.png"))
savefig(p2, joinpath(@__DIR__, "801_structure_sq_rods.png"))
savefig(p3, joinpath(@__DIR__, "801_structure_honeycomb.png"))
savefig(p4, joinpath(@__DIR__, "801_structure_tri_holes.png"))

println("\nSaved individual structure plots.")

# ============================================================================
# Combined 2x2 plot
# ============================================================================
p_all = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 900))
savefig(p_all, joinpath(@__DIR__, "801_structures_all.png"))

println("Saved: 801_structures_all.png")

display(p_all)
