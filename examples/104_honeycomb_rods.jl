# Last-Modified: 2025-12-15T19:00:00+09:00
# Honeycomb Lattice Photonic Crystal Band Structure
# Equivalent to MPB's honey-rods.ctl
# Two dielectric rods per unit cell - both TM and TE polarization

using PhoXonic
using Plots
default(guidefontsize=14, tickfontsize=12, titlefontsize=14, left_margin=10Plots.mm, right_margin=10Plots.mm, top_margin=5Plots.mm, bottom_margin=10Plots.mm)

println("=== Honeycomb Lattice Photonic Crystal ===")
println("Dielectric rods (ε=12) in air, r=0.14a")
println("Two rods per unit cell at (1/6, 1/6) and (-1/6, -1/6)")
println("Equivalent to MPB honey-rods.ctl")
println()

# Triangular lattice (honeycomb is triangular with 2-atom basis)
a = 1.0
lat = hexagonal_lattice(a)

# Geometry: air background with two dielectric rods
# Rod positions in lattice coordinates: (1/6, 1/6) and (-1/6, -1/6)
# Convert to Cartesian: r = s1*a1 + s2*a2
a1 = lat.vectors[1]
a2 = lat.vectors[2]

# Position 1: (1/6, 1/6) in lattice coords
pos1 = (1/6) * a1 + (1/6) * a2
# Position 2: (-1/6, -1/6) in lattice coords
pos2 = (-1/6) * a1 + (-1/6) * a2

println("Lattice vectors:")
println("  a1 = ", a1)
println("  a2 = ", a2)
println("Rod positions (Cartesian):")
println("  Rod 1: ", pos1)
println("  Rod 2: ", pos2)

air = Dielectric(1.0)
rod = Dielectric(12.0)

geo = Geometry(
    lat, air, [(Circle(collect(pos1), 0.14), rod), (Circle(collect(pos2), 0.14), rod)]
)

# Create k-path: Γ → M → K → Γ
kpath = simple_kpath_hexagonal(; a=a, npoints=50)

# ============================================================================
# TM polarization
# ============================================================================
println("\nComputing TM bands...")
solver_tm = Solver(TMWave(), geo, (64, 64); cutoff=7)
bands_tm = compute_bands(solver_tm, kpath; bands=1:14, verbose=false)

println("\nTM Band Gaps:")
gaps_tm = find_all_gaps(bands_tm; threshold=0.005)
if isempty(gaps_tm)
    println("  No significant band gaps found.")
else
    for g in gaps_tm
        println(
            "  Gap between bands ",
            g.bands,
            ": ",
            round(g.gap_ratio*100; digits=1),
            "% gap-to-midgap",
        )
    end
end

# ============================================================================
# TE polarization
# ============================================================================
println("\nComputing TE bands...")
solver_te = Solver(TEWave(), geo, (64, 64); cutoff=7)
bands_te = compute_bands(solver_te, kpath; bands=1:14, verbose=false)

println("\nTE Band Gaps:")
gaps_te = find_all_gaps(bands_te; threshold=0.005)
if isempty(gaps_te)
    println("  No significant band gaps found.")
else
    for g in gaps_te
        println(
            "  Gap between bands ",
            g.bands,
            ": ",
            round(g.gap_ratio*100; digits=1),
            "% gap-to-midgap",
        )
    end
end

# ============================================================================
# Check for complete band gap (overlapping TE and TM)
# ============================================================================
println("\n=== Complete Band Gap Analysis ===")

# Get frequency ranges for each polarization
tm_freqs = bands_tm.frequencies
te_freqs = bands_te.frequencies

# Find overlapping gaps
for n in 1:13
    tm_gap = find_bandgap(bands_tm, n, n+1)
    te_gap = find_bandgap(bands_te, n, n+1)

    if tm_gap.gap > 0 && te_gap.gap > 0
        # Check for overlap
        overlap_min = max(tm_gap.max_lower, te_gap.max_lower)
        overlap_max = min(tm_gap.min_upper, te_gap.min_upper)
        if overlap_max > overlap_min
            midgap = (overlap_min + overlap_max) / 2
            gap_ratio = (overlap_max - overlap_min) / midgap * 100
            println("Complete gap between bands $n and $(n+1):")
            println(
                "  TM: $(round(tm_gap.max_lower, digits=3)) - $(round(tm_gap.min_upper, digits=3))",
            )
            println(
                "  TE: $(round(te_gap.max_lower, digits=3)) - $(round(te_gap.min_upper, digits=3))",
            )
            println(
                "  Overlap: $(round(overlap_min, digits=3)) - $(round(overlap_max, digits=3))",
            )
            println("  Complete gap ratio: $(round(gap_ratio, digits=1))%")
        end
    end
end

# ============================================================================
# Create plots
# ============================================================================
dists = bands_tm.distances
label_positions = [dists[i] for (i, _) in bands_tm.labels]
label_names = [l for (_, l) in bands_tm.labels]

ymax = max(maximum(bands_tm.frequencies), maximum(bands_te.frequencies)) * 1.05
ylims_common = (0, ymax)

# TM plot
p_tm = plot(;
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πc)",
    title="TM Bands (ε=12 rods)",
    legend=false,
    grid=true,
    size=(700, 500),
    ylims=ylims_common,
)
for b in 1:size(bands_tm.frequencies, 2)
    plot!(p_tm, dists, bands_tm.frequencies[:, b]; linewidth=2, color=:blue)
end
vline!(p_tm, label_positions; color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_tm, label_positions, label_names)

# Highlight TM band gap (bands 2-3)
gap_tm_2_3 = find_bandgap(bands_tm, 2, 3)
if gap_tm_2_3.gap > 0
    hspan!(p_tm, [gap_tm_2_3.max_lower, gap_tm_2_3.min_upper]; alpha=0.2, color=:blue, label="")
end

# TE plot
p_te = plot(;
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πc)",
    title="TE Bands (ε=12 rods)",
    legend=false,
    grid=true,
    size=(700, 500),
    ylims=ylims_common,
)
for b in 1:size(bands_te.frequencies, 2)
    plot!(p_te, dists, bands_te.frequencies[:, b]; linewidth=2, color=:red)
end
vline!(p_te, label_positions; color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_te, label_positions, label_names)

# Highlight TE band gap (bands 5-6)
gap_te_5_6 = find_bandgap(bands_te, 5, 6)
if gap_te_5_6.gap > 0
    hspan!(p_te, [gap_te_5_6.max_lower, gap_te_5_6.min_upper]; alpha=0.2, color=:red, label="")
end

# Combined plot (side by side)
p_combined = plot(
    p_tm,
    p_te;
    layout=(1, 2),
    size=(1200, 500),
    plot_title="Honeycomb Lattice Photonic Crystal (r=0.14a)",
)

# Save plots
savefig(p_tm, joinpath(@__DIR__, "104_honeycomb_tm_bands.png"))
savefig(p_te, joinpath(@__DIR__, "104_honeycomb_te_bands.png"))
savefig(p_combined, joinpath(@__DIR__, "104_honeycomb_bands.png"))

println("\nSaved: 104_honeycomb_tm_bands.png")
println("Saved: 104_honeycomb_te_bands.png")
println("Saved: 104_honeycomb_bands.png")

display(p_combined)
