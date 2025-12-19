# Last-Modified: 2025-12-15T19:00:00+09:00
# Triangular Lattice Air Holes in Dielectric
# Equivalent to MPB's tri-holes.ctl
# Air holes in dielectric slab - both TM and TE polarization

using PhoXonic
using Plots

println("=== Triangular Lattice Air Holes ===")
println("Air holes (ε=1) in dielectric (ε=12), r=0.45a")
println("Equivalent to MPB tri-holes.ctl")
println()

# Triangular lattice
a = 1.0
lat = hexagonal_lattice(a)

# Geometry: dielectric background with air hole at center
dielectric = Dielectric(12.0)
air = Dielectric(1.0)

# Large hole radius r=0.45a
geo = Geometry(lat, dielectric, [(Circle([0.0, 0.0], 0.45), air)])

# Create k-path: Γ → M → K → Γ
kpath = simple_kpath_hexagonal(a=a, npoints=50)

# ============================================================================
# TM polarization
# ============================================================================
println("Computing TM bands...")
solver_tm = Solver(TMWave(), geo, (64, 64); cutoff=7)
bands_tm = compute_bands(solver_tm, kpath; bands=1:8, verbose=false)

println("\nTM Band Gaps:")
gaps_tm = find_all_gaps(bands_tm; threshold=0.005)
if isempty(gaps_tm)
    println("  No significant band gaps found.")
else
    for g in gaps_tm
        println("  Gap between bands ", g.bands, ": ",
                round(g.gap_ratio*100, digits=1), "% gap-to-midgap")
    end
end

# ============================================================================
# TE polarization
# ============================================================================
println("\nComputing TE bands...")
solver_te = Solver(TEWave(), geo, (64, 64); cutoff=7)
bands_te = compute_bands(solver_te, kpath; bands=1:8, verbose=false)

println("\nTE Band Gaps:")
gaps_te = find_all_gaps(bands_te; threshold=0.005)
if isempty(gaps_te)
    println("  No significant band gaps found.")
else
    for g in gaps_te
        println("  Gap between bands ", g.bands, ": ",
                round(g.gap_ratio*100, digits=1), "% gap-to-midgap")
    end
end

# ============================================================================
# Check for complete band gap (overlapping TE and TM)
# ============================================================================
println("\n=== Complete Band Gap Analysis ===")

found_complete = false
for n in 1:7
    tm_gap = find_bandgap(bands_tm, n, n+1)
    te_gap = find_bandgap(bands_te, n, n+1)

    if tm_gap.gap > 0 && te_gap.gap > 0
        # Check for overlap
        overlap_min = max(tm_gap.max_lower, te_gap.max_lower)
        overlap_max = min(tm_gap.min_upper, te_gap.min_upper)
        if overlap_max > overlap_min
            global found_complete = true
            midgap = (overlap_min + overlap_max) / 2
            gap_ratio = (overlap_max - overlap_min) / midgap * 100
            println("Complete gap between bands $n and $(n+1):")
            println("  TM: $(round(tm_gap.max_lower, digits=4)) - $(round(tm_gap.min_upper, digits=4))")
            println("  TE: $(round(te_gap.max_lower, digits=4)) - $(round(te_gap.min_upper, digits=4))")
            println("  Overlap: $(round(overlap_min, digits=4)) - $(round(overlap_max, digits=4))")
            println("  Complete gap ratio: $(round(gap_ratio, digits=1))%")
        end
    end
end

if !found_complete
    println("No complete (TM+TE) band gap found in bands 1-8.")
    println("Note: TE gaps are expected for holes in dielectric.")
end

# ============================================================================
# Create plots
# ============================================================================
dists = bands_tm.distances
label_positions = [dists[i] for (i, _) in bands_tm.labels]
label_names = [l for (_, l) in bands_tm.labels]

# Common y-axis range
ymax = max(maximum(bands_tm.frequencies), maximum(bands_te.frequencies)) * 1.05
ylims_common = (0, ymax)

# TM plot
p_tm = plot(
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πc)",
    title="TM Bands (r=0.45a)",
    legend=false,
    grid=true,
    size=(700, 450),
    ylims=ylims_common
)
for b in 1:size(bands_tm.frequencies, 2)
    plot!(p_tm, dists, bands_tm.frequencies[:, b], linewidth=2, color=:blue)
end
vline!(p_tm, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_tm, label_positions, label_names)

# TE plot
p_te = plot(
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πc)",
    title="TE Bands (r=0.45a)",
    legend=false,
    grid=true,
    size=(700, 450),
    ylims=ylims_common
)
for b in 1:size(bands_te.frequencies, 2)
    plot!(p_te, dists, bands_te.frequencies[:, b], linewidth=2, color=:red)
end
vline!(p_te, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_te, label_positions, label_names)

# Highlight TE band gap if exists
te_gap_1_2 = find_bandgap(bands_te, 1, 2)
if te_gap_1_2.gap > 0
    hspan!(p_te, [te_gap_1_2.max_lower, te_gap_1_2.min_upper],
           alpha=0.2, color=:red, label="")
end

# Combined plot (side by side)
p_combined = plot(p_tm, p_te, layout=(1, 2), size=(1200, 450),
    plot_title="Triangular Lattice Air Holes in ε=12")

# Save plots
savefig(p_tm, joinpath(@__DIR__, "111_triholes_tm_bands.png"))
savefig(p_te, joinpath(@__DIR__, "111_triholes_te_bands.png"))
savefig(p_combined, joinpath(@__DIR__, "111_triholes_bands_combined.png"))

println("\nSaved: 111_triholes_tm_bands.png")
println("Saved: 111_triholes_te_bands.png")
println("Saved: 111_triholes_bands_combined.png")

display(p_combined)
