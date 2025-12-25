# Last-Modified: 2025-12-15T19:00:00+09:00
# Joannopoulos Book Chapter 5, Figure 10 Reproduction
#
# Reference: "Photonic Crystals: Molding the Flow of Light" (2nd edition)
# ISBN: 978-0-691-12456-8, http://ab-initio.mit.edu/book/
#            Chapter 5: Two-Dimensional Photonic Crystals
#            Figure 10 (page 76)
#
# This example reproduces the photonic band structure for a triangular array
# of air columns drilled in a dielectric substrate.
#
# Parameters from the book:
#   - Triangular (hexagonal) lattice
#   - Substrate: ε = 13
#   - Air holes: ε = 1
#   - Blue bands = TM modes
#   - Red bands = TE modes
#   - Complete photonic band gap (overlapping TM and TE gaps)

using PhoXonic
using Printf
using Plots
default(guidefontsize=14, tickfontsize=12, titlefontsize=14, left_margin=10Plots.mm, right_margin=10Plots.mm, top_margin=5Plots.mm, bottom_margin=10Plots.mm)

println("=" ^ 70)
println("Joannopoulos Book - Chapter 5, Figure 10")
println("2D Photonic Band Structure: Triangular Array of Air Holes")
println("=" ^ 70)
println()

# ============================================================================
# Parameters (matching the book)
# ============================================================================

a = 1.0  # Lattice constant (normalized)

# Materials (from book caption: ε = 13 substrate)
ε_substrate = 13.0
ε_air = 1.0

substrate = Dielectric(ε_substrate)
air = Dielectric(ε_air)

# Air hole radius (typical value for complete gap)
# The book doesn't specify exact r/a, but r/a ≈ 0.45-0.48 gives complete gap
r = 0.45 * a

println("Parameters (from Joannopoulos Ch.5, Fig.10):")
println("  Lattice: 2D triangular (hexagonal), a = $a")
println("  Structure: Air holes in dielectric substrate")
@printf("  Substrate: ε = %.1f\n", ε_substrate)
@printf("  Air holes: ε = %.1f\n", ε_air)
@printf("  Hole radius: r/a = %.2f\n", r/a)
println()

# ============================================================================
# Geometry: Air holes in dielectric
# ============================================================================

lat = hexagonal_lattice(a)
geo = Geometry(lat, substrate, [(Circle([0.0, 0.0], r), air)])

# ============================================================================
# Computational parameters
# ============================================================================

resolution = (64, 64)
cutoff = 9
nbands = 8

println("Computational parameters:")
@printf("  Resolution: %d × %d\n", resolution...)
println("  Cutoff: $cutoff")
println("  Bands: $nbands")

# Create solvers for TM and TE modes
solver_tm = Solver(TMWave(), geo, resolution; cutoff=cutoff)
solver_te = Solver(TEWave(), geo, resolution; cutoff=cutoff)

println("  Plane waves: $(solver_tm.basis.num_pw)")
println()

# ============================================================================
# K-path: Γ → M → K → Γ (triangular lattice)
# ============================================================================

npoints = 30
kpath = simple_kpath_hexagonal(; a=a, npoints=npoints)

println("Computing band structure...")
println("  TM modes (E_z polarization)...")
bands_tm = compute_bands(solver_tm, kpath; bands=1:nbands)

println("  TE modes (H_z polarization)...")
bands_te = compute_bands(solver_te, kpath; bands=1:nbands)

# ============================================================================
# Convert to normalized frequency (ωa/2πc)
# ============================================================================

freqs_tm = bands_tm.frequencies ./ (2π)
freqs_te = bands_te.frequencies ./ (2π)

# ============================================================================
# Band gap analysis
# ============================================================================

println()
println("-" ^ 70)
println("Band Gap Analysis")
println("-" ^ 70)

# TM gaps
gaps_tm = find_all_gaps(bands_tm; threshold=0.01)
println("\nTM modes:")
tm_gap_lower = Inf
tm_gap_upper = 0.0
if !isempty(gaps_tm)
    for (i, gap) in enumerate(gaps_tm)
        lower = gap.max_lower / (2π)
        upper = gap.min_upper / (2π)
        width = upper - lower
        pct = 100 * width / ((lower + upper) / 2)
        @printf(
            "  Gap %d (bands %d-%d): [%.4f, %.4f], width=%.4f (%.1f%%)\n",
            i,
            gap.bands[1],
            gap.bands[2],
            lower,
            upper,
            width,
            pct
        )
        if i == 1
            global tm_gap_lower = lower
            global tm_gap_upper = upper
        end
    end
else
    println("  No complete TM gap found")
end

# TE gaps
gaps_te = find_all_gaps(bands_te; threshold=0.01)
println("\nTE modes:")
te_gap_lower = Inf
te_gap_upper = 0.0
if !isempty(gaps_te)
    for (i, gap) in enumerate(gaps_te)
        lower = gap.max_lower / (2π)
        upper = gap.min_upper / (2π)
        width = upper - lower
        pct = 100 * width / ((lower + upper) / 2)
        @printf(
            "  Gap %d (bands %d-%d): [%.4f, %.4f], width=%.4f (%.1f%%)\n",
            i,
            gap.bands[1],
            gap.bands[2],
            lower,
            upper,
            width,
            pct
        )
        if i == 1
            global te_gap_lower = lower
            global te_gap_upper = upper
        end
    end
else
    println("  No complete TE gap found")
end

# Check for complete gap (overlap of TM and TE)
println("\nComplete photonic band gap:")
complete_lower = max(tm_gap_lower, te_gap_lower)
complete_upper = min(tm_gap_upper, te_gap_upper)
if complete_upper > complete_lower
    width = complete_upper - complete_lower
    pct = 100 * width / ((complete_lower + complete_upper) / 2)
    @printf(
        "  Complete gap: [%.4f, %.4f], width=%.4f (%.1f%%)\n",
        complete_lower,
        complete_upper,
        width,
        pct
    )
else
    println("  No complete gap (TM and TE gaps do not overlap)")
end

println()
println("Book reference (Fig.10):")
println("  Complete photonic band gap around ω ≈ 0.4-0.5")
println("  Both TM and TE gaps overlap")

# ============================================================================
# Plot
# ============================================================================

println()
println("-" ^ 70)
println("Generating plot...")
println("-" ^ 70)

# Get distances for x-axis
distances = bands_tm.distances
labels = bands_tm.labels
xticks_pos = [distances[idx] for (idx, _) in labels]
xticks_labels = [lbl for (_, lbl) in labels]

# TM plot
p_tm = plot(;
    size=(700, 500),
    legend=false,
    xlabel="Wave vector",
    ylabel="Frequency ωa/2πc",
    title="TM Bands",
    ylims=(0, 0.8),
)
for band in 1:nbands
    plot!(p_tm, distances, freqs_tm[:, band]; color=:blue, lw=1.5)
end
if !isempty(gaps_tm)
    gap = gaps_tm[1]
    lower = gap.max_lower / (2π)
    upper = gap.min_upper / (2π)
    hspan!(p_tm, [lower, upper]; alpha=0.2, color=:blue, label="")
end
xticks!(p_tm, xticks_pos, xticks_labels)
vline!(p_tm, xticks_pos; color=:gray, lw=0.5, label="")

# TE plot
p_te = plot(;
    size=(700, 500),
    legend=false,
    xlabel="Wave vector",
    ylabel="Frequency ωa/2πc",
    title="TE Bands",
    ylims=(0, 0.8),
)
for band in 1:nbands
    plot!(p_te, distances, freqs_te[:, band]; color=:red, lw=1.5)
end
if !isempty(gaps_te)
    gap = gaps_te[1]
    lower = gap.max_lower / (2π)
    upper = gap.min_upper / (2π)
    hspan!(p_te, [lower, upper]; alpha=0.2, color=:red, label="")
end
xticks!(p_te, xticks_pos, xticks_labels)
vline!(p_te, xticks_pos; color=:gray, lw=0.5, label="")

# Combined plot
p_combined = plot(
    p_tm,
    p_te;
    layout=(1, 2),
    size=(1200, 500),
    plot_title="Joannopoulos Ch.5 Fig.10 - Triangular Lattice Air Holes (ε=13)",
)

savefig(p_combined, "examples/912_joannopoulos_ch5_fig10.png")
println("Saved: examples/912_joannopoulos_ch5_fig10.png")

# ============================================================================
# Summary
# ============================================================================

println()
println("=" ^ 70)
println("Summary")
println("=" ^ 70)
println("""

Comparison with Joannopoulos Book (Chapter 5, Figure 10):

Expected (from book):
  - Complete photonic band gap (both TM and TE)
  - Gap approximately at ω ≈ 0.4-0.5

Physical interpretation:
  - Air holes in high-ε substrate creates "veins" of dielectric
  - TE modes can be guided along these connected veins
  - TM modes are also confined, but mechanism differs
  - Triangular lattice with large holes (r/a ≈ 0.45) gives complete gap

Structure comparison:
  - Rods in air (Fig.2): TM gap only
  - Holes in dielectric (Fig.10): Complete gap (TM + TE overlap)
""")

println("Example complete!")
