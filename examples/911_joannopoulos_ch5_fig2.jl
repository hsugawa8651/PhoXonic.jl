# Last-Modified: 2025-12-15T19:00:00+09:00
# Joannopoulos Book Chapter 5, Figure 2 Reproduction
#
# Reference: "Photonic Crystals: Molding the Flow of Light" (2nd edition)
# ISBN: 978-0-691-12456-8, http://ab-initio.mit.edu/book/
#            Chapter 5: Two-Dimensional Photonic Crystals
#            Figure 2 (page 68)
#
# This example reproduces the photonic band structure for a square array
# of dielectric columns (alumina) in air.
#
# Parameters from the book:
#   - Square lattice
#   - Dielectric columns: ε = 8.9 (alumina)
#   - Background: air (ε = 1)
#   - Column radius: r = 0.2a
#   - Blue bands = TM modes
#   - Red bands = TE modes

using PhoXonic
using Printf
using Plots

println("=" ^ 70)
println("Joannopoulos Book - Chapter 5, Figure 2")
println("2D Photonic Band Structure: Square Array of Dielectric Columns")
println("=" ^ 70)
println()

# ============================================================================
# Parameters (matching the book)
# ============================================================================

a = 1.0  # Lattice constant (normalized)

# Materials (from book caption: ε = 8.9, as for alumina)
ε_alumina = 8.9
ε_air = 1.0

alumina = Dielectric(ε_alumina)
air = Dielectric(ε_air)

# Column radius: r = 0.2a
r = 0.2 * a

println("Parameters (from Joannopoulos Ch.5, Fig.2):")
println("  Lattice: 2D square, a = $a")
println("  Structure: Dielectric columns in air")
@printf("  Column material: ε = %.1f (alumina)\n", ε_alumina)
@printf("  Background: ε = %.1f (air)\n", ε_air)
@printf("  Column radius: r/a = %.2f\n", r/a)
println()

# ============================================================================
# Geometry: Dielectric rods in air
# ============================================================================

lat = square_lattice(a)
geo = Geometry(lat, air, [(Circle([0.0, 0.0], r), alumina)])

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
# K-path: Γ → X → M → Γ
# ============================================================================

npoints = 30
kpath = simple_kpath_square(a=a, npoints=npoints)

println("Computing band structure...")
println("  TM modes (E_z polarization)...")
bands_tm = compute_bands(solver_tm, kpath; bands=1:nbands)

println("  TE modes (H_z polarization)...")
bands_te = compute_bands(solver_te, kpath; bands=1:nbands)

# ============================================================================
# Convert to normalized frequency (ωa/2πc)
# ============================================================================

# PhoXonic returns angular frequency ω, divide by 2π to get normalized frequency
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
if !isempty(gaps_tm)
    for (i, gap) in enumerate(gaps_tm)
        lower = gap.max_lower / (2π)
        upper = gap.min_upper / (2π)
        width = upper - lower
        pct = 100 * width / ((lower + upper) / 2)
        @printf("  Gap %d (bands %d-%d): [%.4f, %.4f], width=%.4f (%.1f%%)\n",
                i, gap.bands[1], gap.bands[2], lower, upper, width, pct)
    end
else
    println("  No complete TM gap found")
end

# TE gaps
gaps_te = find_all_gaps(bands_te; threshold=0.01)
println("\nTE modes:")
if !isempty(gaps_te)
    for (i, gap) in enumerate(gaps_te)
        lower = gap.max_lower / (2π)
        upper = gap.min_upper / (2π)
        width = upper - lower
        pct = 100 * width / ((lower + upper) / 2)
        @printf("  Gap %d (bands %d-%d): [%.4f, %.4f], width=%.4f (%.1f%%)\n",
                i, gap.bands[1], gap.bands[2], lower, upper, width, pct)
    end
else
    println("  No complete TE gap found")
end

println()
println("Book reference (Fig.2):")
println("  TM gap: approximately [0.30, 0.45]")
println("  TE gap: none (for rods in air)")

# ============================================================================
# Plot
# ============================================================================

println()
println("-" ^ 70)
println("Generating plot...")
println("-" ^ 70)

# Get distances for x-axis
distances = bands_tm.distances

# Create plot
p = plot(size=(700, 500), legend=:topright,
         xlabel="Wave vector", ylabel="Frequency ωa/2πc",
         title="Joannopoulos Ch.5 Fig.2 - Square Lattice (ε=8.9, r/a=0.2)")

# Plot TM bands (blue)
for band in 1:nbands
    plot!(p, distances, freqs_tm[:, band], color=:blue, lw=1.5,
          label=(band == 1 ? "TM" : ""))
end

# Plot TE bands (red)
for band in 1:nbands
    plot!(p, distances, freqs_te[:, band], color=:red, lw=1.5,
          label=(band == 1 ? "TE" : ""))
end

# Add TM gap shading if exists
if !isempty(gaps_tm)
    gap = gaps_tm[1]
    lower = gap.max_lower / (2π)
    upper = gap.min_upper / (2π)
    hspan!(p, [lower, upper], alpha=0.2, color=:blue, label="TM gap")
end

# Set axis limits
ylims!(p, 0, 0.8)

# Add high-symmetry point labels
labels = bands_tm.labels
xticks_pos = [distances[idx] for (idx, _) in labels]
xticks_labels = [lbl for (_, lbl) in labels]
xticks!(p, xticks_pos, xticks_labels)

# Add vertical lines at high-symmetry points
for (idx, _) in labels
    vline!(p, [distances[idx]], color=:gray, lw=0.5, label="")
end

savefig(p, "examples/911_joannopoulos_ch5_fig2.png")
println("Saved: examples/911_joannopoulos_ch5_fig2.png")

# ============================================================================
# Summary
# ============================================================================

println()
println("=" ^ 70)
println("Summary")
println("=" ^ 70)
println("""

Comparison with Joannopoulos Book (Chapter 5, Figure 2):

Expected (from book):
  - TM modes: Gap between bands 1-2 around ω ≈ 0.30-0.45
  - TE modes: No gap (for dielectric rods in air)

Physical interpretation:
  - TM gap exists because E-field concentrates in high-ε rods
  - TE gap absent because H-field cannot be confined in isolated rods

Note: If PhoXonic results differ significantly from the book,
      this indicates a potential issue with the 2D solver that
      requires investigation (see 94_mpb_benchmark.jl for details).
""")

println("Example complete!")
