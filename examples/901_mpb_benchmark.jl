# Last-Modified: 2025-12-15T19:00:00+09:00
# MPB Benchmark Comparison
#
# Compare PhoXonic band structure results with MIT Photonic Bands (MPB) tutorial.
# Reference: https://mpb.readthedocs.io/en/latest/Scheme_Tutorial/
#
# MPB Tutorial Parameters:
#   - 2D square lattice of dielectric rods in air
#   - Dielectric constant: ε = 12 (similar to Si)
#   - Rod radius: r/a = 0.2
#   - TM polarization
#
# MPB Reference Values (TM mode):
#   - Band 1 top: 0.2826 (ωa/2πc)
#   - Band 2 bottom: 0.4193 (ωa/2πc)
#   - Gap: 38.95%

using PhoXonic
using LinearAlgebra
using Printf

println("=" ^ 70)
println("MPB Benchmark: 2D Square Lattice of Dielectric Rods in Air")
println("=" ^ 70)
println()

# ============================================================================
# MPB Tutorial Parameters (exactly matching)
# ============================================================================

a = 1.0  # Lattice constant (normalized)
epsilon_rod = 12.0  # Dielectric constant (MPB uses 12)
r = 0.2 * a  # Rod radius r/a = 0.2

# Materials
rod_material = Dielectric(epsilon_rod)
air = Dielectric(1.0)

# Geometry: Dielectric rods in air background
lat = square_lattice(a)
geo = Geometry(lat, air, [(Circle([0.0, 0.0], r), rod_material)])

println("Parameters (matching MPB tutorial):")
println("  Lattice: 2D square")
println("  Structure: Dielectric rods in air")
@printf("  Dielectric constant: ε = %.1f\n", epsilon_rod)
@printf("  Rod radius: r/a = %.2f\n", r / a)
println()

# ============================================================================
# Solver Setup
# ============================================================================

# Higher resolution for accuracy
resolution = (32, 32)  # MPB uses 32x32
cutoff = 7

# TM mode (E-field along rod axis)
solver_tm = Solver(TMWave(), geo, resolution; cutoff=cutoff)

# TE mode for comparison
solver_te = Solver(TEWave(), geo, resolution; cutoff=cutoff)

println("Computational parameters:")
@printf("  Resolution: %d × %d\n", resolution...)
@printf("  Plane waves (cutoff=%d): %d\n", cutoff, solver_tm.basis.num_pw)
println()

# ============================================================================
# Band Structure Calculation
# ============================================================================

println("Computing band structure...")
println()

# k-path: Γ → X → M → Γ
kpath = simple_kpath_square(a=a, npoints=30)

# Compute bands
bands_tm = compute_bands(solver_tm, kpath; bands=1:8)
bands_te = compute_bands(solver_te, kpath; bands=1:8)

# ============================================================================
# Extract Band Frequencies at High-Symmetry Points
# ============================================================================

println("-" ^ 70)
println("TM Band Frequencies at High-Symmetry Points")
println("-" ^ 70)
println()

# Find indices for high-symmetry points
# kpath goes: Γ(0) → X(1/3) → M(2/3) → Γ(1)
n_kpts = length(kpath.points)
idx_gamma1 = 1
idx_X = n_kpts ÷ 3
idx_M = 2 * n_kpts ÷ 3
idx_gamma2 = n_kpts

# Extract TM frequencies
# frequencies is (nk × nbands) matrix
function get_frequencies_at_k(bands, idx)
    bands.frequencies[idx, :]
end

freq_tm_gamma = get_frequencies_at_k(bands_tm, idx_gamma1)
freq_tm_X = get_frequencies_at_k(bands_tm, idx_X)
freq_tm_M = get_frequencies_at_k(bands_tm, idx_M)

println("TM mode frequencies (ωa/2πc):")
println()
println("  Band |    Γ     |    X     |    M    ")
println("  -----|----------|----------|----------")
for i in 1:min(6, length(freq_tm_gamma))
    @printf("   %d   | %8.5f | %8.5f | %8.5f\n",
            i, freq_tm_gamma[i], freq_tm_X[i], freq_tm_M[i])
end
println()

# ============================================================================
# Band Gap Analysis (TM mode)
# ============================================================================

println("-" ^ 70)
println("TM Band Gap Analysis")
println("-" ^ 70)
println()

gaps_tm = find_all_gaps(bands_tm; threshold=0.01)

if !isempty(gaps_tm)
    gap = gaps_tm[1]
    lower_band, upper_band = gap.bands
    gap_center = (gap.max_lower + gap.min_upper) / 2
    gap_width = gap.min_upper - gap.max_lower
    gap_percent = 100 * gap_width / gap_center

    println("PhoXonic results:")
    @printf("  Band %d top:     %.6f\n", lower_band, gap.max_lower)
    @printf("  Band %d bottom:  %.6f\n", upper_band, gap.min_upper)
    @printf("  Gap width:       %.4f (%.2f%%)\n", gap_width, gap_percent)
    println()

    # Note on frequency units:
    # PhoXonic uses angular frequency ω = 2πc/λ with c=1, so ω = 2π/λ
    # MPB uses reduced frequency ωa/2πc = a/λ
    # To convert PhoXonic → MPB: divide by 2π
    ω_lower_mpb = gap.max_lower / (2π)
    ω_upper_mpb = gap.min_upper / (2π)

    println("Converted to MPB units (ωa/2πc = a/λ):")
    @printf("  Band %d top:     %.6f\n", lower_band, ω_lower_mpb)
    @printf("  Band %d bottom:  %.6f\n", upper_band, ω_upper_mpb)
    println()

    # MPB reference values
    println("MPB reference values (from tutorial):")
    println("  Band 1 top:     0.282623")
    println("  Band 2 bottom:  0.419335")
    println("  Gap width:      0.1367 (38.95%)")
    println()

    # Relative error (using converted values)
    err_lower = abs(ω_lower_mpb - 0.282623) / 0.282623 * 100
    err_upper = abs(ω_upper_mpb - 0.419335) / 0.419335 * 100

    println("Relative errors (vs MPB):")
    @printf("  Band 1 top:     %.2f%%\n", err_lower)
    @printf("  Band 2 bottom:  %.2f%%\n", err_upper)
else
    println("No complete TM bandgap found (unexpected!)")
    println("Check parameters or increase resolution/cutoff")
end

# ============================================================================
# TE Mode (for reference - no gap expected for rods)
# ============================================================================

println()
println("-" ^ 70)
println("TE Band Gap Analysis (reference)")
println("-" ^ 70)
println()

gaps_te = find_all_gaps(bands_te; threshold=0.01)

if !isempty(gaps_te)
    gap_te = gaps_te[1]
    gap_percent_te = 100 * gap_te.gap_ratio
    # Convert to MPB units
    lower_mpb = gap_te.max_lower / (2π)
    upper_mpb = gap_te.min_upper / (2π)
    @printf("TE gap found: [%.4f, %.4f] (%.2f%%) in MPB units\n",
            lower_mpb, upper_mpb, gap_percent_te)
else
    println("No TE bandgap found (expected for rods in air)")
end

# ============================================================================
# Summary
# ============================================================================

println()
println("=" ^ 70)
println("Summary")
println("=" ^ 70)
println("""

⚠️  DISCREPANCY DETECTED ⚠️

PhoXonic results do not match MPB reference values.
This requires further investigation.

Observed:
  - TM gap position is shifted upward (~0.41 vs ~0.28)
  - Gap width is narrower (18-21% vs 38.95%)
  - Convergence tests show this is not a resolution issue

Possible causes to investigate:
  1. Eigenvalue formulation difference (E-field vs H-field)
  2. Reciprocal lattice vector convention
  3. Boundary of Brillouin zone definition
  4. Dielectric function Fourier coefficients

Key physics (expected):
- Dielectric rods in air → TM gap (field concentrated in rods)
- Air holes in dielectric → TE gap (field concentrated in veins)

References:
  MPB Tutorial: https://mpb.readthedocs.io/en/latest/Scheme_Tutorial/
  Joannopoulos Book: http://ab-initio.mit.edu/book/
""")

println("Benchmark complete!")
