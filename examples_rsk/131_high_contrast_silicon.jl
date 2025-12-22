# Last-Modified: 2025-12-15T23:00:00+09:00
# High-Contrast Si/Air Photonic Crystal - RHSInvMethod Comparison
#
# REQUIRES: ReducedShiftedKrylov.jl (loaded automatically as package extension)
#
# This example demonstrates the new type-based RHSInvMethod API.
#
# Key concepts:
#   - Si/Air has high dielectric contrast (epsilon ratio ~ 12:1)
#   - ApproximateRHSInv(): Fast, element-wise 1/epsilon in real space
#   - CGRHSInv(): Accurate iterative solve, with tunable parameters
#
# New type-based API (recommended):
#   MatrixFreeGF(rhs_inv_method=ApproximateRHSInv())
#   MatrixFreeGF(rhs_inv_method=CGRHSInv(atol=1e-12, rtol=1e-12, maxiter=200))
#
# Deprecated Symbol API (still works with warning):
#   MatrixFreeGF(rhs_inv_method=:approximate)
#   MatrixFreeGF(rhs_inv_method=:cg)

using PhoXonic
using ReducedShiftedKrylov  # Required for MatrixFreeGF
using LinearAlgebra
using StaticArrays
using Printf

println("=== High-Contrast Si/Air Photonic Crystal ===")
println("Demonstrating RHSInvMethod for matrix-free Green's function")
println()

# ============================================================================
# Setup: 2D Square Lattice with Circular Air Holes in Silicon
# ============================================================================

a = 1.0  # Lattice constant

# Materials: High contrast Si/Air
epsilon_Si = 11.56  # Si at optical frequencies
epsilon_air = 1.0
silicon = Dielectric(epsilon_Si)
air = Dielectric(epsilon_air)

# Geometry: Air holes in silicon background
r = 0.3 * a  # Hole radius
lat = square_lattice(a)
geo = Geometry(lat, silicon, [(Circle([0.0, 0.0], r), air)])

@printf("Material contrast: epsilon_Si / epsilon_air = %.1f\n", epsilon_Si / epsilon_air)
println("Hole radius: r/a = $(r/a)")
println()

# Solver setup - use higher resolution to see contrast effects
resolution = (64, 64)
cutoff = 7
solver = Solver(TEWave(), geo, resolution; cutoff=cutoff)

println("Resolution: $resolution")
println("Plane waves: $(solver.basis.num_pw)")
println()

# ============================================================================
# Band Structure (for reference)
# ============================================================================

println("Computing band structure...")
kpath = simple_kpath_square(; a=a, npoints=20)
bands = compute_bands(solver, kpath; bands=1:6)

# Find bandgap
gaps = find_all_gaps(bands; threshold=0.01)
if !isempty(gaps)
    gap = gaps[1]
    @printf("TE bandgap: [%.4f, %.4f]\n", gap.max_lower, gap.min_upper)
    ω_center = (gap.max_lower + gap.min_upper) / 2
    ω_width = gap.min_upper - gap.max_lower
else
    println("No complete bandgap found, using frequency range near band edge")
    ω_center = 0.30
    ω_width = 0.10
end

# ============================================================================
# Warmup (JIT compilation)
# ============================================================================

println("\nWarmup (JIT compilation)...")
k_warmup = [0.0, 0.0]
ω_warmup = [ω_center]
source_warmup = zeros(ComplexF64, solver.basis.num_pw)
source_warmup[1] = 1.0

compute_greens_function(solver, k_warmup, ω_warmup, source_warmup, DirectGF(); η=0.01)
compute_greens_function(
    solver,
    k_warmup,
    ω_warmup,
    source_warmup,
    MatrixFreeGF(; rhs_inv_method=ApproximateRHSInv());
    η=0.01,
)
compute_greens_function(
    solver,
    k_warmup,
    ω_warmup,
    source_warmup,
    MatrixFreeGF(; rhs_inv_method=CGRHSInv());
    η=0.01,
)
println("Warmup complete.")

# ============================================================================
# Green's Function: RHSInvMethod Comparison
# ============================================================================

println("\n" * "="^60)
println("Green's Function: RHSInvMethod Comparison")
println("="^60)

# Parameters
k = [0.1, 0.1]  # Off-center k-point
n_freq = 30
ω_values = collect(range(ω_center - ω_width, ω_center + ω_width; length=n_freq))
η = 0.005  # Small broadening for sharper features

# Source vector (point source at origin)
source = zeros(ComplexF64, solver.basis.num_pw)
source[1] = 1.0  # G=0 component

println("\nk-point: $k")
println("Frequency points: $n_freq")
println("Broadening: η = $η")

# --- Reference: DirectGF ---
println("\n--- DirectGF (reference) ---")
t_direct = @elapsed begin
    G_direct = compute_greens_function(solver, k, ω_values, source, DirectGF(); η=η)
end
ldos_direct = [-imag(dot(source, G)) / π for G in G_direct]
@printf("  Time: %.3f sec\n", t_direct)

# --- ApproximateRHSInv ---
println("\n--- MatrixFreeGF(ApproximateRHSInv()) ---")
t_approx = @elapsed begin
    G_approx = compute_greens_function(
        solver, k, ω_values, source, MatrixFreeGF(rhs_inv_method=ApproximateRHSInv()); η=η
    )
end
ldos_approx = [-imag(dot(source, G)) / π for G in G_approx]
@printf("  Time: %.3f sec\n", t_approx)

# --- CGRHSInv with default parameters ---
println("\n--- MatrixFreeGF(CGRHSInv()) [default] ---")
t_cg_default = @elapsed begin
    G_cg_default = compute_greens_function(
        solver, k, ω_values, source, MatrixFreeGF(rhs_inv_method=CGRHSInv()); η=η
    )
end
ldos_cg_default = [-imag(dot(source, G)) / π for G in G_cg_default]
@printf("  Time: %.3f sec\n", t_cg_default)

# --- CGRHSInv with tighter tolerance ---
println("\n--- MatrixFreeGF(CGRHSInv(atol=1e-12)) [tight] ---")
t_cg_tight = @elapsed begin
    G_cg_tight = compute_greens_function(
        solver,
        k,
        ω_values,
        source,
        MatrixFreeGF(rhs_inv_method=CGRHSInv(atol=1e-12, rtol=1e-12, maxiter=200));
        η=η,
    )
end
ldos_cg_tight = [-imag(dot(source, G)) / π for G in G_cg_tight]
@printf("  Time: %.3f sec\n", t_cg_tight)

# ============================================================================
# Error Analysis
# ============================================================================

println("\n" * "="^60)
println("Error Analysis (vs DirectGF reference)")
println("="^60)

function relative_error(test, ref)
    norm(test - ref) / norm(ref)
end

err_approx = relative_error(ldos_approx, ldos_direct)
err_cg_default = relative_error(ldos_cg_default, ldos_direct)
err_cg_tight = relative_error(ldos_cg_tight, ldos_direct)

println()
@printf("ApproximateRHSInv:    error = %.2e\n", err_approx)
@printf("CGRHSInv (default):   error = %.2e\n", err_cg_default)
@printf("CGRHSInv (tight):     error = %.2e\n", err_cg_tight)

println()
println("Timing comparison:")
@printf("  DirectGF:           %.3f sec (reference)\n", t_direct)
@printf("  ApproximateRHSInv:  %.3f sec (%.1fx)\n", t_approx, t_direct / t_approx)
@printf("  CGRHSInv (default): %.3f sec (%.1fx)\n", t_cg_default, t_direct / t_cg_default)
@printf("  CGRHSInv (tight):   %.3f sec (%.1fx)\n", t_cg_tight, t_direct / t_cg_tight)

# ============================================================================
# Summary
# ============================================================================

println("\n" * "="^60)
println("Summary: New Type-Based RHSInvMethod API")
println("="^60)

println("""

Usage examples:
```julia
using PhoXonic

# Default: ApproximateRHSInv (fast)
method = MatrixFreeGF()

# Explicit ApproximateRHSInv
method = MatrixFreeGF(rhs_inv_method=ApproximateRHSInv())

# CGRHSInv with default parameters
method = MatrixFreeGF(rhs_inv_method=CGRHSInv())

# CGRHSInv with custom parameters
method = MatrixFreeGF(rhs_inv_method=CGRHSInv(
    atol=1e-12,    # Absolute tolerance
    rtol=1e-12,    # Relative tolerance
    maxiter=200    # Maximum iterations
))

# Then use with compute_greens_function or compute_ldos
G = compute_greens_function(solver, k, ω_values, source, method; η=1e-2)
ldos = compute_ldos(solver, position, ω_values, k_points, method; η=1e-2)
```

CGRHSInv parameters for different scenarios:
  - Standard:     CGRHSInv()                              # atol=1e-10, maxiter=100
  - High accuracy: CGRHSInv(atol=1e-12, maxiter=200)
  - Fast approx:  CGRHSInv(atol=1e-6, maxiter=20)
""")

println("Example complete!")
