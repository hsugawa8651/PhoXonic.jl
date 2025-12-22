# Last-Modified: 2025-12-15T23:00:00+09:00
# Defect Mode Analysis with Matrix-Free RSCG
# Comparison of Dense vs Matrix-free LDOS calculation
#
# REQUIRES: ReducedShiftedKrylov.jl (loaded automatically as package extension)
#
# This example demonstrates the unified Green's function API:
#   compute_ldos(solver, pos, ω_values, k_points, method; η=...)
#
# Available methods (GFMethod):
#   - DirectGF()      - Dense direct solve (most accurate)
#   - RSKGF()         - Dense RSK using ReducedShiftedKrylov.jl
#   - MatrixFreeGF()  - Matrix-free RSCG (O(N) memory)
#   - MatrixFreeGF(rhs_inv_method=CGRHSInv()) - Matrix-free with exact RHS⁻¹
#
# IMPORTANT: Practical Usage Guidelines
# =====================================
# MatrixFreeGF() is recommended for:
#   - Peak/resonance frequency detection (accurate even without full RSCG convergence)
#   - Large-scale calculations where Dense method is memory-prohibitive
#   - Qualitative spectral analysis
#
# DirectGF() is recommended for:
#   - Accurate absolute LDOS/DOS values
#   - Quantitative analysis requiring high precision
#   - Smaller systems where O(N²) memory is acceptable

using PhoXonic
using ReducedShiftedKrylov  # Required for MatrixFreeGF, RSKGF
using LinearAlgebra
using StaticArrays
using Plots
using Printf

println("=== Defect Mode: Dense vs Matrix-Free Comparison ===")
println()

# ============================================================================
# Setup: Same as 91_defect_mode.jl but with comparison
# ============================================================================

a = 1.0  # Lattice constant
lat = square_lattice(a)

# Materials: Low contrast for matrix-free approximation accuracy
# (High contrast like Si/air would have larger approximation error)
air = Dielectric(1.0)
rod_material = Dielectric(2.5)  # Low contrast: ε=2.5 vs ε=1.0
r = 0.3  # Rod radius

# Perfect crystal for reference bandgap
geo_perfect = Geometry(lat, air, [(Circle([0.0, 0.0], r), rod_material)])
solver_perfect = Solver(TMWave(), geo_perfect, (32, 32); cutoff=5)
kpath = simple_kpath_square(; a=a, npoints=30)
bands_perfect = compute_bands(solver_perfect, kpath; bands=1:8)

gaps = find_all_gaps(bands_perfect; threshold=0.01)
if !isempty(gaps)
    gap = gaps[1]
    ω_min = gap.max_lower
    ω_max = gap.min_upper
    println("Bandgap: [$(round(ω_min, digits=4)), $(round(ω_max, digits=4))]")
else
    ω_min, ω_max = 0.25, 0.45
    println("Using fallback gap range")
end

# ============================================================================
# Function to create supercell with defect
# ============================================================================

function create_defect_supercell(N_super, a, r, rod_material, air)
    lat_super = Lattice(SVector(N_super * a, 0.0), SVector(0.0, N_super * a))

    inclusions = Tuple{Circle,Dielectric}[]
    center_idx = (N_super + 1) / 2

    for i in 1:N_super
        for j in 1:N_super
            # Skip center rod (point defect)
            if i == ceil(Int, center_idx) && j == ceil(Int, center_idx)
                continue
            end
            x = (i - 0.5) * a
            y = (j - 0.5) * a
            push!(inclusions, (Circle([x, y], r), rod_material))
        end
    end

    return Geometry(lat_super, air, inclusions)
end

# ============================================================================
# Part 1: Small supercell - Dense vs Matrix-free comparison
# ============================================================================

println("\n" * "="^60)
println("Part 1: 5x5 Supercell - Dense vs Matrix-free Comparison")
println("="^60)

N_super = 5
geo_defect = create_defect_supercell(N_super, a, r, rod_material, air)

resolution_per_cell = 16
resolution = (N_super * resolution_per_cell, N_super * resolution_per_cell)
cutoff = 5

solver = Solver(TMWave(), geo_defect, resolution; cutoff=cutoff)

println("\nSupercell: $(N_super)x$(N_super)")
println("Resolution: $resolution")
println("Plane waves: $(solver.basis.num_pw)")
println("Matrix size: $(solver.basis.num_pw) x $(solver.basis.num_pw)")

# Parameters
defect_pos = [N_super * a / 2, N_super * a / 2]
n_freq = 30
ω_range = collect(range(ω_min * 0.9, ω_max * 1.1; length=n_freq))
k_points = [[0.0, 0.0]]
η = 0.02

println("\nFrequency points: $n_freq")
println("Broadening η: $η")

# --- Dense LDOS (DirectGF) ---
println("\n--- DirectGF() ---")
t_dense = @elapsed begin
    ldos_dense = compute_ldos(solver, defect_pos, ω_range, k_points, DirectGF(); η=η)
end
@printf("  Time: %.2f sec\n", t_dense)

# --- Matrix-free LDOS (approximate) ---
println("\n--- MatrixFreeGF() ---")
t_mf_approx = @elapsed begin
    ldos_mf_approx = compute_ldos(
        solver, defect_pos, ω_range, k_points, MatrixFreeGF(); η=η
    )
end
@printf("  Time: %.2f sec\n", t_mf_approx)

# --- Matrix-free LDOS (CG) ---
println("\n--- MatrixFreeGF(rhs_inv_method=CGRHSInv()) ---")
t_mf_cg = @elapsed begin
    ldos_mf_cg = compute_ldos(
        solver, defect_pos, ω_range, k_points, MatrixFreeGF(rhs_inv_method=CGRHSInv()); η=η
    )
end
@printf("  Time: %.2f sec\n", t_mf_cg)

# --- Comparison ---
rel_error_approx = norm(ldos_dense - ldos_mf_approx) / norm(ldos_dense)
rel_error_cg = norm(ldos_dense - ldos_mf_cg) / norm(ldos_dense)
@printf("\nRelative error (MatrixFreeGF):           %.2e\n", rel_error_approx)
@printf("Relative error (MatrixFreeGF :cg):       %.2e\n", rel_error_cg)
@printf("Speedup (MatrixFreeGF):     %.1fx\n", t_dense / t_mf_approx)
@printf("Speedup (MatrixFreeGF :cg): %.1fx\n", t_dense / t_mf_cg)

# Find defect mode
idx_max_dense = argmax(ldos_dense)
idx_max_approx = argmax(ldos_mf_approx)
idx_max_cg = argmax(ldos_mf_cg)
ω_defect_dense = ω_range[idx_max_dense]
ω_defect_approx = ω_range[idx_max_approx]
ω_defect_cg = ω_range[idx_max_cg]

println("\nDefect mode frequency:")
println("  Dense:       $(round(ω_defect_dense, digits=4))")
println("  MF approx:   $(round(ω_defect_approx, digits=4))")
println("  MF CG:       $(round(ω_defect_cg, digits=4))")

# ============================================================================
# Part 2: Larger supercell - Matrix-free advantage
# ============================================================================

println("\n" * "="^60)
println("Part 2: 7x7 Supercell - Matrix-free for Large Systems")
println("="^60)

N_super_large = 7
geo_defect_large = create_defect_supercell(N_super_large, a, r, rod_material, air)

resolution_large = (
    N_super_large * resolution_per_cell, N_super_large * resolution_per_cell
)
solver_large = Solver(TMWave(), geo_defect_large, resolution_large; cutoff=cutoff)

println("\nSupercell: $(N_super_large)x$(N_super_large)")
println("Resolution: $resolution_large")
println("Plane waves: $(solver_large.basis.num_pw)")

N_pw = solver_large.basis.num_pw
dense_memory_mb = N_pw^2 * 16 / 1e6  # ComplexF64 = 16 bytes
mf_memory_mb = N_pw * 16 / 1e6 * 10   # ~10 arrays of size N
println("\nEstimated memory:")
@printf("  Dense matrices: %.1f MB\n", dense_memory_mb)
@printf("  Matrix-free:    %.1f MB\n", mf_memory_mb)

defect_pos_large = [N_super_large * a / 2, N_super_large * a / 2]

println("\n--- MatrixFreeGF (7x7 supercell, CGRHSInv) ---")
t_mf_large = @elapsed begin
    ldos_mf_large = compute_ldos(
        solver_large,
        defect_pos_large,
        ω_range,
        k_points,
        MatrixFreeGF(rhs_inv_method=CGRHSInv());
        η=η,
    )
end
@printf("  Time: %.2f sec\n", t_mf_large)

idx_max_large = argmax(ldos_mf_large)
ω_defect_large = ω_range[idx_max_large]
println("  Defect mode: $(round(ω_defect_large, digits=4))")

# ============================================================================
# Visualization
# ============================================================================

println("\n" * "="^60)
println("Creating Plots")
println("="^60)

# Plot 1: Dense vs Matrix-free comparison (5x5)
p1 = plot(;
    xlabel="Frequency (2πc/a)",
    ylabel="LDOS (arb. units)",
    title="5×5 Supercell: Dense vs Matrix-free",
    legend=:topright,
    grid=true,
    size=(600, 400),
)

plot!(p1, ω_range * a / (2π), ldos_dense; label="Dense", linewidth=2, color=:blue)
plot!(
    p1,
    ω_range * a / (2π),
    ldos_mf_approx;
    label="MF approx",
    linewidth=2,
    color=:red,
    linestyle=:dash,
)
plot!(
    p1,
    ω_range * a / (2π),
    ldos_mf_cg;
    label="MF CG",
    linewidth=2,
    color=:green,
    linestyle=:dot,
)

vspan!(p1, [ω_min * a / (2π), ω_max * a / (2π)]; alpha=0.2, color=:yellow, label="Bandgap")

# Plot 2: Large supercell (7x7)
p2 = plot(;
    xlabel="Frequency (2πc/a)",
    ylabel="LDOS (arb. units)",
    title="7×7 Supercell: Matrix-free LDOS",
    legend=:topright,
    grid=true,
    size=(600, 400),
)

plot!(
    p2,
    ω_range * a / (2π),
    ldos_mf_large;
    label="Matrix-free (7×7)",
    linewidth=2,
    color=:green,
)

vspan!(p2, [ω_min * a / (2π), ω_max * a / (2π)]; alpha=0.2, color=:yellow, label="Bandgap")

vline!(
    p2,
    [ω_defect_large * a / (2π)];
    color=:red,
    linestyle=:dash,
    linewidth=2,
    label="Defect mode",
)

# Plot 3: Supercell size comparison
p3 = plot(;
    xlabel="Frequency (2πc/a)",
    ylabel="LDOS (normalized)",
    title="Supercell Size Comparison",
    legend=:topright,
    grid=true,
    size=(600, 400),
)

# Normalize for comparison
ldos_mf_cg_norm = ldos_mf_cg / maximum(ldos_mf_cg)
ldos_mf_large_norm = ldos_mf_large / maximum(ldos_mf_large)

plot!(
    p3,
    ω_range * a / (2π),
    ldos_mf_cg_norm;
    label="5×5 supercell (CG)",
    linewidth=2,
    color=:blue,
)
plot!(
    p3,
    ω_range * a / (2π),
    ldos_mf_large_norm;
    label="7×7 supercell (CG)",
    linewidth=2,
    color=:green,
)

vspan!(p3, [ω_min * a / (2π), ω_max * a / (2π)]; alpha=0.2, color=:yellow, label="Bandgap")

# Combined plot
p_combined = plot(p1, p2, p3; layout=(1, 3), size=(1500, 400))

# Save plots
savefig(p1, joinpath(@__DIR__, "502_dense_vs_matrixfree.png"))
savefig(p2, joinpath(@__DIR__, "502_large_supercell.png"))
savefig(p3, joinpath(@__DIR__, "502_supercell_comparison.png"))
savefig(p_combined, joinpath(@__DIR__, "502_combined.png"))

println("\nSaved: 502_dense_vs_matrixfree.png")
println("Saved: 502_large_supercell.png")
println("Saved: 502_supercell_comparison.png")
println("Saved: 502_combined.png")

# ============================================================================
# Summary
# ============================================================================

println("\n" * "="^60)
println("Summary")
println("="^60)

println("\n5×5 Supercell:")
@printf("  DirectGF():             %.2f sec\n", t_dense)
@printf("  MatrixFreeGF():         %.2f sec (error: %.2e)\n", t_mf_approx, rel_error_approx)
@printf("  MatrixFreeGF(:cg):      %.2f sec (error: %.2e)\n", t_mf_cg, rel_error_cg)

println("\n7×7 Supercell:")
@printf("  MatrixFreeGF(:cg): %.2f sec\n", t_mf_large)
println("  (Dense would require ~$(round(Int, dense_memory_mb)) MB for matrices)")

println("\nDefect mode frequencies:")
println("  5×5 DirectGF:       $(round(ω_defect_dense, digits=4))")
println("  5×5 MatrixFreeGF:   $(round(ω_defect_approx, digits=4))")
println("  5×5 MatrixFreeGF(:cg): $(round(ω_defect_cg, digits=4))")
println("  7×7 MatrixFreeGF(:cg): $(round(ω_defect_large, digits=4))")

println("\n✅ Matrix-free LDOS example complete!")
