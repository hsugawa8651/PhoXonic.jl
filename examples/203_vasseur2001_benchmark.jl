# Last-Modified: 2025-12-19
# Vasseur 2001: Steel/Epoxy phononic crystal (hexagonal lattice)
# Reference: M. Vasseur et al., J. Phys.: Condens. Matter 13, 8759 (2001)
# DOI: https://doi.org/10.1088/0953-8984/13/41/305
#
# Demonstrates:
# 1. Dense vs LOBPCG result consistency
# 2. LOBPCG warm start speedup

using PhoXonic
using Printf
using Plots

println("=" ^ 60)
println("Vasseur 2001: Steel/Epoxy Phononic Crystal")
println("Hexagonal lattice, P-SV wave")
println("=" ^ 60)
println()

# ============================================================================
# Material properties (from Vasseur 2001)
# ============================================================================
# Steel (inclusion)
ρ_steel = 7780.0      # kg/m³
μ_steel = 8.1e10      # Pa (shear modulus)
λ_steel = 26.4e10 - 2*μ_steel  # Pa (first Lame parameter)

# Epoxy (matrix)
ρ_epoxy = 1142.0      # kg/m³
μ_epoxy = 0.148e10    # Pa
λ_epoxy = 0.754e10 - 2*μ_epoxy  # Pa

steel = IsotropicElastic(ρ=ρ_steel, λ=λ_steel, μ=μ_steel)
epoxy = IsotropicElastic(ρ=ρ_epoxy, λ=λ_epoxy, μ=μ_epoxy)

println("Materials:")
println("  Steel: ρ = $ρ_steel kg/m³, μ = $(μ_steel/1e9) GPa")
println("  Epoxy: ρ = $ρ_epoxy kg/m³, μ = $(μ_epoxy/1e9) GPa")

# ============================================================================
# Geometry: hexagonal lattice with steel cylinders
# ============================================================================
a = 1.0  # lattice constant (normalized)
lat = hexagonal_lattice(a)

# Filling fraction f = 0.4
f_target = 0.4
r = sqrt(f_target * a^2 * sqrt(3) / (2π))

geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], r), steel)])

println("\nGeometry:")
println("  Hexagonal lattice, a = $a")
println("  Cylinder radius r = $(round(r, digits=4))")
println("  Filling fraction f = $(round(π*r^2 / (a^2 * sqrt(3)/2), digits=2))")

# ============================================================================
# K-path: Γ → M → K → Γ
# ============================================================================
npoints = 30
kpath = simple_kpath_hexagonal(a=a, npoints=npoints)
println("\nK-path: $(length(kpath.points)) points")

# ============================================================================
# Solver parameters
# ============================================================================
cutoff = 12  # Moderate size for demonstration
nbands = 1:15

println("\nSolver parameters:")
println("  Cutoff = $cutoff")
println("  Bands = $nbands")

# ============================================================================
# 1. Dense solver (reference)
# ============================================================================
println("\n" * "=" ^ 60)
println("1. Dense Solver (reference)")
println("=" ^ 60)

solver_dense = Solver(PSVWave(), geo, (64, 64); cutoff=cutoff)
println("Matrix size: $(matrix_dimension(solver_dense)) × $(matrix_dimension(solver_dense))")

t_dense = @elapsed bands_dense = compute_bands(solver_dense, kpath; bands=nbands)
println("Time: $(round(t_dense, digits=2)) s")

# ============================================================================
# 2. LOBPCG with warm start (default)
# ============================================================================
println("\n" * "=" ^ 60)
println("2. LOBPCG with Warm Start")
println("   (warm_start=true, scale=true, first_dense=true)")
println("=" ^ 60)

solver_lobpcg = Solver(PSVWave(), geo, (64, 64),
                       LOBPCGMethod(warm_start=true, scale=true, first_dense=true);
                       cutoff=cutoff)

t_lobpcg = @elapsed bands_lobpcg = compute_bands(solver_lobpcg, kpath; bands=nbands)
println("Time: $(round(t_lobpcg, digits=2)) s")

# ============================================================================
# 3. LOBPCG without warm start (for comparison)
# ============================================================================
println("\n" * "=" ^ 60)
println("3. LOBPCG without Warm Start")
println("   (warm_start=false, scale=false, first_dense=false)")
println("=" ^ 60)

solver_lobpcg_no_ws = Solver(PSVWave(), geo, (64, 64),
                              LOBPCGMethod(warm_start=false, scale=false, first_dense=false);
                              cutoff=cutoff)

t_lobpcg_no_ws = @elapsed bands_lobpcg_no_ws = compute_bands(solver_lobpcg_no_ws, kpath; bands=nbands)
println("Time: $(round(t_lobpcg_no_ws, digits=2)) s")

# ============================================================================
# Result comparison
# ============================================================================
println("\n" * "=" ^ 60)
println("Result Comparison")
println("=" ^ 60)

# Dense vs LOBPCG (warm start)
max_error_ws = maximum(abs.(bands_dense.frequencies - bands_lobpcg.frequencies))
rel_error_ws = max_error_ws / maximum(bands_dense.frequencies)

# Dense vs LOBPCG (no warm start)
max_error_no_ws = maximum(abs.(bands_dense.frequencies - bands_lobpcg_no_ws.frequencies))
rel_error_no_ws = max_error_no_ws / maximum(bands_dense.frequencies)

println("\nDense vs LOBPCG (warm start):")
println("  Max absolute error: $(@sprintf("%.2e", max_error_ws)) rad/s")
println("  Relative error: $(@sprintf("%.2e", rel_error_ws))")

println("\nDense vs LOBPCG (no warm start):")
println("  Max absolute error: $(@sprintf("%.2e", max_error_no_ws)) rad/s")
println("  Relative error: $(@sprintf("%.2e", rel_error_no_ws))")

# ============================================================================
# Speedup summary
# ============================================================================
println("\n" * "=" ^ 60)
println("Speedup Summary")
println("=" ^ 60)
println()
println("| Method                  | Time (s) | Speedup |")
println("|-------------------------|----------|---------|")
println("| Dense (reference)       | $(@sprintf("%8.2f", t_dense)) | 1.0x    |")
println("| LOBPCG (warm start)     | $(@sprintf("%8.2f", t_lobpcg)) | $(@sprintf("%.1f", t_dense/t_lobpcg))x    |")
println("| LOBPCG (no warm start)  | $(@sprintf("%8.2f", t_lobpcg_no_ws)) | $(@sprintf("%.1f", t_dense/t_lobpcg_no_ws))x    |")
println()

if t_lobpcg < t_dense
    println("Warm start speedup: $(round(t_dense/t_lobpcg, digits=1))x faster than Dense")
else
    println("Note: For small problems, Dense can be faster due to BLAS optimization")
end

# ============================================================================
# SH Wave for comparison
# ============================================================================
println("\n" * "=" ^ 60)
println("4. SH Wave (Dense)")
println("=" ^ 60)

solver_sh = Solver(SHWave(), geo, (64, 64); cutoff=cutoff)
bands_sh = compute_bands(solver_sh, kpath; bands=1:10)
println("SH bands computed")

# ============================================================================
# Band gap analysis
# ============================================================================
println("\n" * "=" ^ 60)
println("Band Gap Analysis")
println("=" ^ 60)

gaps_psv = find_all_gaps(bands_dense; threshold=0.01)
gaps_sh = find_all_gaps(bands_sh; threshold=0.01)

println("\nP-SV Band Gaps:")
if isempty(gaps_psv)
    println("  No significant band gaps found.")
else
    for g in gaps_psv
        println("  Bands $(g.bands): $(round(g.gap_ratio*100, digits=1))% gap-to-midgap")
    end
end

println("\nSH Band Gaps:")
if isempty(gaps_sh)
    println("  No significant band gaps found.")
else
    for g in gaps_sh
        println("  Bands $(g.bands): $(round(g.gap_ratio*100, digits=1))% gap-to-midgap")
    end
end

# ============================================================================
# Plot band structure
# ============================================================================
println("\n" * "=" ^ 60)
println("Creating plots...")
println("=" ^ 60)

# Normalize frequencies
v_t = sqrt(epoxy.C44 / epoxy.ρ)  # transverse velocity in matrix
norm_factor = a / (2π * v_t)
freqs_psv_norm = bands_dense.frequencies * norm_factor
freqs_sh_norm = bands_sh.frequencies * norm_factor

dists = bands_dense.distances
label_positions = [dists[i] for (i, _) in bands_dense.labels]
label_names = [l for (_, l) in bands_dense.labels]

# Common y-axis range
ymax = max(maximum(freqs_psv_norm), maximum(freqs_sh_norm)) * 1.05
ylims_common = (0, ymax)

# P-SV bands plot
p_psv = plot(
    xlabel = "Wave vector",
    ylabel = "Normalized frequency (ωa/2πvₜ)",
    title = "Vasseur 2001: Steel/Epoxy (P-SV)",
    legend = false,
    grid = true,
    size = (700, 500),
    ylims = ylims_common
)

for b in 1:size(freqs_psv_norm, 2)
    plot!(p_psv, dists, freqs_psv_norm[:, b], linewidth=2, color=:red)
end

# Add P-SV gap shading (normalized)
for g in gaps_psv
    gmin = g.max_lower * norm_factor
    gmax = g.min_upper * norm_factor
    hspan!(p_psv, [gmin, gmax], alpha=0.2, color=:red, label="")
end

vline!(p_psv, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_psv, label_positions, label_names)

# SH bands plot
p_sh = plot(
    xlabel = "Wave vector",
    ylabel = "Normalized frequency (ωa/2πvₜ)",
    title = "Vasseur 2001: Steel/Epoxy (SH)",
    legend = false,
    grid = true,
    size = (700, 500),
    ylims = ylims_common
)

for b in 1:size(freqs_sh_norm, 2)
    plot!(p_sh, dists, freqs_sh_norm[:, b], linewidth=2, color=:blue)
end

# Add SH gap shading (normalized)
for g in gaps_sh
    gmin = g.max_lower * norm_factor
    gmax = g.min_upper * norm_factor
    hspan!(p_sh, [gmin, gmax], alpha=0.2, color=:blue, label="")
end

vline!(p_sh, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_sh, label_positions, label_names)

# Combined plot (side by side)
p_combined = plot(p_sh, p_psv, layout=(1, 2), size=(1200, 500),
    plot_title="Vasseur 2001: Steel/Epoxy Phononic Crystal")

# Save plots
savefig(p_psv, joinpath(@__DIR__, "203_vasseur2001_psv_bands.png"))
savefig(p_sh, joinpath(@__DIR__, "203_vasseur2001_sh_bands.png"))
savefig(p_combined, joinpath(@__DIR__, "203_vasseur2001_combined.png"))

println("\nSaved: 203_vasseur2001_psv_bands.png")
println("Saved: 203_vasseur2001_sh_bands.png")
println("Saved: 203_vasseur2001_combined.png")

display(p_combined)

println("\nDone!")
