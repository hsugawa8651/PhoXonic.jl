# Last-Modified: 2025-12-19
# Kushwaha 1993: Phononic crystal band structure
# Reference: M.S. Kushwaha et al., Phys. Rev. Lett. 71, 2022 (1993)
# DOI: https://doi.org/10.1103/PhysRevLett.71.2022
#
# Demonstrates:
# 1. Dense vs LOBPCG result consistency
# 2. LOBPCG warm start speedup
# 3. SH and P-SV wave comparison

using PhoXonic
using Printf
using Plots

println("=" ^ 60)
println("Kushwaha 1993: Phononic Crystal")
println("Square lattice, Steel cylinders in Epoxy matrix")
println("=" ^ 60)
println()

# ============================================================================
# Material properties (typical steel/epoxy values)
# ============================================================================
# Steel (inclusion)
ρ_steel = 7800.0      # kg/m³
λ_steel = 115e9       # Pa
μ_steel = 82e9        # Pa

# Epoxy (matrix)
ρ_epoxy = 1180.0      # kg/m³
λ_epoxy = 4.43e9      # Pa
μ_epoxy = 1.59e9      # Pa

steel = IsotropicElastic(ρ=ρ_steel, λ=λ_steel, μ=μ_steel)
epoxy = IsotropicElastic(ρ=ρ_epoxy, λ=λ_epoxy, μ=μ_epoxy)

println("Materials:")
println("  Steel: ρ = $ρ_steel kg/m³, C44 = $(steel.C44/1e9) GPa")
println("  Epoxy: ρ = $ρ_epoxy kg/m³, C44 = $(epoxy.C44/1e9) GPa")

# ============================================================================
# Geometry: square lattice with steel cylinders
# ============================================================================
a = 1.0  # lattice constant (normalized)
lat = square_lattice(a)

r = 0.4  # cylinder radius, filling fraction ~ π*0.4² ≈ 0.50
geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], r), steel)])

println("\nGeometry:")
println("  Square lattice, a = $a")
println("  Cylinder radius r = $r")
println("  Filling fraction f = $(round(π*r^2, digits=2))")

# ============================================================================
# K-path: Γ → X → M → Γ
# ============================================================================
npoints = 40
kpath = simple_kpath_square(a=a, npoints=npoints)
println("\nK-path: $(length(kpath.points)) points")

# ============================================================================
# Solver parameters
# ============================================================================
cutoff = 10
nbands_sh = 1:8
nbands_psv = 1:12

println("\nSolver parameters:")
println("  Cutoff = $cutoff")
println("  SH bands = $nbands_sh")
println("  P-SV bands = $nbands_psv")

# ============================================================================
# SH Wave (scalar, out-of-plane)
# ============================================================================
println("\n" * "=" ^ 60)
println("SH Wave (out-of-plane, scalar)")
println("=" ^ 60)

# Dense
solver_sh_dense = Solver(SHWave(), geo, (64, 64); cutoff=cutoff)
println("Matrix size: $(matrix_dimension(solver_sh_dense))")

t_sh_dense = @elapsed bands_sh_dense = compute_bands(solver_sh_dense, kpath; bands=nbands_sh)
println("Dense: $(round(t_sh_dense, digits=2)) s")

# LOBPCG with warm start
solver_sh_lobpcg = Solver(SHWave(), geo, (64, 64),
                          LOBPCGMethod(warm_start=true); cutoff=cutoff)

t_sh_lobpcg = @elapsed bands_sh_lobpcg = compute_bands(solver_sh_lobpcg, kpath; bands=nbands_sh)
println("LOBPCG (warm start): $(round(t_sh_lobpcg, digits=2)) s")

# Error check
max_error_sh = maximum(abs.(bands_sh_dense.frequencies - bands_sh_lobpcg.frequencies))
rel_error_sh = max_error_sh / maximum(bands_sh_dense.frequencies)
println("Max error: $(@sprintf("%.2e", max_error_sh)) rad/s (rel: $(@sprintf("%.2e", rel_error_sh)))")

# ============================================================================
# P-SV Wave (vector, in-plane)
# ============================================================================
println("\n" * "=" ^ 60)
println("P-SV Wave (in-plane, vector)")
println("=" ^ 60)

# Dense
solver_psv_dense = Solver(PSVWave(), geo, (64, 64); cutoff=cutoff)
println("Matrix size: $(matrix_dimension(solver_psv_dense))")

t_psv_dense = @elapsed bands_psv_dense = compute_bands(solver_psv_dense, kpath; bands=nbands_psv)
println("Dense: $(round(t_psv_dense, digits=2)) s")

# LOBPCG with warm start
solver_psv_lobpcg = Solver(PSVWave(), geo, (64, 64),
                           LOBPCGMethod(warm_start=true); cutoff=cutoff)

t_psv_lobpcg = @elapsed bands_psv_lobpcg = compute_bands(solver_psv_lobpcg, kpath; bands=nbands_psv)
println("LOBPCG (warm start): $(round(t_psv_lobpcg, digits=2)) s")

# Error check
max_error_psv = maximum(abs.(bands_psv_dense.frequencies - bands_psv_lobpcg.frequencies))
rel_error_psv = max_error_psv / maximum(bands_psv_dense.frequencies)
println("Max error: $(@sprintf("%.2e", max_error_psv)) rad/s (rel: $(@sprintf("%.2e", rel_error_psv)))")

# ============================================================================
# Summary
# ============================================================================
println("\n" * "=" ^ 60)
println("Summary")
println("=" ^ 60)
println()
println("| Wave | Method | Time (s) | Speedup | Max Error |")
println("|------|--------|----------|---------|-----------|")
println("| SH   | Dense  | $(@sprintf("%8.2f", t_sh_dense)) | 1.0x    | (ref)     |")
println("| SH   | LOBPCG | $(@sprintf("%8.2f", t_sh_lobpcg)) | $(@sprintf("%.1f", t_sh_dense/t_sh_lobpcg))x    | $(@sprintf("%.1e", max_error_sh)) |")
println("| PSV  | Dense  | $(@sprintf("%8.2f", t_psv_dense)) | 1.0x    | (ref)     |")
println("| PSV  | LOBPCG | $(@sprintf("%8.2f", t_psv_lobpcg)) | $(@sprintf("%.1f", t_psv_dense/t_psv_lobpcg))x    | $(@sprintf("%.1e", max_error_psv)) |")
println()

# ============================================================================
# Plot band structure
# ============================================================================
println("Creating plots...")

# Normalize frequencies
v_t = sqrt(epoxy.C44 / epoxy.ρ)  # transverse velocity in matrix
norm_factor = a / (2π * v_t)

freqs_sh_norm = bands_sh_dense.frequencies * norm_factor
freqs_psv_norm = bands_psv_dense.frequencies * norm_factor

dists = bands_sh_dense.distances
label_positions = [dists[i] for (i, _) in bands_sh_dense.labels]
label_names = [l for (_, l) in bands_sh_dense.labels]

# Combined plot
p = plot(
    xlabel = "Wave vector",
    ylabel = "Normalized frequency (ωa/2πvₜ)",
    title = "Kushwaha 1993: Steel/Epoxy\nSH (blue) and P-SV (red)",
    legend = false,
    grid = true,
    size = (700, 500)
)

for b in 1:size(freqs_sh_norm, 2)
    plot!(p, dists, freqs_sh_norm[:, b], linewidth=2, color=:blue)
end
for b in 1:size(freqs_psv_norm, 2)
    plot!(p, dists, freqs_psv_norm[:, b], linewidth=2, color=:red, linestyle=:dash)
end

vline!(p, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p, label_positions, label_names)

savefig(p, joinpath(@__DIR__, "205_kushwaha1993_bands.png"))
println("\nSaved: 205_kushwaha1993_bands.png")

display(p)

println("\nDone!")
