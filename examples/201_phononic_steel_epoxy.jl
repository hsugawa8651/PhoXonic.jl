# Last-Modified: 2025-12-15T19:00:00+09:00
# Phononic Crystal Band Structure - P-SV Wave (In-plane elastic)
# Steel cylinders in epoxy matrix
# Reference: Kushwaha et al., Phys. Rev. Lett. 71, 2022 (1993)

using PhoXonic
using Plots
using Printf
using LinearAlgebra: PosDefException

println("=== Phononic Crystal: P-SV Wave ===")
println("Steel cylinders in epoxy matrix")
println()

# Square lattice
a = 1.0  # Lattice constant (normalized)
lat = square_lattice(a)

# Material properties
# Epoxy (matrix)
# ρ = 1180 kg/m³, λ = 4.43 GPa, μ = 1.59 GPa
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)

# Steel (inclusion)
# ρ = 7800 kg/m³, λ = 115 GPa, μ = 82 GPa
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)

println("Materials:")
println("  Epoxy: ρ=$(epoxy.ρ) kg/m³, C11=$(epoxy.C11/1e9) GPa, C44=$(epoxy.C44/1e9) GPa")
println("  Steel: ρ=$(steel.ρ) kg/m³, C11=$(steel.C11/1e9) GPa, C44=$(steel.C44/1e9) GPa")

# Geometry: steel cylinder in epoxy
r = 0.4  # Filling fraction ~ π*0.4² ≈ 0.5
geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], r), steel)])

println("\nGeometry: Steel cylinder, r = $(r)a, filling fraction = $(round(π*r^2, digits=2))")

# Create k-path: Γ → X → M → Γ
kpath = simple_kpath_square(a=a, npoints=40)

# ============================================================================
# SH wave (out-of-plane, scalar) for comparison
# ============================================================================
println("\nComputing SH bands (out-of-plane)...")
solver_sh = Solver(SHWave(), geo, (64, 64); cutoff=7)
bands_sh = compute_bands(solver_sh, kpath; bands=1:8, verbose=false)

println("SH Band Gaps:")
gaps_sh = find_all_gaps(bands_sh; threshold=0.01)
if isempty(gaps_sh)
    println("  No significant band gaps found.")
else
    for g in gaps_sh
        # Convert to frequency in Hz (for a = 1 mm)
        println("  Gap between bands ", g.bands, ": ",
                round(g.gap_ratio*100, digits=1), "% gap-to-midgap")
    end
end

# ============================================================================
# P-SV wave (in-plane, vector)
# ============================================================================
println("\nComputing P-SV bands (in-plane)...")
solver_psv = Solver(PSVWave(), geo, (64, 64); cutoff=7)
bands_psv = compute_bands(solver_psv, kpath; bands=1:12, verbose=false)

println("P-SV Band Gaps:")
gaps_psv = find_all_gaps(bands_psv; threshold=0.01)
if isempty(gaps_psv)
    println("  No significant band gaps found.")
else
    for g in gaps_psv
        println("  Gap between bands ", g.bands, ": ",
                round(g.gap_ratio*100, digits=1), "% gap-to-midgap")
    end
end

# ============================================================================
# Iterative solver comparison (KrylovKit with automatic scaling)
# ============================================================================
println("\n=== Iterative Solver Comparison (KrylovKit) ===")
println("(Phononic eigenvalues ω² ~ 10¹⁰ are automatically scaled)")

println("\nComputing SH bands with KrylovKitMethod...")
solver_sh_krylov = Solver(SHWave(), geo, (64, 64), KrylovKitMethod())
bands_sh_krylov = compute_bands(solver_sh_krylov, kpath; bands=1:8, verbose=false)

max_diff_sh_krylov = maximum(abs.(bands_sh.frequencies .- bands_sh_krylov.frequencies))
rel_diff_sh_krylov = max_diff_sh_krylov / maximum(bands_sh.frequencies)
println("  SH: max |Dense - KrylovKit| = ", @sprintf("%.2e", max_diff_sh_krylov), " rad/s")
println("      relative error = ", @sprintf("%.2e", rel_diff_sh_krylov))

println("\nComputing P-SV bands with KrylovKitMethod...")
solver_psv_krylov = Solver(PSVWave(), geo, (64, 64), KrylovKitMethod())
bands_psv_krylov = compute_bands(solver_psv_krylov, kpath; bands=1:12, verbose=false)

max_diff_psv_krylov = maximum(abs.(bands_psv.frequencies .- bands_psv_krylov.frequencies))
rel_diff_psv_krylov = max_diff_psv_krylov / maximum(bands_psv.frequencies)
println("  PSV: max |Dense - KrylovKit| = ", @sprintf("%.2e", max_diff_psv_krylov), " rad/s")
println("       relative error = ", @sprintf("%.2e", rel_diff_psv_krylov))

# ============================================================================
# Iterative solver comparison (LOBPCG)
# ============================================================================
println("\n=== Iterative Solver Comparison (LOBPCG) ===")
println("(LOBPCG does not require eigenvalue scaling for phononic problems)")

try
    println("\nComputing SH bands with LOBPCGMethod...")
    solver_sh_lobpcg = Solver(SHWave(), geo, (64, 64), LOBPCGMethod())
    bands_sh_lobpcg = compute_bands(solver_sh_lobpcg, kpath; bands=1:8, verbose=false)

    max_diff_sh_lobpcg = maximum(abs.(bands_sh.frequencies .- bands_sh_lobpcg.frequencies))
    rel_diff_sh_lobpcg = max_diff_sh_lobpcg / maximum(bands_sh.frequencies)
    println("  SH: max |Dense - LOBPCG| = ", @sprintf("%.2e", max_diff_sh_lobpcg), " rad/s")
    println("      relative error = ", @sprintf("%.2e", rel_diff_sh_lobpcg))

    println("\nComputing P-SV bands with LOBPCGMethod...")
    solver_psv_lobpcg = Solver(PSVWave(), geo, (64, 64), LOBPCGMethod())
    bands_psv_lobpcg = compute_bands(solver_psv_lobpcg, kpath; bands=1:12, verbose=false)

    max_diff_psv_lobpcg = maximum(abs.(bands_psv.frequencies .- bands_psv_lobpcg.frequencies))
    rel_diff_psv_lobpcg = max_diff_psv_lobpcg / maximum(bands_psv.frequencies)
    println("  PSV: max |Dense - LOBPCG| = ", @sprintf("%.2e", max_diff_psv_lobpcg), " rad/s")
    println("       relative error = ", @sprintf("%.2e", rel_diff_psv_lobpcg))
catch e
    if e isa PosDefException
        println("  LOBPCG skipped: PosDefException (use KrylovKitMethod instead)")
    else
        rethrow(e)
    end
end

# ============================================================================
# Normalize frequencies for plotting
# ============================================================================
# For phononic crystals: ω has units of rad/s
# Normalize: ω̃ = ωa/(2π*v_t) where v_t = √(C44/ρ) is transverse velocity in matrix
v_t_epoxy = sqrt(epoxy.C44 / epoxy.ρ)
println("\nTransverse velocity in epoxy: v_t = $(round(v_t_epoxy, digits=1)) m/s")

# Normalized frequency
norm_factor = a / (2π * v_t_epoxy)

freqs_sh_norm = bands_sh.frequencies * norm_factor
freqs_psv_norm = bands_psv.frequencies * norm_factor

# ============================================================================
# Create plots
# ============================================================================
dists = bands_sh.distances
label_positions = [dists[i] for (i, _) in bands_sh.labels]
label_names = [l for (_, l) in bands_sh.labels]

# Common y-axis range
ymax = max(maximum(freqs_sh_norm), maximum(freqs_psv_norm)) * 1.05
ylims_common = (0, ymax)

# SH bands plot
p_sh = plot(
    xlabel="Wave vector",
    ylabel="Normalized frequency (ωa/2πvₜ)",
    title="SH Bands (out-of-plane)",
    legend=false,
    grid=true,
    size=(600, 450),
    ylims=ylims_common
)
for b in 1:size(freqs_sh_norm, 2)
    plot!(p_sh, dists, freqs_sh_norm[:, b], linewidth=2, color=:blue)
end
vline!(p_sh, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_sh, label_positions, label_names)

# P-SV bands plot
p_psv = plot(
    xlabel="Wave vector",
    ylabel="Normalized frequency (ωa/2πvₜ)",
    title="P-SV Bands (in-plane)",
    legend=false,
    grid=true,
    size=(600, 450),
    ylims=ylims_common
)
for b in 1:size(freqs_psv_norm, 2)
    plot!(p_psv, dists, freqs_psv_norm[:, b], linewidth=2, color=:red)
end
vline!(p_psv, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_psv, label_positions, label_names)

# Combined plot (side by side)
p_combined = plot(p_sh, p_psv, layout=(1, 2), size=(1200, 450),
    plot_title="Steel/Epoxy Phononic Crystal")

# Save plots
savefig(p_sh, joinpath(@__DIR__, "201_phononic_sh_bands.png"))
savefig(p_psv, joinpath(@__DIR__, "201_phononic_psv_bands.png"))
savefig(p_combined, joinpath(@__DIR__, "201_phononic_bands_combined.png"))

println("\nSaved: 201_phononic_sh_bands.png")
println("Saved: 201_phononic_psv_bands.png")
println("Saved: 201_phononic_bands_combined.png")

display(p_combined)
