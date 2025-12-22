# Last-Modified: 2025-12-15T19:00:00+09:00
# 1D Elastic Superlattice - Longitudinal Wave
# Alternating layers of two materials (e.g., Steel/Epoxy)

using PhoXonic
using Plots
using Printf

println("=== 1D Elastic Superlattice ===")
println("Longitudinal wave in Steel/Epoxy multilayer")
println()

# Material properties
# Steel
ρ_steel = 7800.0   # kg/m³
E_steel = 200e9    # Young's modulus (Pa)
# For 1D longitudinal wave: C11 = E (in thin rod approximation)
# Or use: C11 = λ + 2μ for bulk
C11_steel = E_steel

# Epoxy
ρ_epoxy = 1180.0
E_epoxy = 3.5e9
C11_epoxy = E_epoxy

# Wave velocities
v_steel = sqrt(C11_steel / ρ_steel)
v_epoxy = sqrt(C11_epoxy / ρ_epoxy)

println("Materials:")
println(
    "  Steel: ρ=$(ρ_steel) kg/m³, C11=$(C11_steel/1e9) GPa, v=$(round(v_steel, digits=0)) m/s",
)
println(
    "  Epoxy: ρ=$(ρ_epoxy) kg/m³, C11=$(C11_epoxy/1e9) GPa, v=$(round(v_epoxy, digits=0)) m/s",
)

# Acoustic impedance
Z_steel = ρ_steel * v_steel
Z_epoxy = ρ_epoxy * v_epoxy
println("\nAcoustic impedance:")
println("  Steel: Z=$(round(Z_steel/1e6, digits=2)) MRayl")
println("  Epoxy: Z=$(round(Z_epoxy/1e6, digits=2)) MRayl")
println("  Ratio: $(round(Z_steel/Z_epoxy, digits=1))")

# Layer thicknesses (equal thickness)
d_steel = 0.5  # Normalized to unit cell
d_epoxy = 0.5

# 1D Lattice with period a=1
a = 1.0
lat = lattice_1d(a)

# Geometry: epoxy background, steel layer from 0 to d_steel
steel = IsotropicElastic(ρ_steel, C11_steel, 0.0, 0.0)  # C12=C44=0 for 1D
epoxy = IsotropicElastic(ρ_epoxy, C11_epoxy, 0.0, 0.0)

geo = Geometry(lat, epoxy, [(Segment(0.0, d_steel), steel)])

# Create 1D solver for longitudinal wave
println("\nCreating 1D solver (resolution=128, cutoff=15)...")
solver = Solver(Longitudinal1D(), geo, (128,); cutoff=15)
println("Number of plane waves: ", solver.basis.num_pw)

# K-path: 0 to π/a (first Brillouin zone)
nk = 50
k_points = range(0, π/a; length=nk)

# Compute bands
println("\nComputing bands...")
nbands = 8
frequencies = zeros(nk, nbands)

for (ik, k) in enumerate(k_points)
    freqs, _ = solve(solver, k; bands=1:nbands)
    frequencies[ik, :] = freqs
end

println("Done.")

# ============================================================================
# Iterative solver comparison
# ============================================================================
println("\n=== Iterative Solver Comparison ===")

# KrylovKitMethod
println("\nComputing with KrylovKitMethod...")
solver_krylov = Solver(Longitudinal1D(), geo, 128, KrylovKitMethod(); cutoff=15)
frequencies_krylov = zeros(nk, nbands)
for (ik, k) in enumerate(k_points)
    freqs, _ = solve(solver_krylov, k; bands=1:nbands)
    frequencies_krylov[ik, :] = freqs
end
max_diff_krylov = maximum(abs.(frequencies .- frequencies_krylov))
rel_diff_krylov = max_diff_krylov / maximum(frequencies)
println("  max |Dense - KrylovKit| = ", @sprintf("%.2e", max_diff_krylov), " rad/s")
println("  relative error = ", @sprintf("%.2e", rel_diff_krylov))

# LOBPCGMethod (requires shift for 1D to ensure positive-definite B matrix)
# Note: For phononic problems, shift should be scaled appropriately
# Typical ω² ~ 10^8 for this problem, so shift ~ 10^6 is appropriate
println("\nComputing with LOBPCGMethod (shift=1e6)...")
solver_lobpcg = Solver(Longitudinal1D(), geo, 128, LOBPCGMethod(; shift=1e6); cutoff=15)
frequencies_lobpcg = zeros(nk, nbands)
for (ik, k) in enumerate(k_points)
    freqs, _ = solve(solver_lobpcg, k; bands=1:nbands)
    frequencies_lobpcg[ik, :] = freqs
end
max_diff_lobpcg = maximum(abs.(frequencies .- frequencies_lobpcg))
rel_diff_lobpcg = max_diff_lobpcg / maximum(frequencies)
println("  max |Dense - LOBPCG| = ", @sprintf("%.2e", max_diff_lobpcg), " rad/s")
println("  relative error = ", @sprintf("%.2e", rel_diff_lobpcg))

# ============================================================================
# Find band gaps
# ============================================================================
println("\n=== Band Gaps ===")
for b in 1:(nbands - 1)
    max_lower = maximum(frequencies[:, b])
    min_upper = minimum(frequencies[:, b + 1])
    if min_upper > max_lower
        gap = min_upper - max_lower
        midgap = (max_lower + min_upper) / 2
        gap_ratio = gap / midgap
        println("Gap between bands $b and $(b+1):")
        println(
            "  Range: $(round(max_lower, digits=2)) - $(round(min_upper, digits=2)) rad/s"
        )
        println("  Gap-to-midgap: $(round(gap_ratio * 100, digits=1))%")
    end
end

# ============================================================================
# Normalize frequencies for plotting
# ============================================================================
# Normalize by reference frequency: ω₀ = v_epoxy / a
ω_ref = v_epoxy / a
freq_norm = frequencies ./ ω_ref

# ============================================================================
# Plot band structure
# ============================================================================
p = plot(;
    xlabel="Wave vector (ka/π)",
    ylabel="Normalized frequency (ωa/v_epoxy)",
    title="1D Phononic Crystal: Steel/Epoxy Superlattice",
    legend=false,
    grid=true,
    size=(600, 450),
)

k_norm = k_points ./ (π/a)

for b in 1:nbands
    plot!(p, k_norm, freq_norm[:, b]; linewidth=2, color=:blue)
end

# Highlight band gaps
for b in 1:(nbands - 1)
    max_lower = maximum(freq_norm[:, b])
    min_upper = minimum(freq_norm[:, b + 1])
    if min_upper > max_lower
        hspan!(p, [max_lower, min_upper]; alpha=0.2, color=:yellow, label="")
    end
end

savefig(p, joinpath(@__DIR__, "302_elastic_superlattice_bands.png"))
println("\nSaved: 302_elastic_superlattice_bands.png")

# ============================================================================
# Plot structure (acoustic impedance profile)
# ============================================================================
x_plot = range(0, 2a; length=200)
Z_plot = [x < d_steel || (a <= x < a + d_steel) ? Z_steel : Z_epoxy for x in x_plot]

p_struct = plot(
    x_plot ./ a,
    Z_plot ./ 1e6;
    xlabel="Position (x/a)",
    ylabel="Acoustic impedance (MRayl)",
    title="Steel/Epoxy Superlattice Structure",
    legend=false,
    linewidth=2,
    fill=(0, 0.3, :blue),
    size=(600, 250),
)
vline!(p_struct, [1.0]; color=:gray, linestyle=:dash)

savefig(p_struct, joinpath(@__DIR__, "302_elastic_superlattice_structure.png"))
println("Saved: 302_elastic_superlattice_structure.png")

display(p)
