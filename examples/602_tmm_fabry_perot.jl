# Example 602: Fabry-Pérot Cavity using TMM
#
# Demonstrates resonance behavior of a Fabry-Pérot cavity:
# a dielectric slab between two reflecting surfaces.
# Requires: Plots.jl

using PhoXonic
using Plots

# ============================================================================
# Simple Fabry-Pérot: Dielectric Slab in Air
# ============================================================================

println("=== Simple Fabry-Pérot Cavity ===\n")

# Cavity parameters
n_cavity = 2.0   # Refractive index of cavity
d_cavity = 1.0   # Cavity thickness (normalized units)
n_air = 1.0

mat_cavity = Dielectric(n_cavity^2)
mat_air = Dielectric(n_air^2)

# Create structure: air | cavity | air
ml = Multilayer([Layer(mat_cavity, d_cavity)], mat_air, mat_air)
solver = TMMSolver(Photonic1D(), ml)

# Resonance condition: m*λ = 2*n*d, so λ_m = 2*n*d/m
println("Resonance wavelengths λ_m = 2nd/m:")
for m in 1:5
    λ_res = 2 * n_cavity * d_cavity / m
    println("  m = $m: λ = $(round(λ_res, digits=3))")
end

# Compute spectrum
λ_values = range(0.5, 5.0; length=501)
R, T = tmm_spectrum(solver, collect(λ_values))

# Find transmission peaks
println("\nTransmission peaks (T > 0.99):")
for i in 2:(length(T) - 1)
    if T[i] > T[i - 1] && T[i] > T[i + 1] && T[i] > 0.99
        println("  λ = $(round(λ_values[i], digits=3)), T = $(round(T[i], digits=4))")
    end
end

# ============================================================================
# High-Finesse Cavity with Bragg Mirrors
# ============================================================================

println("\n=== High-Finesse Cavity with Bragg Mirrors ===\n")

# Design wavelength
λ0 = 1.55

# Mirror materials
n_hi, n_lo = 3.0, 1.5
mat_hi = Dielectric(n_hi^2)
mat_lo = Dielectric(n_lo^2)

# Quarter-wave layers for mirrors
d_hi = λ0 / (4 * n_hi)
d_lo = λ0 / (4 * n_lo)

# Mirror unit cells
left_mirror_cell = [Layer(mat_hi, d_hi), Layer(mat_lo, d_lo)]   # (HL)
right_mirror_cell = [Layer(mat_lo, d_lo), Layer(mat_hi, d_hi)]  # (LH)

# Defect: double-thickness high-index layer (full wavelength optical path)
d_defect = 2 * d_hi  # λ/2 physical = λ optical in high-n material

# Build symmetric structure: (HL)^N HH (LH)^N
n_mirror_pairs = 5
left_layers = vcat([l for _ in 1:n_mirror_pairs for l in left_mirror_cell]...)
defect_layer = Layer(mat_hi, d_defect)
right_layers = vcat([l for _ in 1:n_mirror_pairs for l in right_mirror_cell]...)

all_layers = vcat(left_layers, [defect_layer], right_layers)
ml_cavity = Multilayer(all_layers, mat_lo, mat_lo)
solver_cavity = TMMSolver(Photonic1D(), ml_cavity)

println(
    "Structure: $(n_mirror_pairs) mirror pairs | defect (λ/2) | $(n_mirror_pairs) mirror pairs",
)
println("Total layers: $(length(all_layers))")

# Compute spectrum around design wavelength
λ_range = range(0.9*λ0, 1.1*λ0; length=1001)
R_cav, T_cav = tmm_spectrum(solver_cavity, collect(λ_range))

# Find resonance (transmission peak in stopband)
idx_peak = argmax(T_cav)
λ_peak = λ_range[idx_peak]
T_peak = T_cav[idx_peak]

println("\nCavity resonance:")
println("  Peak wavelength: λ = $(round(λ_peak, digits=4))")
println("  Peak transmission: T = $(round(T_peak, digits=4))")
println("  Design wavelength: λ₀ = $λ0")
println("  Deviation: Δλ/λ₀ = $(round(100*(λ_peak - λ0)/λ0, digits=2))%")

# Estimate Q-factor from linewidth
half_max = T_peak / 2
above_half = T_cav .> half_max
if any(above_half)
    first_idx = findfirst(above_half)
    last_idx = findlast(above_half)
    FWHM = λ_range[last_idx] - λ_range[first_idx]
    Q = λ_peak / FWHM
    println("  FWHM: $(round(FWHM, digits=5))")
    println("  Q-factor: $(round(Q, digits=0))")
end

# ============================================================================
# Finesse vs Mirror Reflectivity
# ============================================================================

println("\n=== Finesse vs Number of Mirror Pairs ===\n")

for n_pairs in [3, 5, 7, 10]
    left_layers_n = vcat([l for _ in 1:n_pairs for l in left_mirror_cell]...)
    right_layers_n = vcat([l for _ in 1:n_pairs for l in right_mirror_cell]...)
    all_layers_n = vcat(left_layers_n, [defect_layer], right_layers_n)
    ml_n = Multilayer(all_layers_n, mat_lo, mat_lo)
    solver_n = TMMSolver(Photonic1D(), ml_n)

    R_n, T_n = tmm_spectrum(solver_n, collect(λ_range))

    idx_peak_n = argmax(T_n)
    T_peak_n = T_n[idx_peak_n]
    λ_peak_n = λ_range[idx_peak_n]

    # Estimate FWHM
    half_max_n = T_peak_n / 2
    above_half_n = T_n .> half_max_n
    if any(above_half_n)
        first_idx_n = findfirst(above_half_n)
        last_idx_n = findlast(above_half_n)
        FWHM_n = λ_range[last_idx_n] - λ_range[first_idx_n]
        Q_n = λ_peak_n / FWHM_n
        println(
            "  $n_pairs pairs: T_max = $(round(T_peak_n, digits=3)), Q ≈ $(round(Q_n, digits=0))",
        )
    else
        println("  $n_pairs pairs: T_max = $(round(T_peak_n, digits=3)), Q = N/A")
    end
end

# ============================================================================
# Plot 1: Simple Fabry-Pérot Transmission Spectrum
# ============================================================================

println("\nGenerating plots...")

p1 = plot(;
    xlabel="Wavelength λ",
    ylabel="Transmittance T",
    title="Simple Fabry-Pérot Cavity (n=$n_cavity, d=$d_cavity)",
    legend=false,
    grid=true,
    size=(800, 500),
    ylim=(0, 1.05),
)

plot!(p1, collect(λ_values), T; linewidth=2, color=:blue)

# Mark theoretical resonances
for m in 1:5
    λ_res = 2 * n_cavity * d_cavity / m
    if λ_res >= minimum(λ_values) && λ_res <= maximum(λ_values)
        vline!(p1, [λ_res]; color=:red, linestyle=:dash, alpha=0.5)
    end
end

savefig(p1, joinpath(@__DIR__, "602_fabry_perot_simple.png"))
println("Saved: 602_fabry_perot_simple.png")

# ============================================================================
# Plot 2: High-Finesse Cavity Spectrum
# ============================================================================

p2 = plot(;
    xlabel="Wavelength λ/λ₀",
    ylabel="Transmittance T",
    title="High-Finesse Fabry-Pérot Cavity ($n_mirror_pairs pairs)",
    legend=false,
    grid=true,
    size=(800, 500),
    ylim=(0, 1.05),
)

plot!(p2, collect(λ_range) ./ λ0, T_cav; linewidth=2, color=:blue)

# Mark resonance peak
scatter!(p2, [λ_peak/λ0], [T_peak]; markersize=8, color=:red, label="")

savefig(p2, joinpath(@__DIR__, "602_fabry_perot_cavity.png"))
println("Saved: 602_fabry_perot_cavity.png")

# ============================================================================
# Plot 3: Finesse Comparison
# ============================================================================

p3 = plot(;
    xlabel="Wavelength λ/λ₀",
    ylabel="Transmittance T",
    title="Fabry-Pérot Cavity: Effect of Mirror Pairs",
    legend=:topright,
    grid=true,
    size=(800, 500),
    ylim=(0, 1.05),
)

colors_fp = [:blue, :red, :green, :purple]
for (i, n_pairs) in enumerate([3, 5, 7, 10])
    left_layers_n = vcat([l for _ in 1:n_pairs for l in left_mirror_cell]...)
    right_layers_n = vcat([l for _ in 1:n_pairs for l in right_mirror_cell]...)
    all_layers_n = vcat(left_layers_n, [defect_layer], right_layers_n)
    ml_n = Multilayer(all_layers_n, mat_lo, mat_lo)
    solver_n = TMMSolver(Photonic1D(), ml_n)

    R_n, T_n = tmm_spectrum(solver_n, collect(λ_range))
    plot!(
        p3,
        collect(λ_range) ./ λ0,
        T_n;
        label="$n_pairs pairs",
        linewidth=2,
        color=colors_fp[i],
    )
end

savefig(p3, joinpath(@__DIR__, "602_fabry_perot_finesse.png"))
println("Saved: 602_fabry_perot_finesse.png")

display(p1)
println("\nDone!")
