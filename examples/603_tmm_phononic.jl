# Example 603: Phononic Crystal TMM - Steel/Epoxy Superlattice
#
# Demonstrates TMM for longitudinal elastic waves in a periodic
# steel/epoxy multilayer structure.
# Requires: Plots.jl

using PhoXonic
using Plots
default(guidefontsize=14, tickfontsize=12, titlefontsize=14, left_margin=10Plots.mm, right_margin=10Plots.mm, top_margin=5Plots.mm, bottom_margin=10Plots.mm)

# ============================================================================
# Material Properties
# ============================================================================

# Steel: high impedance
steel = IsotropicElastic(; ρ=7800.0, λ=115e9, μ=82e9)

# Epoxy: low impedance
epoxy = IsotropicElastic(; ρ=1180.0, λ=4.43e9, μ=1.59e9)

# Material parameters
c_steel = longitudinal_velocity(steel)
c_epoxy = longitudinal_velocity(epoxy)
Z_steel = PhoXonic.acoustic_impedance(steel)
Z_epoxy = PhoXonic.acoustic_impedance(epoxy)

println("=== Material Properties ===\n")
println("Steel:")
println("  Density: $(PhoXonic.density(steel)) kg/m³")
println("  Longitudinal velocity: $(round(c_steel, digits=1)) m/s")
println("  Acoustic impedance: $(round(Z_steel/1e6, digits=2)) MRayl")

println("\nEpoxy:")
println("  Density: $(PhoXonic.density(epoxy)) kg/m³")
println("  Longitudinal velocity: $(round(c_epoxy, digits=1)) m/s")
println("  Acoustic impedance: $(round(Z_epoxy/1e6, digits=2)) MRayl")

println("\nImpedance ratio Z_steel/Z_epoxy: $(round(Z_steel/Z_epoxy, digits=1))")

# ============================================================================
# Single Interface Reflection
# ============================================================================

println("\n=== Single Interface Reflection ===\n")

r, t = PhoXonic.acoustic_fresnel(Z_steel, Z_epoxy)
R = abs2(r)
T_coeff = 4 * Z_steel * Z_epoxy / (Z_steel + Z_epoxy)^2

println("Steel → Epoxy interface:")
println("  Reflection coefficient r = $(round(r, digits=3))")
println("  Transmission coefficient t = $(round(t, digits=3))")
println("  Reflectivity R = $(round(100*R, digits=1))%")
println("  Transmissivity T = $(round(100*T_coeff, digits=1))%")

# ============================================================================
# Phononic Bragg Mirror
# ============================================================================

println("\n=== Phononic Bragg Mirror ===\n")

# Design frequency
f0 = 100e3  # 100 kHz
λ_steel = c_steel / f0
λ_epoxy = c_epoxy / f0

# Quarter-wave thicknesses
d_steel = λ_steel / 4
d_epoxy = λ_epoxy / 4

println("Design frequency: $(f0/1e3) kHz")
println("Quarter-wave thicknesses:")
println("  Steel: $(round(d_steel*1e3, digits=2)) mm")
println("  Epoxy: $(round(d_epoxy*1e3, digits=2)) mm")

# Unit cell
unit_cell = [Layer(steel, d_steel), Layer(epoxy, d_epoxy)]

# Create mirrors with different numbers of pairs
println("\nReflectivity at design frequency:")
for n_pairs in [3, 5, 10, 15]
    ml = periodic_multilayer(unit_cell, n_pairs)
    solver = TMMSolver(Longitudinal1D(), ml)

    result = tmm_spectrum(solver, λ_steel)
    println("  $n_pairs pairs: R = $(round(100*result.R, digits=2))%")
end

# ============================================================================
# Transmission Spectrum
# ============================================================================

println("\n=== Transmission Spectrum ===\n")

n_pairs = 10
ml = periodic_multilayer(unit_cell, n_pairs)
solver = TMMSolver(Longitudinal1D(), ml)

# Frequency range (50 kHz to 200 kHz)
f_values = range(50e3, 200e3; length=151)
λ_values = c_steel ./ f_values

R_spectrum = Float64[]
T_spectrum = Float64[]

for λ in λ_values
    result = tmm_spectrum(solver, λ)
    push!(R_spectrum, result.R)
    push!(T_spectrum, result.T)
end

# Find stopband (R > 90%)
in_stopband = R_spectrum .> 0.90
if any(in_stopband)
    first_idx = findfirst(in_stopband)
    last_idx = findlast(in_stopband)
    f_start = f_values[first_idx]
    f_stop = f_values[last_idx]
    println(
        "Stopband (R > 90%): $(round(f_start/1e3, digits=1)) - $(round(f_stop/1e3, digits=1)) kHz",
    )
    println("Stopband width: $(round((f_stop - f_start)/1e3, digits=1)) kHz")
    println("Relative width: $(round(100*(f_stop - f_start)/f0, digits=1))%")
end

# ============================================================================
# Band Structure
# ============================================================================

println("\n=== Band Structure ===\n")

# Single period for band structure
ml_period = Multilayer(unit_cell, epoxy, epoxy)
solver_period = TMMSolver(Longitudinal1D(), ml_period)

# Compute bands
bands = tmm_bandstructure(solver_period; k_points=51, bands=1:4)

# Period
a = d_steel + d_epoxy

println("Period: $(round(a*1e3, digits=2)) mm")
println("First Brillouin zone: k ∈ [0, π/a] = [0, $(round(π/a, digits=1))] rad/m")

# Report band edges
println("\nBand edges (at zone boundary k = π/a):")
for i in 1:4
    ω_edge = bands.frequencies[end, i]
    f_edge = ω_edge / (2π)
    println(
        "  Band $i: ω = $(round(ω_edge, digits=1)) rad/s = $(round(f_edge/1e3, digits=1)) kHz",
    )
end

# Check for bandgaps
println("\nBandgaps:")
for i in 1:3
    band_i_max = maximum(bands.frequencies[:, i])
    band_i1_min = minimum(bands.frequencies[:, i + 1])

    if band_i1_min > band_i_max
        gap_width = band_i1_min - band_i_max
        gap_center = (band_i_max + band_i1_min) / 2
        f_center = gap_center / (2π)
        f_width = gap_width / (2π)
        println("  Gap between bands $i-$(i+1):")
        println("    Center: $(round(f_center/1e3, digits=1)) kHz")
        println("    Width: $(round(f_width/1e3, digits=1)) kHz")
        println("    Relative: $(round(100*gap_width/gap_center, digits=1))%")
    end
end

# ============================================================================
# Plot 1: Transmission Spectrum
# ============================================================================

println("\nGenerating plots...")

p1 = plot(;
    xlabel="Frequency (kHz)",
    ylabel="Reflectance / Transmittance",
    title="Phononic Bragg Mirror: Steel/Epoxy ($n_pairs pairs)",
    legend=:right,
    grid=true,
    size=(800, 500),
    ylim=(0, 1.05),
)

plot!(p1, collect(f_values) ./ 1e3, R_spectrum; label="R", linewidth=2, color=:red)
plot!(p1, collect(f_values) ./ 1e3, T_spectrum; label="T", linewidth=2, color=:blue)

# Mark design frequency
vline!(p1, [f0/1e3]; color=:gray, linestyle=:dash, label="f₀", alpha=0.7)

# Highlight stopband
in_stopband = R_spectrum .> 0.90
if any(in_stopband)
    first_idx = findfirst(in_stopband)
    last_idx = findlast(in_stopband)
    f_start = f_values[first_idx] / 1e3
    f_stop = f_values[last_idx] / 1e3
    vspan!(p1, [f_start, f_stop]; alpha=0.15, color=:yellow, label="Stopband")
end

savefig(p1, joinpath(@__DIR__, "603_phononic_spectrum.png"))
println("Saved: 603_phononic_spectrum.png")

# ============================================================================
# Plot 2: Band Structure
# ============================================================================

p2 = plot(;
    xlabel="Wave vector k (rad/m)",
    ylabel="Frequency (kHz)",
    title="Phononic Band Structure: Steel/Epoxy",
    legend=false,
    grid=true,
    size=(800, 500),
)

k_values = bands.distances

for b in 1:size(bands.frequencies, 2)
    # Convert ω to frequency in kHz
    freqs_kHz = bands.frequencies[:, b] ./ (2π * 1e3)
    plot!(p2, k_values, freqs_kHz; linewidth=2, color=:blue)
end

# Highlight band gaps
for i in 1:3
    band_i_max = maximum(bands.frequencies[:, i])
    band_i1_min = minimum(bands.frequencies[:, i + 1])
    if band_i1_min > band_i_max
        f_lower = band_i_max / (2π * 1e3)
        f_upper = band_i1_min / (2π * 1e3)
        hspan!(p2, [f_lower, f_upper]; alpha=0.2, color=:yellow)
    end
end

# Add zone boundary labels
xticks!(p2, [0, π/a], ["Γ", "X"])

savefig(p2, joinpath(@__DIR__, "603_phononic_bands.png"))
println("Saved: 603_phononic_bands.png")

# ============================================================================
# Plot 3: Reflectivity vs Number of Pairs
# ============================================================================

p3 = plot(;
    xlabel="Number of Pairs",
    ylabel="Reflectivity at f₀ (%)",
    title="Phononic Bragg Mirror: Reflectivity vs Pairs",
    legend=false,
    grid=true,
    size=(800, 500),
    ylim=(0, 105),
)

pairs_list = 1:20
R_vs_pairs = Float64[]

for np in pairs_list
    ml_np = periodic_multilayer(unit_cell, np)
    solver_np = TMMSolver(Longitudinal1D(), ml_np)
    result_np = tmm_spectrum(solver_np, λ_steel)
    push!(R_vs_pairs, 100 * result_np.R)
end

plot!(
    p3,
    collect(pairs_list),
    R_vs_pairs;
    linewidth=2,
    color=:blue,
    marker=:circle,
    markersize=4,
)

savefig(p3, joinpath(@__DIR__, "603_phononic_pairs.png"))
println("Saved: 603_phononic_pairs.png")

display(p1)
println("\nDone!")
