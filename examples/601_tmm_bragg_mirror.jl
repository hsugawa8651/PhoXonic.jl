# Example 601: Bragg Mirror Transmission Spectrum using TMM
#
# Demonstrates the Transfer Matrix Method for computing transmission
# and reflection spectra of a dielectric Bragg mirror (1D photonic crystal).
# Requires: Plots.jl

using PhoXonic
using Plots

# ============================================================================
# Parameters
# ============================================================================

# Design wavelength
λ0 = 1.55  # 1550 nm (telecom wavelength), normalized units

# Materials: high/low refractive index
n_hi = 2.5   # e.g., TiO2
n_lo = 1.45  # e.g., SiO2

mat_hi = Dielectric(n_hi^2)
mat_lo = Dielectric(n_lo^2)
mat_air = Dielectric(1.0)

# Quarter-wave thicknesses at design wavelength
d_hi = λ0 / (4 * n_hi)
d_lo = λ0 / (4 * n_lo)

println("Bragg Mirror Design:")
println("  Design wavelength: λ₀ = $λ0")
println("  High-index layer: n = $n_hi, d = $(round(d_hi, digits=4))")
println("  Low-index layer:  n = $n_lo, d = $(round(d_lo, digits=4))")

# ============================================================================
# Build Multilayer Structure
# ============================================================================

# Unit cell: [high, low]
unit_cell = [Layer(mat_hi, d_hi), Layer(mat_lo, d_lo)]

# Create Bragg mirrors with different number of pairs
n_pairs_list = [5, 10, 20]
solvers = Dict{Int,TMMSolver}()

for n_pairs in n_pairs_list
    ml = periodic_multilayer(unit_cell, n_pairs)
    solvers[n_pairs] = TMMSolver(Photonic1D(), ml)
end

# ============================================================================
# Compute Transmission Spectrum
# ============================================================================

# Wavelength range
λ_min, λ_max = 0.8 * λ0, 1.2 * λ0
λ_values = range(λ_min, λ_max; length=201)

println("\nComputing transmission spectra...")

results = Dict{Int,Tuple{Vector{Float64},Vector{Float64}}}()
for n_pairs in n_pairs_list
    R, T = tmm_spectrum(solvers[n_pairs], collect(λ_values))
    results[n_pairs] = (R, T)

    # Find reflectivity at design wavelength
    idx_λ0 = argmin(abs.(collect(λ_values) .- λ0))
    R_λ0 = R[idx_λ0]
    println("  $n_pairs pairs: R(λ₀) = $(round(100*R_λ0, digits=2))%")
end

# ============================================================================
# Stopband Analysis
# ============================================================================

println("\nStopband Analysis (R > 99%):")

# Theoretical stopband width
Δω_theory = (4/π) * asin((n_hi - n_lo) / (n_hi + n_lo))
println("  Theoretical Δω/ω₀ = $(round(100*Δω_theory, digits=1))%")

for n_pairs in n_pairs_list
    R, _ = results[n_pairs]

    # Find stopband edges (R > 0.99)
    in_stopband = R .> 0.99
    if any(in_stopband)
        first_idx = findfirst(in_stopband)
        last_idx = findlast(in_stopband)
        λ_start = λ_values[first_idx]
        λ_stop = λ_values[last_idx]
        Δλ = λ_stop - λ_start
        println(
            "  $n_pairs pairs: stopband $(round(λ_start, digits=3)) - $(round(λ_stop, digits=3)), Δλ/λ₀ = $(round(100*Δλ/λ0, digits=1))%",
        )
    else
        println("  $n_pairs pairs: no R > 99% region")
    end
end

# ============================================================================
# Oblique Incidence (Optional)
# ============================================================================

println("\nOblique Incidence (20 pairs, λ = λ₀):")

solver_20 = solvers[20]
angles = [0, 15, 30, 45, 60]

for θ_deg in angles
    θ = deg2rad(θ_deg)
    result_te = tmm_spectrum(solver_20, λ0; angle=θ, polarization=:TE)
    result_tm = tmm_spectrum(solver_20, λ0; angle=θ, polarization=:TM)

    println(
        "  θ = $(θ_deg)°: R_TE = $(round(100*result_te.R, digits=1))%, R_TM = $(round(100*result_tm.R, digits=1))%",
    )
end

# ============================================================================
# Plot 1: Transmission Spectrum for Different Numbers of Pairs
# ============================================================================

println("\nGenerating plots...")

p1 = plot(;
    xlabel="Wavelength λ/λ₀",
    ylabel="Transmittance T",
    title="Bragg Mirror: Transmittance vs Number of Pairs",
    legend=:topright,
    grid=true,
    size=(800, 500),
    ylim=(0, 1.05),
)

colors = [:blue, :red, :green]
for (i, n_pairs) in enumerate(n_pairs_list)
    R, T = results[n_pairs]
    plot!(
        p1,
        collect(λ_values) ./ λ0,
        T;
        label="$(n_pairs) pairs",
        linewidth=2,
        color=colors[i],
    )
end

# Mark design wavelength
vline!(p1, [1.0]; color=:gray, linestyle=:dash, label="λ₀", alpha=0.7)

savefig(p1, joinpath(@__DIR__, "601_bragg_transmittance.png"))
println("Saved: 601_bragg_transmittance.png")

# ============================================================================
# Plot 2: Reflectance Spectrum (20 pairs)
# ============================================================================

R_20, T_20 = results[20]

p2 = plot(;
    xlabel="Wavelength λ/λ₀",
    ylabel="Reflectance / Transmittance",
    title="Bragg Mirror (20 pairs): R and T Spectra",
    legend=:right,
    grid=true,
    size=(800, 500),
    ylim=(0, 1.05),
)

plot!(p2, collect(λ_values) ./ λ0, R_20; label="R", linewidth=2, color=:red)
plot!(p2, collect(λ_values) ./ λ0, T_20; label="T", linewidth=2, color=:blue)

# Highlight stopband
in_stopband = R_20 .> 0.99
if any(in_stopband)
    first_idx = findfirst(in_stopband)
    last_idx = findlast(in_stopband)
    λ_start = λ_values[first_idx] / λ0
    λ_stop = λ_values[last_idx] / λ0
    vspan!(p2, [λ_start, λ_stop]; alpha=0.15, color=:yellow, label="Stopband")
end

savefig(p2, joinpath(@__DIR__, "601_bragg_spectrum.png"))
println("Saved: 601_bragg_spectrum.png")

# ============================================================================
# Plot 3: Oblique Incidence - Angular Dependence
# ============================================================================

angles_fine = 0:2:80
R_te = Float64[]
R_tm = Float64[]

for θ_deg in angles_fine
    θ = deg2rad(θ_deg)
    result_te = tmm_spectrum(solver_20, λ0; angle=θ, polarization=:TE)
    result_tm = tmm_spectrum(solver_20, λ0; angle=θ, polarization=:TM)
    push!(R_te, result_te.R)
    push!(R_tm, result_tm.R)
end

p3 = plot(;
    xlabel="Incident Angle θ (degrees)",
    ylabel="Reflectance at λ₀",
    title="Bragg Mirror (20 pairs): Angular Dependence",
    legend=:bottomleft,
    grid=true,
    size=(800, 500),
    ylim=(0, 1.05),
)

plot!(p3, collect(angles_fine), R_te; label="TE (s-pol)", linewidth=2, color=:blue)
plot!(p3, collect(angles_fine), R_tm; label="TM (p-pol)", linewidth=2, color=:red)

savefig(p3, joinpath(@__DIR__, "601_bragg_angular.png"))
println("Saved: 601_bragg_angular.png")

display(p1)
println("\nDone!")
