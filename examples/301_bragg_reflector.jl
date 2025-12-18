# Last-Modified: 2025-12-15T19:00:00+09:00
# 1D Bragg Reflector (Quarter-Wave Stack)
# Equivalent to MPB's bragg.ctl
# Alternating layers of high/low refractive index materials

using PhoXonic
using Plots
using Printf

println("=== 1D Bragg Reflector ===")
println("Quarter-wave stack: n_hi=3.0, n_lo=1.0")
println()

# Parameters (from MPB bragg.ctl)
n_lo = 1.0
n_hi = 3.0
ε_lo = n_lo^2  # = 1.0
ε_hi = n_hi^2  # = 9.0

# Quarter-wave condition: n_hi * d_hi = n_lo * d_lo = λ/4
# For unit cell: d_hi + d_lo = 1 (normalized)
# d_hi = n_lo / (n_hi + n_lo) = 0.25
# d_lo = n_hi / (n_hi + n_lo) = 0.75
d_hi = n_lo / (n_hi + n_lo)  # 0.25
d_lo = n_hi / (n_hi + n_lo)  # 0.75

println("Layer thicknesses: d_hi = $(d_hi), d_lo = $(d_lo)")
println("Dielectric constants: ε_hi = $(ε_hi), ε_lo = $(ε_lo)")

# 1D Lattice with period a=1
a = 1.0
lat = lattice_1d(a)

# Geometry: low-index background, high-index layer at center
# High-index layer from 0 to d_hi (centered at d_hi/2)
background = Dielectric(ε_lo)
layer_hi = Dielectric(ε_hi)

# Segment from 0 to d_hi
geo = Geometry(lat, background, [(Segment(0.0, d_hi), layer_hi)])

# Create 1D solver
println("\nCreating 1D solver (resolution=128, cutoff=15)...")
solver = Solver(Photonic1D(), geo, (128,); cutoff=15)
println("Number of plane waves: ", solver.basis.num_pw)

# K-path: 0 to π/a (first Brillouin zone)
nk = 50
k_points = range(0, π/a, length=nk)

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
solver_krylov = Solver(Photonic1D(), geo, 128, KrylovKitMethod(); cutoff=15)
frequencies_krylov = zeros(nk, nbands)
for (ik, k) in enumerate(k_points)
    freqs, _ = solve(solver_krylov, k; bands=1:nbands)
    frequencies_krylov[ik, :] = freqs
end
max_diff_krylov = maximum(abs.(frequencies .- frequencies_krylov))
println("  max |Dense - KrylovKit| = ", @sprintf("%.2e", max_diff_krylov))

# LOBPCGMethod (requires shift for 1D to ensure positive-definite B matrix)
println("\nComputing with LOBPCGMethod (shift=0.01)...")
solver_lobpcg = Solver(Photonic1D(), geo, 128, LOBPCGMethod(shift=0.01); cutoff=15)
frequencies_lobpcg = zeros(nk, nbands)
for (ik, k) in enumerate(k_points)
    freqs, _ = solve(solver_lobpcg, k; bands=1:nbands)
    frequencies_lobpcg[ik, :] = freqs
end
max_diff_lobpcg = maximum(abs.(frequencies .- frequencies_lobpcg))
println("  max |Dense - LOBPCG| = ", @sprintf("%.2e", max_diff_lobpcg))

# ============================================================================
# Analytical band gap location
# ============================================================================
# For quarter-wave stack, fundamental gap is at ω₀ = π c / (2 n_eff a)
# where n_eff ≈ (n_hi + n_lo) / 2
# Gap-to-midgap ratio ≈ 4/π * arcsin((n_hi - n_lo)/(n_hi + n_lo))

n_eff = (n_hi + n_lo) / 2
ω_center = π / (2 * n_eff * a) * (2π)  # Normalized: ωa/2πc
gap_ratio_theory = 4/π * asin((n_hi - n_lo)/(n_hi + n_lo))

println("\n=== Theoretical Values ===")
println("Center frequency (ωa/2πc): ", round(ω_center / (2π), digits=3))
println("Gap-to-midgap ratio: ", round(gap_ratio_theory * 100, digits=1), "%")

# ============================================================================
# Find band gaps from computed data
# ============================================================================
println("\n=== Computed Band Gaps ===")
for b in 1:(nbands-1)
    max_lower = maximum(frequencies[:, b])
    min_upper = minimum(frequencies[:, b+1])
    if min_upper > max_lower
        gap = min_upper - max_lower
        midgap = (max_lower + min_upper) / 2
        gap_ratio = gap / midgap
        println("Gap between bands $b and $(b+1):")
        println("  Range: $(round(max_lower, digits=4)) - $(round(min_upper, digits=4))")
        println("  Gap-to-midgap: $(round(gap_ratio * 100, digits=1))%")
    end
end

# ============================================================================
# Plot
# ============================================================================
p = plot(
    xlabel="Wave vector (ka/π)",
    ylabel="Frequency (ωa/2πc)",
    title="1D Bragg Reflector (n=1/3 quarter-wave stack)",
    legend=false,
    grid=true,
    size=(600, 450)
)

k_norm = k_points ./ (π/a)  # Normalize to [0, 1]

for b in 1:nbands
    plot!(p, k_norm, frequencies[:, b], linewidth=2, color=:blue)
end

# Add light line for comparison
# plot!(p, k_norm, k_points, linewidth=1, color=:gray, linestyle=:dash)

# Highlight band gaps
for b in 1:(nbands-1)
    max_lower = maximum(frequencies[:, b])
    min_upper = minimum(frequencies[:, b+1])
    if min_upper > max_lower
        hspan!(p, [max_lower, min_upper], alpha=0.2, color=:yellow, label="")
    end
end

savefig(p, joinpath(@__DIR__, "301_bragg_bands.png"))
println("\nSaved: 301_bragg_bands.png")

# ============================================================================
# Plot structure
# ============================================================================
x_plot = range(0, 2a, length=200)
ε_plot = [x < d_hi || (a <= x < a + d_hi) ? ε_hi : ε_lo for x in x_plot]

p_struct = plot(
    x_plot ./ a, ε_plot,
    xlabel="Position (x/a)",
    ylabel="Dielectric constant ε",
    title="Bragg Reflector Structure",
    legend=false,
    linewidth=2,
    fill=(0, 0.3, :blue),
    size=(600, 250)
)
vline!(p_struct, [1.0], color=:gray, linestyle=:dash)

savefig(p_struct, joinpath(@__DIR__, "301_bragg_structure.png"))
println("Saved: 301_bragg_structure.png")

display(p)
