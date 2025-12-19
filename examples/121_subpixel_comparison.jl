# Last-Modified: 2025-12-15T19:00:00+09:00
# Comparison: SimpleGrid vs SubpixelAverage discretization
# Shows convergence improvement with subpixel averaging

using PhoXonic
using Plots
using Statistics

println("=== Discretization Method Comparison ===")
println("Triangular lattice, ε=12 rods in air, r=0.2a")
println()

# Setup geometry
a = 1.0
lat = hexagonal_lattice(a)
air = Dielectric(1.0)
rod = Dielectric(12.0)
geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])

# K-path
kpath = simple_kpath_hexagonal(a=a, npoints=30)

# Reference: high resolution SimpleGrid
println("Computing reference (128x128 SimpleGrid)...")
solver_ref = Solver(TMWave(), geo, (128, 128); cutoff=9)
bands_ref = compute_bands(solver_ref, kpath; bands=1:6, verbose=false)

# Test resolutions
resolutions = [16, 24, 32, 48, 64]

errors_simple = Float64[]
errors_subpix = Float64[]

for res in resolutions
    println("Testing resolution $(res)x$(res)...")

    # SimpleGrid
    solver_s = Solver(TMWave(), geo, (res, res); cutoff=7, discretization=SimpleGrid())
    bands_s = compute_bands(solver_s, kpath; bands=1:6, verbose=false)

    # SubpixelAverage
    solver_p = Solver(TMWave(), geo, (res, res); cutoff=7, discretization=SubpixelAverage(4))
    bands_p = compute_bands(solver_p, kpath; bands=1:6, verbose=false)

    # Calculate RMS error for first band gap edge (band 1 max)
    err_s = sqrt(mean((bands_s.frequencies[:, 1] .- bands_ref.frequencies[:, 1]).^2))
    err_p = sqrt(mean((bands_p.frequencies[:, 1] .- bands_ref.frequencies[:, 1]).^2))

    push!(errors_simple, err_s)
    push!(errors_subpix, err_p)
end

# Print results
println("\n=== Convergence Results ===")
println("Resolution  SimpleGrid   SubpixelAvg  Improvement")
for (i, res) in enumerate(resolutions)
    improvement = errors_simple[i] / errors_subpix[i]
    println("  $(res)x$(res)      $(round(errors_simple[i], digits=5))     $(round(errors_subpix[i], digits=5))      $(round(improvement, digits=1))x")
end

# ============================================================================
# Plot convergence
# ============================================================================
p1 = plot(
    resolutions, errors_simple,
    label="SimpleGrid",
    marker=:circle,
    linewidth=2,
    xlabel="Resolution (N×N)",
    ylabel="RMS Error (band 1)",
    title="Convergence: SimpleGrid vs SubpixelAverage",
    yscale=:log10,
    legend=:topright,
    size=(600, 400)
)
plot!(p1, resolutions, errors_subpix,
      label="SubpixelAverage",
      marker=:square,
      linewidth=2)

savefig(p1, joinpath(@__DIR__, "121_convergence.png"))
println("\nSaved: 121_convergence.png")

# ============================================================================
# Plot band structure comparison at low resolution
# ============================================================================
res_low = 24
solver_s = Solver(TMWave(), geo, (res_low, res_low); cutoff=7, discretization=SimpleGrid())
solver_p = Solver(TMWave(), geo, (res_low, res_low); cutoff=7, discretization=SubpixelAverage(4))

bands_s = compute_bands(solver_s, kpath; bands=1:6, verbose=false)
bands_p = compute_bands(solver_p, kpath; bands=1:6, verbose=false)

dists = bands_ref.distances
label_positions = [dists[i] for (i, _) in bands_ref.labels]
label_names = [l for (_, l) in bands_ref.labels]

# Common y-axis range
ymax = max(maximum(bands_ref.frequencies), maximum(bands_s.frequencies), maximum(bands_p.frequencies)) * 1.05
ylims_common = (0, ymax)

# SimpleGrid plot
p_simple = plot(
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πc)",
    title="SimpleGrid $(res_low)×$(res_low)",
    legend=false,
    grid=true,
    size=(600, 450),
    ylims=ylims_common
)

# Reference (black dashed)
for b in 1:6
    plot!(p_simple, dists, bands_ref.frequencies[:, b], linewidth=1, color=:black, linestyle=:dash)
end
# SimpleGrid (blue)
for b in 1:6
    plot!(p_simple, dists, bands_s.frequencies[:, b], linewidth=2, color=:blue)
end
vline!(p_simple, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_simple, label_positions, label_names)

# SubpixelAverage plot
p_subpix = plot(
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πc)",
    title="SubpixelAverage $(res_low)×$(res_low)",
    legend=false,
    grid=true,
    size=(600, 450),
    ylims=ylims_common
)

# Reference (black dashed)
for b in 1:6
    plot!(p_subpix, dists, bands_ref.frequencies[:, b], linewidth=1, color=:black, linestyle=:dash)
end
# SubpixelAverage (red)
for b in 1:6
    plot!(p_subpix, dists, bands_p.frequencies[:, b], linewidth=2, color=:red)
end
vline!(p_subpix, label_positions, color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_subpix, label_positions, label_names)

# Combined plot (side by side)
p2 = plot(p_simple, p_subpix, layout=(1, 2), size=(1200, 450),
    plot_title="Discretization Comparison (black dashed = 128×128 reference)")

savefig(p2, joinpath(@__DIR__, "121_bands_comparison.png"))
println("Saved: 121_bands_comparison.png")

display(p2)
