# LDOS calculation at defect site using RSK
# Shows localized mode detection via LDOS peak
#
# Uses ReducedShiftedKrylov.jl extension

using PhoXonic
using ReducedShiftedKrylov  # Load RSK extension
using Plots

default(
    guidefontsize=14,
    tickfontsize=12,
    titlefontsize=14,
    left_margin=10Plots.mm,
    right_margin=10Plots.mm,
)

println("Example 503: LDOS at Defect Site")
println("=" ^ 50)

# Unit cell geometry
lat = square_lattice(1.0)
geo_unit = Geometry(lat, Dielectric(1.0), [
    (Circle([0.5, 0.5], 0.2), Dielectric(9.0))
])

# 3x3 supercell with center rod removed (defect at position (1,1))
geo = create_supercell(geo_unit, (3, 3); point_defects=[(1, 1)])

solver = Solver(TMWave(), geo, (48, 48); cutoff=5)

println("Structure: 3x3 supercell with center rod removed")
println("Solver: TMWave, 48x48 grid, cutoff=5")
println("Plane waves: $(solver.basis.num_pw)")

# Frequency range around expected defect mode
ω_values = range(0.2, 0.6, length=50)

# LDOS at defect center and away from defect (bulk site)
pos_defect = [1.5, 1.5]  # Center of supercell (defect site)
pos_bulk = [0.5, 0.5]    # Corner site (bulk-like)

# Use Gamma point only for supercell
k_points = [[0.0, 0.0]]

println("\nComputing LDOS at defect site...")
t_defect = @elapsed ldos_defect = compute_ldos(
    solver, pos_defect, ω_values, k_points, MatrixFreeGF();
    η=0.01
)
println("  Time: $(round(t_defect, digits=2)) s")

println("Computing LDOS at bulk site...")
t_bulk = @elapsed ldos_bulk = compute_ldos(
    solver, pos_bulk, ω_values, k_points, MatrixFreeGF();
    η=0.01
)
println("  Time: $(round(t_bulk, digits=2)) s")

# Find peak in defect LDOS
peak_idx = argmax(ldos_defect)
peak_ω = ω_values[peak_idx]

# Plot
p = plot(size=(700, 500))
plot!(p, ω_values, ldos_defect;
    label="Defect site (1.5, 1.5)",
    linewidth=2,
    color=:red,
)
plot!(p, ω_values, ldos_bulk;
    label="Bulk site (0.5, 0.5)",
    linewidth=2,
    linestyle=:dash,
    color=:blue,
)

# Mark the peak
vline!(p, [peak_ω];
    linestyle=:dot,
    color=:gray,
    label="Peak at ω = $(round(peak_ω, digits=3))",
)

plot!(p;
    xlabel="Frequency ω",
    ylabel="LDOS (arb. units)",
    title="LDOS: Defect Mode Detection in 2D Photonic Crystal",
    legend=:topright,
)

output_file = joinpath(@__DIR__, "503_ldos_rsk.png")
savefig(p, output_file)
println("\nSaved: $output_file")

println("\n=== Results ===")
println("Defect mode frequency: ω ≈ $(round(peak_ω, digits=3))")
println("Peak LDOS at defect: $(round(maximum(ldos_defect), digits=2))")
println("Peak LDOS at bulk:   $(round(maximum(ldos_bulk), digits=2))")
println("Enhancement ratio:   $(round(maximum(ldos_defect) / maximum(ldos_bulk), digits=1))x")
