# Example 811: 1D Field Visualization
#
# Demonstrates field reconstruction and visualization for 1D photonic crystals.
# Output: 811_field_1d.png

using PhoXonic
using Plots

println("Example 811: 1D Field Visualization")
println("=" ^ 50)

# Create 1D photonic crystal: air + dielectric slab
lat = lattice_1d(1.0)
geo = Geometry(lat, Dielectric(1.0), [
    (Segment(0.0, 0.3), Dielectric(9.0))  # 30% filling fraction
])

# Create solver
solver = Solver(Photonic1D(), geo, (128,); cutoff=10)
println("Resolution: 128, Cutoff: 10")
println("Number of plane waves: $(solver.basis.num_pw)")

# Solve at Gamma point (k = 0)
freqs, vecs = solve_at_k_with_vectors(solver, [0.0], DenseMethod())
println("\nFirst 4 frequencies: ", round.(freqs[1:4], digits=4))

# Create multi-panel plot
p = plot(layout=(2, 2), size=(1000, 800),
    tickfontsize=12, guidefontsize=14, titlefontsize=14,
    left_margin=5Plots.mm, bottom_margin=5Plots.mm, top_margin=-2Plots.mm,
    legend=:topright)

# Get epsilon for reference
eps = get_epsilon_field(solver)
x = range(0, 1, length=128)

# Plot first 4 modes
for i in 1:4
    field = reconstruct_field(solver, vecs[:, i])
    field_fixed = fix_phase(field)

    plot!(p, x, real.(field_fixed);
        subplot=i,
        title="Band $i (ω = $(round(freqs[i], digits=3)))",
        xlabel="x/a",
        ylabel="Field",
        label="Field",
        linewidth=2
    )

    # Overlay scaled epsilon as filled area
    eps_scaled = eps ./ maximum(eps) .* maximum(abs.(field_fixed)) .* 0.3
    plot!(p, x, eps_scaled;
        subplot=i,
        label="ε (scaled)",
        fillrange=0,
        fillalpha=0.3,
        fillcolor=:orange,
        linecolor=:orange,
        alpha=0.5
    )
end

savefig(p, "811_field_1d.png")
println("\nSaved: 811_field_1d.png")
