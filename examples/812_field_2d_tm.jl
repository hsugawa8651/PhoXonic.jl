# Example 812: 2D TM Field Visualization
#
# Demonstrates field reconstruction for 2D TM modes (Ez field).
# Output: 812_field_2d_tm.png

using PhoXonic
using Plots

println("Example 812: 2D TM Field Visualization")
println("=" ^ 50)

# Create 2D photonic crystal: square lattice with dielectric rods
lat = square_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [
    (Circle([0.5, 0.5], 0.2), Dielectric(9.0))
])

# Create solver
solver = Solver(TMWave(), geo, (128, 128); cutoff=7)
println("Resolution: 128x128, Cutoff: 7")
println("Number of plane waves: $(solver.basis.num_pw)")

# Solve at Gamma point
freqs, vecs = solve_at_k_with_vectors(solver, [0.0, 0.0], DenseMethod())
println("\nFirst 5 frequencies: ", round.(freqs[1:5], digits=4))

# Create 2x3 panel plot
p = plot(layout=(2, 3), size=(1200, 800),
    tickfontsize=12, guidefontsize=14, titlefontsize=14,
    bottom_margin=5Plots.mm, top_margin=-2Plots.mm)

# Plot epsilon
eps = get_epsilon_field(solver)
heatmap!(p, eps';
    subplot=1,
    c=:grays,
    title="ε(x,y)",
    aspect_ratio=:equal
)

# Plot first 5 bands
for i in 1:5
    field = reconstruct_field(solver, vecs[:, i])
    field_fixed = fix_phase(field)
    data = real.(field_fixed)

    # Symmetric colormap
    maxval = maximum(abs.(data))

    heatmap!(p, data';
        subplot=i+1,
        c=:RdBu,
        clims=(-maxval, maxval),
        title="Band $i (ω = $(round(freqs[i], digits=3)))",
        aspect_ratio=:equal
    )
end

savefig(p, "812_field_2d_tm.png")
println("\nSaved: 812_field_2d_tm.png")
