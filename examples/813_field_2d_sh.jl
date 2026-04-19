# Example 813: 2D SH Phononic Field Visualization
#
# Demonstrates field reconstruction for 2D SH (shear horizontal) phononic modes.
# Output: 813_field_2d_sh.png

using PhoXonic
using Plots

println("Example 813: 2D SH Phononic Field Visualization")
println("=" ^ 50)

# Create 2D phononic crystal: square lattice with steel cylinders in epoxy
lat = square_lattice(1.0)

# Materials: epoxy matrix with steel inclusions
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)

geo = Geometry(lat, epoxy, [
    (Circle([0.5, 0.5], 0.3), steel)  # r/a = 0.3
])

# Create solver
solver = Solver(SHWave(), geo, (128, 128); cutoff=7)
println("Resolution: 128x128, Cutoff: 7")
println("Number of plane waves: $(solver.basis.num_pw)")

# Solve at Gamma point
freqs, vecs = solve_at_k_with_vectors(solver, [0.0, 0.0], DenseMethod())
println("\nFirst 5 frequencies (rad/s): ", round.(freqs[1:5], digits=1))

# Create 2x3 panel plot
p = plot(layout=(2, 3), size=(1200, 800),
    tickfontsize=12, guidefontsize=14, titlefontsize=14,
    bottom_margin=5Plots.mm, top_margin=-2Plots.mm, right_margin=5Plots.mm)

# Plot density distribution
rho = get_material_field(solver, :ρ)
heatmap!(p, rho';
    subplot=1,
    c=:viridis,
    title="ρ(x,y) [kg/m³]",
    aspect_ratio=:equal
)

# Plot first 5 bands (displacement field uz)
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
        title="Band $i (ω = $(round(freqs[i], digits=1)))",
        aspect_ratio=:equal
    )
end

savefig(p, "813_field_2d_sh.png")
println("\nSaved: 813_field_2d_sh.png")
