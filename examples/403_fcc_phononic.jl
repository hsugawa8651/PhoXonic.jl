# 3D FCC Phononic Crystal - Tungsten Carbide spheres in Epoxy
# Reference: Lu et al., Scientific Reports 7, 43407 (2017)
#            "3-D phononic crystals with ultra-wide band gaps"
#
# Material: WC spheres (ρ=13800, E=387.56 GPa, ν=0.346) in Epoxy matrix
# Expected result: ~93.3% normalized band gap between bands 6-7 at f=0.588

using PhoXonic
using LinearAlgebra
using Plots
default(guidefontsize=14, tickfontsize=12, titlefontsize=14, left_margin=10Plots.mm, right_margin=10Plots.mm, top_margin=5Plots.mm, bottom_margin=10Plots.mm)

println("=" ^ 60)
println("3D FCC Phononic Crystal: WC spheres in Epoxy")
println("Reference: Lu et al., Sci. Rep. 7, 43407 (2017)")
println("=" ^ 60)
println()

# FCC lattice
a = 1.0  # conventional lattice constant (normalized)
lat = fcc_lattice(a)

println("Structure: FCC lattice")
println("  Conventional lattice constant: a = $a")
println("  Primitive vectors:")
for (i, v) in enumerate(lat.vectors)
    println("    a$i = ", round.(v; digits=4))
end

# Material parameters from Lu et al. Scientific Reports 7, 43407 (2017)
# Tungsten Carbide (WC): ρ = 13800 kg/m³, E = 387.5559 GPa, ν = 0.3459
# Epoxy:                 ρ = 1180 kg/m³,  E = 4.3438 GPa,   ν = 0.3679
#
# Normalize by epoxy properties
ρ_epoxy_real = 1180.0      # kg/m³
E_epoxy_real = 4.3438e9    # Pa

# Normalized values (epoxy as reference)
ρ_WC = 13800.0 / ρ_epoxy_real    # ~11.69
ρ_epoxy = 1.0

# Create materials using from_E_ν (exact values from paper)
WC = from_E_ν(ρ_WC, 387.5559e9 / E_epoxy_real, 0.3459)
epoxy = from_E_ν(ρ_epoxy, 1.0, 0.3679)

println("\nMaterial parameters (normalized by epoxy):")
println("  WC sphere:   ρ = $(round(WC.ρ, digits=2)), C11 = $(round(WC.C11, digits=2)), C44 = $(round(WC.C44, digits=2))")
println("  Epoxy matrix: ρ = $(round(epoxy.ρ, digits=2)), C11 = $(round(epoxy.C11, digits=2)), C44 = $(round(epoxy.C44, digits=2))")

# Contrast ratio
println("\nContrast ratios:")
println("  ρ_WC / ρ_epoxy = $(round(WC.ρ / epoxy.ρ, digits=1))")
println("  C11_WC / C11_epoxy = $(round(WC.C11 / epoxy.C11, digits=1))")
println("  C44_WC / C44_epoxy = $(round(WC.C44 / epoxy.C44, digits=1))")

# Wave velocities
v_L_WC = sqrt(WC.C11 / WC.ρ)
v_T_WC = sqrt(WC.C44 / WC.ρ)
v_L_epoxy = sqrt(epoxy.C11 / epoxy.ρ)
v_T_epoxy = sqrt(epoxy.C44 / epoxy.ρ)
println("\nWave velocities (normalized):")
println("  WC:    v_L = $(round(v_L_WC, digits=2)), v_T = $(round(v_T_WC, digits=2))")
println("  Epoxy: v_L = $(round(v_L_epoxy, digits=2)), v_T = $(round(v_T_epoxy, digits=2))")

# Optimal filling fraction f = 0.588 for maximum band gap (Lu et al. Table 1)
# For FCC with sphere at origin: f = (4π/3)r³ / V_primitive
# V_primitive = a³/4 for FCC
# f = (4π/3)r³ / (a³/4) = (16π/3)(r/a)³
# For f = 0.588: r/a = (3f/16π)^(1/3) ≈ 0.327
filling_target = 0.588
radius = a * (3 * filling_target / (16π))^(1 / 3)  # ≈ 0.327a

# Volume fraction calculation
V_primitive = a^3 / 4  # FCC primitive cell volume
V_sphere = (4π / 3) * radius^3
filling_fraction = V_sphere / V_primitive

println("\nGeometry:")
println("  Sphere radius: r/a = $(round(radius/a, digits=3))")
println("  Filling fraction: f = $(round(filling_fraction, digits=3))")

# WC spheres in epoxy matrix
geo = Geometry(lat, epoxy, [(Sphere([0.0, 0.0, 0.0], radius), WC)])

# Solver setup (higher resolution for 3D convergence)
resolution = (16, 16, 16)
cutoff = 7
shift = 0.1

println("\nSolver configuration:")
println("  Wave type: FullElastic (3N DOF)")
println("  Resolution: $resolution")
println("  Cutoff: $cutoff")

solver = Solver(FullElastic(), geo, resolution, KrylovKitMethod(; shift=shift); cutoff=cutoff)
println("  Plane waves: ", solver.basis.num_pw)
println("  Matrix size: ", 3 * solver.basis.num_pw, " × ", 3 * solver.basis.num_pw)

# K-path for FCC: Γ → X → W → L → Γ → K
kpath = simple_kpath_fcc(; a=a, npoints=12)
println("\nK-path: Γ → X → W → L → Γ → K")
println("K-points: ", length(kpath))

# Compute bands (need at least 8 bands to see gap between 6-7)
nbands = 10
println("\nComputing $nbands bands...")
println("(This may take several minutes for 3D elastic calculations)")
# Use track_bands=true for smooth band tracking across k-points
bands_result = compute_bands(solver, kpath; bands=1:nbands, verbose=true, track_bands=true)

# Find band gaps
println("\n=== Band Gaps ===")
gaps = find_all_gaps(bands_result; threshold=0.01)
if isempty(gaps)
    println("No complete band gap found.")
    println("Note: May need higher resolution or different parameters")
else
    for g in gaps
        gap_pct = round(g.gap_ratio * 100; digits=1)
        println("Gap between bands $(g.bands): $(gap_pct)% gap-to-midgap")
    end
end

println("\n=== Expected (Lu et al. 2017, Table 1) ===")
println("~93.3% gap-to-midgap for FCC WC/Epoxy spherical inclusions at f=0.588")
println("(Gap between bands 6-7)")

# Normalize frequencies
freqs_normalized = bands_result.frequencies
dists = bands_result.distances

# Plot
p = plot(;
    size=(700, 500),
    xlabel="Wave vector",
    ylabel="Frequency ω (normalized)",
    title="3D FCC Phononic: WC/Epoxy (f=0.588)",
    legend=false,
    grid=true,
    framestyle=:box,
)

# Find and highlight band gap if exists
if !isempty(gaps)
    g = gaps[1]
    gap_bottom = maximum(freqs_normalized[:, g.bands[1]])
    gap_top = minimum(freqs_normalized[:, g.bands[2]])
    plot!(
        p,
        [dists[1], dists[end], dists[end], dists[1], dists[1]],
        [gap_bottom, gap_bottom, gap_top, gap_top, gap_bottom];
        fillrange=gap_bottom,
        fillalpha=0.3,
        fillcolor=:lightgreen,
        linecolor=:green,
        linewidth=1,
        label="",
    )
    annotate!(
        p,
        (dists[1] + dists[end]) / 2,
        (gap_bottom + gap_top) / 2,
        text("Complete Band Gap", 8),
    )
end

# Plot bands as lines
for b in 1:size(freqs_normalized, 2)
    plot!(p, dists, freqs_normalized[:, b]; linewidth=1.5, color=:blue, label="")
end

# Add vertical lines at high-symmetry points
for (idx, label) in bands_result.labels
    vline!(p, [dists[idx]]; color=:gray, linestyle=:solid, alpha=0.5, linewidth=0.5)
end

# Add x-axis labels
label_positions = [dists[idx] for (idx, _) in bands_result.labels]
label_names = [label for (_, label) in bands_result.labels]
plot!(p; xticks=(label_positions, label_names))

output_file = joinpath(@__DIR__, "403_fcc_phononic_bands.png")
savefig(p, output_file)
println("\nSaved: $output_file")
