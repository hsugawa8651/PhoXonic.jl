# 212_maldovan2006_phoxonic.jl
# Simultaneous photonic and phononic band gaps in the same structure
#
# Reference: Maldovan & Thomas, Appl. Phys. B 83, 595 (2006)
# DOI: https://doi.org/10.1007/s00340-006-2241-y
# "Simultaneous complete elastic and electromagnetic band gaps in periodic structures"
#
# Structure: Air holes in Si, triangular lattice, r/a = 0.46
# Figure 3: Simultaneous photonic (a) and phononic (b) band gaps

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhoXonic
using Plots

println("=" ^ 70)
println("Maldovan 2006: Simultaneous Photonic-Phononic Band Gaps")
println("=" ^ 70)

# =============================================================================
# Material Parameters: Silicon
# =============================================================================

# Silicon - Photonic properties
ε_si = 13.0  # Relative permittivity (n ≈ 3.6)

# Silicon - Phononic properties
ρ_si = 2330.0   # kg/m³
c_T = 5360.0    # m/s (transverse velocity)
c_L = 8950.0    # m/s (longitudinal velocity)

# Convert to Lamé parameters
μ_si = ρ_si * c_T^2
λ_si = ρ_si * (c_L^2 - 2 * c_T^2)

# Create materials
si_photonic = Dielectric(ε_si)
si_phononic = IsotropicElastic(; ρ=ρ_si, λ=λ_si, μ=μ_si)

air_photonic = Dielectric(1.0)
void_phononic = ElasticVoid()

println("\nMaterials:")
println("  Si (photonic): ε = $ε_si")
println("  Si (phononic): ρ = $ρ_si kg/m³, c_T = $c_T m/s, c_L = $c_L m/s")
println("  Air/Void: ε = 1.0, ElasticVoid (Tanaka limit)")

# =============================================================================
# Geometry: Triangular lattice, r/a = 0.46
# =============================================================================

a = 1.0
r_over_a = 0.46

lat = hexagonal_lattice(a)
geo_photonic = Geometry(
    lat, si_photonic, [(Circle([0.0, 0.0], r_over_a * a), air_photonic)]
)
geo_phononic = Geometry(
    lat, si_phononic, [(Circle([0.0, 0.0], r_over_a * a), void_phononic)]
)

println("\nGeometry:")
println("  Lattice: Triangular (hexagonal)")
println("  r/a = $r_over_a")
println("  Structure: Air/vacuum holes in Si matrix")

# =============================================================================
# Compute band structures
# =============================================================================

resolution = (64, 64)
cutoff = 8
n_bands = 10

println("\nSolver: resolution=$resolution, cutoff=$cutoff")

kpath = simple_kpath_hexagonal(; npoints=50)

# --- Photonic ---
println("\n--- Photonic Bands ---")
println("  TE...")
solver_te = Solver(TEWave(), geo_photonic, resolution; cutoff=cutoff)
@time bands_te = compute_bands(solver_te, kpath; bands=1:n_bands)

println("  TM...")
solver_tm = Solver(TMWave(), geo_photonic, resolution; cutoff=cutoff)
@time bands_tm = compute_bands(solver_tm, kpath; bands=1:n_bands)

# --- Phononic ---
println("\n--- Phononic Bands ---")
println("  SH (out-of-plane)...")
solver_sh = Solver(SHWave(), geo_phononic, resolution; cutoff=cutoff)
@time bands_sh = compute_bands(solver_sh, kpath; bands=1:n_bands)

println("  PSV (in-plane)...")
solver_psv = Solver(PSVWave(), geo_phononic, resolution; cutoff=cutoff)
@time bands_psv = compute_bands(solver_psv, kpath; bands=1:n_bands)

# =============================================================================
# Normalize frequencies
# =============================================================================

# Photonic: Convert from ωa/c to ωa/2πc
freqs_te = bands_te.frequencies / (2π)
freqs_tm = bands_tm.frequencies / (2π)

# Phononic: ωa/2πc_T
norm_phononic = 1.0 / (2π * c_T)
freqs_sh = bands_sh.frequencies * norm_phononic
freqs_psv = bands_psv.frequencies * norm_phononic

# =============================================================================
# Band Gap Analysis
# =============================================================================

println("\n" * "=" ^ 70)
println("Band Gap Analysis")
println("=" ^ 70)

# Photonic gaps (need to normalize by 2π)
println("\n--- Photonic Gaps (ωa/2πc) ---")
te_gaps = find_all_gaps(bands_te)
tm_gaps = find_all_gaps(bands_tm)

println("TE gaps:")
for gap in te_gaps
    lower = gap.max_lower / (2π)
    upper = gap.min_upper / (2π)
    if lower < 1.0
        println(
            "  Band $(gap.bands): $(round(lower, digits=3)) - $(round(upper, digits=3))"
        )
    end
end

println("TM gaps:")
for gap in tm_gaps
    lower = gap.max_lower / (2π)
    upper = gap.min_upper / (2π)
    if lower < 1.0
        println(
            "  Band $(gap.bands): $(round(lower, digits=3)) - $(round(upper, digits=3))"
        )
    end
end

# Find complete photonic gap (TE ∩ TM overlap)
te_gap_1 =
    isempty(te_gaps) ? nothing : (te_gaps[1].max_lower / (2π), te_gaps[1].min_upper / (2π))
tm_gap_2 = if length(tm_gaps) >= 2
    (tm_gaps[2].max_lower / (2π), tm_gaps[2].min_upper / (2π))
else
    nothing
end
if te_gap_1 !== nothing && tm_gap_2 !== nothing
    overlap_lower = max(te_gap_1[1], tm_gap_2[1])
    overlap_upper = min(te_gap_1[2], tm_gap_2[2])
    if overlap_lower < overlap_upper
        println(
            "\nComplete photonic gap (TE∩TM): $(round(overlap_lower, digits=3)) - $(round(overlap_upper, digits=3))",
        )
    end
end

println("Expected (Maldovan 2006 Fig.3a): Complete photonic gap = 0.408 - 0.464")

# Phononic gaps
println("\n--- Phononic Gaps (ωa/2πcT) ---")
sh_gaps = find_all_gaps(bands_sh)
psv_gaps = find_all_gaps(bands_psv)

println("SH (out-of-plane) gaps:")
for gap in sh_gaps
    lower = gap.max_lower * norm_phononic
    upper = gap.min_upper * norm_phononic
    if lower < 2.0
        println(
            "  Band $(gap.bands): $(round(lower, digits=3)) - $(round(upper, digits=3))"
        )
    end
end

println("PSV (in-plane) gaps:")
for gap in psv_gaps
    lower = gap.max_lower * norm_phononic
    upper = gap.min_upper * norm_phononic
    if lower < 2.0
        println(
            "  Band $(gap.bands): $(round(lower, digits=3)) - $(round(upper, digits=3))"
        )
    end
end

println("\nExpected (Maldovan 2006 Fig.3b): Complete phononic gap = 0.830 - 0.963")

# =============================================================================
# Helper function to draw structure
# =============================================================================
function draw_circle!(p, cx, cy, r; color=:blue, fillalpha=1.0, n=100)
    θ = range(0, 2π; length=n)
    xs = cx .+ r .* cos.(θ)
    ys = cy .+ r .* sin.(θ)
    plot!(
        p,
        xs,
        ys;
        seriestype=:shape,
        fillcolor=color,
        fillalpha=fillalpha,
        linecolor=color,
        linewidth=0.5,
        label="",
    )
end

# =============================================================================
# Plot: Structure + 2x2 band diagrams
# =============================================================================

dists = bands_te.distances
labels = bands_te.labels
label_pos = [dists[i] for (i, _) in labels]
label_names = [l for (_, l) in labels]

# Expected gap regions from paper
photonic_gap = (0.408, 0.464)  # Complete photonic gap
phononic_gap = (0.830, 0.963)  # Complete phononic gap

# --- Structure plot ---
a1 = lat.vectors[1]
a2 = lat.vectors[2]

p_struct = plot(;
    aspect_ratio=:equal,
    xlabel="x/a",
    ylabel="y/a",
    title="Structure",
    xlim=(-0.6, 1.6),
    ylim=(-0.6, 1.2),
    legend=false,
    grid=false,
    background_color_inside=:steelblue,  # Si background
)

# Draw air holes at visible unit cells
for di in -1:2
    for dj in -1:2
        offset = di * a1 + dj * a2
        cx = 0.0 + offset[1]
        cy = 0.0 + offset[2]
        draw_circle!(p_struct, cx, cy, r_over_a; color=:white, fillalpha=1.0)
    end
end

# Draw unit cell boundary
corners = [[0.0, 0.0], collect(a1), collect(a1 + a2), collect(a2), [0.0, 0.0]]
plot!(
    p_struct,
    [c[1] for c in corners],
    [c[2] for c in corners];
    color=:black,
    linewidth=2,
    linestyle=:dash,
    label="",
)

# Add material labels
annotate!(p_struct, 0.8, 0.9, text("Si", :white, 12))
annotate!(p_struct, 0.0, 0.0, text("Air", :gray, 10))

# --- TE subplot ---
p1 = plot(; ylabel="ωa/2πc", title="TE (photonic)", legend=false, grid=true, ylims=(0, 0.8))
for b in 1:size(freqs_te, 2)
    plot!(p1, dists, freqs_te[:, b]; lw=2, color=:blue, label="")
end
hspan!(p1, [photonic_gap...]; alpha=0.3, color=:green, label="")
vline!(p1, label_pos; color=:gray, ls=:dash, alpha=0.5, label="")
xticks!(p1, label_pos, label_names)

# --- TM subplot ---
p2 = plot(; title="TM (photonic)", legend=false, grid=true, ylims=(0, 0.8))
for b in 1:size(freqs_tm, 2)
    plot!(p2, dists, freqs_tm[:, b]; lw=2, color=:red, label="")
end
hspan!(p2, [photonic_gap...]; alpha=0.3, color=:green, label="")
vline!(p2, label_pos; color=:gray, ls=:dash, alpha=0.5, label="")
xticks!(p2, label_pos, label_names)

# --- SH subplot ---
p3 = plot(;
    ylabel="ωa/2πcT", title="SH (phononic)", legend=false, grid=true, ylims=(0, 1.5)
)
for b in 1:size(freqs_sh, 2)
    plot!(p3, dists, freqs_sh[:, b]; lw=2, color=:forestgreen, label="")
end
hspan!(p3, [phononic_gap...]; alpha=0.3, color=:orange, label="")
vline!(p3, label_pos; color=:gray, ls=:dash, alpha=0.5, label="")
xticks!(p3, label_pos, label_names)

# --- PSV subplot ---
p4 = plot(; title="PSV (phononic)", legend=false, grid=true, ylims=(0, 1.5))
for b in 1:size(freqs_psv, 2)
    plot!(p4, dists, freqs_psv[:, b]; lw=2, color=:darkorange, label="")
end
hspan!(p4, [phononic_gap...]; alpha=0.3, color=:orange, label="")
vline!(p4, label_pos; color=:gray, ls=:dash, alpha=0.5, label="")
xticks!(p4, label_pos, label_names)

# --- Combined layout: structure on left, 2x2 bands on right ---
# Use @layout macro for custom layout
lay = @layout [a{0.35w} grid(2, 2)]
p = plot(
    p_struct,
    p1,
    p2,
    p3,
    p4;
    layout=lay,
    size=(1400, 700),
    plot_title="Maldovan 2006: PhoXonic Crystal - Si with Air Holes (r/a=0.46)",
    left_margin=5Plots.mm,
    bottom_margin=5Plots.mm,
)

savefig(p, joinpath(@__DIR__, "212_maldovan2006_phoxonic.png"))
println("\nSaved: 212_maldovan2006_phoxonic.png")

println("\n" * "=" ^ 70)
println("This structure is 'deaf and blind' - it has simultaneous")
println("complete photonic AND phononic band gaps!")
println("=" ^ 70)
println("\nReference: Maldovan & Thomas, Appl. Phys. B 83, 595 (2006)")
