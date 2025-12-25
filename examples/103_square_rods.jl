# Last-Modified: 2025-12-15T23:30:00+09:00
# Square Lattice Photonic Crystal Band Structure
# Equivalent to MPB's sq-rods.ctl
# Dielectric rods in air - both TM and TE polarization

using PhoXonic
using Plots
default(guidefontsize=14, tickfontsize=12, titlefontsize=14, left_margin=10Plots.mm, right_margin=10Plots.mm, top_margin=5Plots.mm, bottom_margin=10Plots.mm)
using Printf
using LinearAlgebra: PosDefException

println("=== Square Lattice Photonic Crystal ===")
println("Dielectric rods (ε=11.56) in air, r=0.2a")
println("Equivalent to MPB sq-rods.ctl")
println()

# Square lattice
a = 1.0
lat = square_lattice(a)

# Geometry: air background with dielectric rod at center
air = Dielectric(1.0)
rod = Dielectric(11.56)  # Same as MPB default
geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])

# Create k-path: Γ → X → M → Γ
kpath = simple_kpath_square(; a=a, npoints=50)

# ============================================================================
# TM polarization (E field parallel to rods)
# ============================================================================
println("Computing TM bands...")
solver_tm = Solver(TMWave(), geo, (64, 64); cutoff=7)
bands_tm = compute_bands(solver_tm, kpath; bands=1:8, verbose=false)

# Find TM band gaps
println("\nTM Band Gaps:")
gaps_tm = find_all_gaps(bands_tm; threshold=0.001)
if isempty(gaps_tm)
    println("  No significant band gaps found.")
else
    for g in gaps_tm
        println(
            "  Gap between bands ",
            g.bands,
            ": ",
            round(g.gap_ratio*100; digits=1),
            "% gap-to-midgap",
        )
    end
end

# ============================================================================
# TE polarization (H field parallel to rods)
# ============================================================================
println("\nComputing TE bands...")
solver_te = Solver(TEWave(), geo, (64, 64); cutoff=7)
bands_te = compute_bands(solver_te, kpath; bands=1:8, verbose=false)

# Find TE band gaps
println("\nTE Band Gaps:")
gaps_te = find_all_gaps(bands_te; threshold=0.001)
if isempty(gaps_te)
    println("  No significant band gaps found.")
else
    for g in gaps_te
        println(
            "  Gap between bands ",
            g.bands,
            ": ",
            round(g.gap_ratio*100; digits=1),
            "% gap-to-midgap",
        )
    end
end

# ============================================================================
# Iterative solver comparison (KrylovKit)
# ============================================================================
println("\n=== Iterative Solver Comparison (KrylovKit) ===")

println("\nComputing TM bands with KrylovKitMethod...")
solver_tm_krylov = Solver(TMWave(), geo, (64, 64), KrylovKitMethod())
bands_tm_krylov = compute_bands(solver_tm_krylov, kpath; bands=1:8, verbose=false)

# Compare Dense vs KrylovKit
max_diff_tm_krylov = maximum(abs.(bands_tm.frequencies .- bands_tm_krylov.frequencies))
println("  TM: max |Dense - KrylovKit| = ", @sprintf("%.2e", max_diff_tm_krylov))

println("\nComputing TE bands with KrylovKitMethod...")
solver_te_krylov = Solver(TEWave(), geo, (64, 64), KrylovKitMethod())
bands_te_krylov = compute_bands(solver_te_krylov, kpath; bands=1:8, verbose=false)

max_diff_te_krylov = maximum(abs.(bands_te.frequencies .- bands_te_krylov.frequencies))
println("  TE: max |Dense - KrylovKit| = ", @sprintf("%.2e", max_diff_te_krylov))

# ============================================================================
# Iterative solver comparison (LOBPCG)
# Note: LOBPCG may fail with PosDefException for some problems.
#       In such cases, use KrylovKitMethod instead.
# ============================================================================
println("\n=== Iterative Solver Comparison (LOBPCG) ===")

try
    println("\nComputing TM bands with LOBPCGMethod...")
    solver_tm_lobpcg = Solver(TMWave(), geo, (64, 64), LOBPCGMethod())
    bands_tm_lobpcg = compute_bands(solver_tm_lobpcg, kpath; bands=1:8, verbose=false)

    max_diff_tm_lobpcg = maximum(abs.(bands_tm.frequencies .- bands_tm_lobpcg.frequencies))
    println("  TM: max |Dense - LOBPCG| = ", @sprintf("%.2e", max_diff_tm_lobpcg))

    println("\nComputing TE bands with LOBPCGMethod...")
    solver_te_lobpcg = Solver(TEWave(), geo, (64, 64), LOBPCGMethod())
    bands_te_lobpcg = compute_bands(solver_te_lobpcg, kpath; bands=1:8, verbose=false)

    max_diff_te_lobpcg = maximum(abs.(bands_te.frequencies .- bands_te_lobpcg.frequencies))
    println("  TE: max |Dense - LOBPCG| = ", @sprintf("%.2e", max_diff_te_lobpcg))
catch e
    if e isa PosDefException
        println("  LOBPCG skipped: PosDefException (use KrylovKitMethod instead)")
    else
        rethrow(e)
    end
end

# ============================================================================
# Create plots
# ============================================================================
dists = bands_tm.distances
label_positions = [dists[i] for (i, _) in bands_tm.labels]
label_names = [l for (_, l) in bands_tm.labels]

# Common y-axis range
ymax = max(maximum(bands_tm.frequencies), maximum(bands_te.frequencies)) * 1.05
ylims_common = (0, ymax)

# TM plot
p_tm = plot(;
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πc)",
    title="TM Bands (ε=11.56 rods)",
    legend=false,
    grid=true,
    size=(700, 500),
    ylims=ylims_common,
    left_margin=5Plots.mm,
)
for b in 1:size(bands_tm.frequencies, 2)
    plot!(p_tm, dists, bands_tm.frequencies[:, b]; linewidth=2, color=:blue)
end
vline!(p_tm, label_positions; color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_tm, label_positions, label_names)

# Highlight TM band gap if exists
gap_tm_1_2 = find_bandgap(bands_tm, 1, 2)
if gap_tm_1_2.gap > 0
    hspan!(
        p_tm, [gap_tm_1_2.max_lower, gap_tm_1_2.min_upper]; alpha=0.2, color=:blue, label=""
    )
end

# TE plot
p_te = plot(;
    xlabel="Wave vector",
    ylabel="Frequency (ωa/2πc)",
    title="TE Bands (ε=11.56 rods)",
    legend=false,
    grid=true,
    size=(700, 500),
    ylims=ylims_common,
)
for b in 1:size(bands_te.frequencies, 2)
    plot!(p_te, dists, bands_te.frequencies[:, b]; linewidth=2, color=:red)
end
vline!(p_te, label_positions; color=:gray, linestyle=:dash, alpha=0.5)
xticks!(p_te, label_positions, label_names)

# Combined plot (side by side)
p_combined = plot(
    p_tm,
    p_te;
    layout=(1, 2),
    size=(1200, 500),
    plot_title="Square Lattice Photonic Crystal (r=0.2a)",
)

# Save plots
savefig(p_tm, joinpath(@__DIR__, "103_square_tm_bands.png"))
savefig(p_te, joinpath(@__DIR__, "103_square_te_bands.png"))
savefig(p_combined, joinpath(@__DIR__, "103_square_bands_combined.png"))

println("\nSaved: 103_square_tm_bands.png")
println("Saved: 103_square_te_bands.png")
println("Saved: 103_square_bands_combined.png")

display(p_combined)
