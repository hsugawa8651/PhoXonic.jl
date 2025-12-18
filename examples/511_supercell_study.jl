# Last-Modified: 2025-12-15T19:00:00+09:00
# Supercell Study for Non-Periodic Boundary Conditions
#
# This script explores how to implement supercells in PhoXonic
# for simulating finite structures and non-periodic boundaries.

using PhoXonic
using LinearAlgebra
using Printf

println("=" ^ 70)
println("Supercell Study")
println("=" ^ 70)
println()

# ============================================================================
# 1. Basic Supercell Concept
# ============================================================================
println("-" ^ 70)
println("1. Basic Supercell Concept")
println("-" ^ 70)
println()

# Original unit cell
a = 1.0
lat_original = square_lattice(a)
println("Original lattice vectors:")
println("  a1 = $(lat_original.vectors[1])")
println("  a2 = $(lat_original.vectors[2])")
println("  Brillouin zone: k ∈ [-π/a, π/a]")
println()

# 3×3 supercell
Nx, Ny = 3, 3
a1_super = Nx * lat_original.vectors[1]
a2_super = Ny * lat_original.vectors[2]
lat_super = Lattice(a1_super, a2_super)

println("$(Nx)×$(Ny) Supercell lattice vectors:")
println("  A1 = $(lat_super.vectors[1])")
println("  A2 = $(lat_super.vectors[2])")
println("  Brillouin zone: K ∈ [-π/$(Nx)a, π/$(Nx)a]")
println()

# Reciprocal space relation
println("Reciprocal lattice vectors:")
println("  Original: b1 = $(lat_original.reciprocal[1])")
println("  Supercell: B1 = $(lat_super.reciprocal[1])")
println("  Ratio: B1/b1 = $(norm(lat_super.reciprocal[1])/norm(lat_original.reciprocal[1]))")
println()

# ============================================================================
# 2. Replicate Inclusions in Supercell
# ============================================================================
println("-" ^ 70)
println("2. Replicate Inclusions in Supercell")
println("-" ^ 70)
println()

# Original unit cell: one rod at center
r = 0.2
air = Dielectric(1.0)
rod = Dielectric(8.9)

# Replicate the rod at each unit cell position within the supercell
function create_supercell_inclusions_2d(original_inclusions, lat_original, Nx, Ny)
    a1, a2 = lat_original.vectors
    new_inclusions = Tuple{PhoXonic.Shape{PhoXonic.Dim2}, PhoXonic.Material}[]

    for (shape, mat) in original_inclusions
        for i in 0:(Nx-1)
            for j in 0:(Ny-1)
                # Offset for this unit cell
                offset = i * a1 + j * a2

                # Create shifted shape
                if shape isa Circle
                    new_center = shape.center + offset
                    push!(new_inclusions, (Circle(new_center, shape.radius), mat))
                elseif shape isa Rectangle
                    new_center = shape.center + offset
                    push!(new_inclusions, (Rectangle(new_center, shape.size), mat))
                end
            end
        end
    end

    return new_inclusions
end

# Original geometry: rod at (0.5, 0.5) in unit cell
original_inclusions = [(Circle([0.5, 0.5], r), rod)]
supercell_inclusions = create_supercell_inclusions_2d(original_inclusions, lat_original, Nx, Ny)

println("Original: 1 inclusion")
println("Supercell: $(length(supercell_inclusions)) inclusions")
println()

# Show inclusion positions
println("Inclusion positions in supercell:")
for (i, (shape, _)) in enumerate(supercell_inclusions)
    if shape isa Circle
        @printf("  %d: center = (%.2f, %.2f)\n", i, shape.center[1], shape.center[2])
    end
end
println()

# ============================================================================
# 3. Verify Discretization
# ============================================================================
println("-" ^ 70)
println("3. Verify Discretization")
println("-" ^ 70)
println()

# Create supercell geometry
geo_super = Geometry(lat_super, air, supercell_inclusions)

# Discretize
resolution = (Nx * 32, Ny * 32)
ε_super = PhoXonic.discretize(geo_super, resolution, :ε)

# Check fill fraction
rod_pixels = sum(ε_super .> 2.0)
total_pixels = prod(resolution)
fill_fraction = rod_pixels / total_pixels
theoretical_fill = (Nx * Ny * π * r^2) / (Nx * a * Ny * a)

println("Resolution: $(resolution)")
println("Rod pixels: $rod_pixels / $total_pixels")
@printf("Fill fraction: %.4f (theoretical: %.4f)\n", fill_fraction, theoretical_fill)
@printf("Ratio: %.4f\n", fill_fraction / theoretical_fill)
println()

# ============================================================================
# 4. Create Defect in Supercell (Point Defect / Cavity)
# ============================================================================
println("-" ^ 70)
println("4. Create Defect (Cavity)")
println("-" ^ 70)
println()

function create_supercell_with_defect_2d(original_inclusions, lat_original, Nx, Ny;
                                          defect_positions=[])
    """
    Create supercell with defects (missing inclusions).
    defect_positions: list of (i, j) tuples indicating which unit cells have defects
    """
    a1, a2 = lat_original.vectors
    new_inclusions = Tuple{PhoXonic.Shape{PhoXonic.Dim2}, PhoXonic.Material}[]

    for (shape, mat) in original_inclusions
        for i in 0:(Nx-1)
            for j in 0:(Ny-1)
                # Skip defect positions
                if (i, j) in defect_positions
                    continue
                end

                offset = i * a1 + j * a2

                if shape isa Circle
                    new_center = shape.center + offset
                    push!(new_inclusions, (Circle(new_center, shape.radius), mat))
                elseif shape isa Rectangle
                    new_center = shape.center + offset
                    push!(new_inclusions, (Rectangle(new_center, shape.size), mat))
                end
            end
        end
    end

    return new_inclusions
end

# Create 5×5 supercell with center defect (cavity)
Nx_def, Ny_def = 5, 5
center_i, center_j = Nx_def ÷ 2, Ny_def ÷ 2
defect_inclusions = create_supercell_with_defect_2d(
    original_inclusions, lat_original, Nx_def, Ny_def;
    defect_positions=[(center_i, center_j)]
)

println("$(Nx_def)×$(Ny_def) Supercell with point defect:")
println("  Total unit cells: $(Nx_def * Ny_def)")
println("  Defect at: ($center_i, $center_j)")
println("  Inclusions: $(length(defect_inclusions)) (missing 1)")
println()

# ============================================================================
# 5. Create Line Defect (Waveguide)
# ============================================================================
println("-" ^ 70)
println("5. Create Line Defect (Waveguide)")
println("-" ^ 70)
println()

# Create 7×5 supercell with line defect along x (middle row)
Nx_wg, Ny_wg = 7, 5
center_row = Ny_wg ÷ 2
waveguide_defects = [(i, center_row) for i in 0:(Nx_wg-1)]

waveguide_inclusions = create_supercell_with_defect_2d(
    original_inclusions, lat_original, Nx_wg, Ny_wg;
    defect_positions=waveguide_defects
)

println("$(Nx_wg)×$(Ny_wg) Supercell with line defect (waveguide):")
println("  Total unit cells: $(Nx_wg * Ny_wg)")
println("  Defect row: j = $center_row")
println("  Inclusions: $(length(waveguide_inclusions)) (missing $(length(waveguide_defects)))")
println()

# ============================================================================
# 6. Add Absorbing Layer (Simple PML-like)
# ============================================================================
println("-" ^ 70)
println("6. Absorbing Boundary Layer Concept")
println("-" ^ 70)
println()

println("""
For absorbing boundaries, we can add:

1. Simple approach: High-loss material at edges
   - Add imaginary part to ε: ε → ε + i*σ/ω
   - Gradually increase σ towards boundary

2. PML (Perfectly Matched Layer):
   - Requires coordinate stretching in frequency domain
   - More complex to implement in PWE framework

3. Practical approach for PWE:
   - Large supercell + gradual impedance matching
   - Tapered material properties at boundaries

Example implementation concept:
""")

# Concept for adding absorbing layers
function add_absorbing_layer_concept()
    # Inner region: normal crystal
    # Outer region: gradually increasing loss

    # For photonic crystals with complex ε:
    # ε(x) = ε_real + i * σ(x) / ω
    # where σ(x) increases towards boundaries

    println("""
    # Pseudo-code for absorbing layer:

    function create_pml_supercell(inner_geo, n_pml_layers, σ_max)
        # 1. Create larger lattice with PML layers
        total_x = inner_geo.size_x + 2 * n_pml_layers * a
        total_y = inner_geo.size_y + 2 * n_pml_layers * a

        # 2. Copy inner inclusions
        inclusions = copy_to_center(inner_geo.inclusions)

        # 3. Add absorbing material at edges
        for layer in 1:n_pml_layers
            # Conductivity profile (polynomial grading)
            σ = σ_max * (layer / n_pml_layers)^2

            # Add absorbing elements at this layer distance from edge
            # (Implementation depends on how complex ε is handled)
        end

        return Geometry(large_lattice, background, inclusions)
    end
    """)
end

add_absorbing_layer_concept()

# ============================================================================
# 7. Band Folding Effect
# ============================================================================
println("-" ^ 70)
println("7. Band Folding in Supercell")
println("-" ^ 70)
println()

println("""
When using a supercell, the Brillouin zone is folded:

Original BZ: k ∈ [-π/a, π/a]
Supercell BZ (N×): K ∈ [-π/(Na), π/(Na)]

This causes band folding:
- N×N supercell → N² times more bands appear
- Many bands are "folded" copies of original bands
- True defect states appear in the gap region

To identify defect states:
1. Compare with bulk crystal bands
2. Look for states in the band gap
3. Check field localization at defect site

Example: For a perfect supercell (no defect), bands at Γ include
contributions from k = 0, 2π/Na, 4π/Na, ..., (N-1)2π/Na of original BZ.
""")

# ============================================================================
# Summary: Proposed Supercell API
# ============================================================================
println("=" ^ 70)
println("Proposed Supercell API")
println("=" ^ 70)
println()

println("""
# 1. Create supercell from unit cell
supercell = create_supercell(unit_cell_geo, (Nx, Ny))

# 2. Create supercell with defects
supercell = create_supercell(unit_cell_geo, (Nx, Ny);
    point_defects = [(2, 2)],           # Remove inclusion at (2,2)
    line_defects = [(0:6, 2)],          # Remove row j=2
    modified = [(3, 3) => new_shape]    # Replace inclusion
)

# 3. Add boundary layers
supercell = add_boundary_layer(supercell;
    type = :vacuum,      # :vacuum, :absorbing, :pml
    thickness = 3,       # number of unit cells
    grading = :quadratic # for absorbing
)

# 4. Solve for defect modes
solver = Solver(TMWave(), supercell, (Nx*32, Ny*32); cutoff=5)
ω, modes = solve(solver, [0.0, 0.0]; bands=1:20)

# 5. Identify defect states (in gap)
gap_min, gap_max = 0.32, 0.44  # from bulk calculation
defect_modes = findall(gap_min .< ω .< gap_max)
""")
