# Field Visualization

This page describes how to visualize electromagnetic/acoustic fields in PhoXonic.jl.

## Quick Start

```julia
using PhoXonic
using Plots

# Create geometry and solver
lat = square_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [
    (Circle([0.5, 0.5], 0.2), Dielectric(9.0))
])
solver = Solver(TMWave(), geo, (64, 64); cutoff=7)

# Solve at k-point
freqs, vecs = solve_at_k_with_vectors(solver, [0.0, 0.0], DenseMethod())

# Visualize field
plot_field(solver, vecs[:, 1])
```

## Core Functions

### reconstruct_field

Convert an eigenvector (plane wave coefficients) to a real-space field distribution.

```julia
# Basic usage
field = reconstruct_field(solver, vecs[:, 1])

# Custom output grid resolution
field_hires = reconstruct_field(solver, vecs[:, 1]; grid=(128, 128))
```

The return type depends on the wave type:
- **Scalar fields** (TEWave, TMWave, SHWave, Photonic1D, Longitudinal1D): Returns an array of dimension D
- **Vector fields** (PSVWave, TransverseEM, FullElastic): Returns a tuple of arrays, one per component

### get_epsilon_field

Get the permittivity distribution on the real-space grid.

```julia
eps = get_epsilon_field(solver)
```

### get_material_field

Get an arbitrary material property distribution.

```julia
# For phononic: get density
rho = get_material_field(solver, :ρ)

# For photonic: get permeability
mu = get_material_field(solver, :μ)

# Available properties depend on wave type:
# - Photonic: :ε, :μ, :ε_inv, :μ_inv
# - Phononic: :ρ, :C11, :C12, :C44
```

### fix_phase

Normalize the global phase of a complex field for visualization.

```julia
field_fixed = fix_phase(field; method=:max)
```

Available methods:
- `:max` (default): Make the maximum amplitude point real and positive
- `:center`: Make the center point real and positive
- `:mean`: Rotate to minimize imaginary part (mean phase = 0)

After calling `fix_phase`, plotting `real(field)` gives meaningful results.

### field_energy

Compute the energy density of a mode.

```julia
energy = field_energy(solver, vecs[:, 1])
```

Returns the energy density on the real-space grid (always real, non-negative):
- For photonic modes: proportional to |H|² (or |E|² for TM)
- For phononic modes: proportional to ρ|u|²

## Plotting Functions

### plot_field

Visualize a field distribution. Requires Plots.jl.

**1D fields:**
```julia
plot_field(solver, vecs[:, 1];
    quantity=:real,      # :real, :imag, :abs, :phase
    fix_phase_flag=true,
    title="Field",
    xlabel="x/a",
    ylabel="Field"
)
```

**2D fields:**
```julia
plot_field(solver, vecs[:, 1];
    component=:auto,     # For vector fields: :auto, :x, :y, :z
    quantity=:real,
    colormap=:RdBu,
    fix_phase_flag=true,
    title="Field",
    xlabel="x/a",
    ylabel="y/a"
)
```

### plot_epsilon

Visualize the permittivity distribution.

**1D:**
```julia
plot_epsilon(solver;
    title="Permittivity",
    xlabel="x/a",
    ylabel="ε"
)
```

**2D:**
```julia
plot_epsilon(solver;
    colormap=:grays,
    title="Permittivity",
    xlabel="x/a",
    ylabel="y/a"
)
```

## Examples

### Example 811: 1D Field Visualization

```julia
using PhoXonic
using Plots

# Create 1D photonic crystal
lat = lattice_1d(1.0)
geo = Geometry(lat, Dielectric(1.0), [
    (Segment(0.0, 0.3), Dielectric(9.0))
])
solver = Solver(Photonic1D(), geo, (128,); cutoff=10)

# Solve at Gamma point
freqs, vecs = solve_at_k_with_vectors(solver, [0.0], DenseMethod())

# Plot first 4 modes
p = plot(layout=(2, 2), size=(800, 600))
x = range(0, 1, length=128)
eps = get_epsilon_field(solver)

for i in 1:4
    field = reconstruct_field(solver, vecs[:, i])
    field_fixed = fix_phase(field)

    plot!(p, x, real.(field_fixed);
        subplot=i,
        title="Band $i (ω = $(round(freqs[i], digits=3)))",
        label="Field", linewidth=2
    )
    plot!(p, x, eps ./ 10;
        subplot=i, label="ε/10", linestyle=:dash, alpha=0.5
    )
end
```

### Example 812: 2D TM Field Visualization

```julia
using PhoXonic
using Plots

# Create 2D photonic crystal
lat = square_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [
    (Circle([0.5, 0.5], 0.2), Dielectric(9.0))
])
solver = Solver(TMWave(), geo, (64, 64); cutoff=7)

# Solve at Gamma point
freqs, vecs = solve_at_k_with_vectors(solver, [0.0, 0.0], DenseMethod())

# Visualize first mode
field = reconstruct_field(solver, vecs[:, 1])
field_fixed = fix_phase(field)
heatmap(real.(field_fixed)'; c=:RdBu, aspect_ratio=:equal)
```

### Example 813: 2D SH Phononic Field Visualization

```julia
using PhoXonic
using Plots

# Create 2D phononic crystal: steel in epoxy
lat = square_lattice(1.0)
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)

geo = Geometry(lat, epoxy, [
    (Circle([0.5, 0.5], 0.3), steel)
])
solver = Solver(SHWave(), geo, (64, 64); cutoff=7)

# Solve at Gamma point
freqs, vecs = solve_at_k_with_vectors(solver, [0.0, 0.0], DenseMethod())

# Visualize displacement field uz
field = reconstruct_field(solver, vecs[:, 1])
field_fixed = fix_phase(field)
heatmap(real.(field_fixed)'; c=:RdBu, aspect_ratio=:equal)
```

## API Reference

See [Advanced API - Field Visualization](api-advanced.md#Field-Visualization).
