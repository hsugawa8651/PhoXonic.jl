# Analysis and Post-Processing

## Convergence Check

```julia
for res in [32, 64, 128]
    solver = Solver(TMWave(), geo, (res, res); cutoff=7)
    bands = compute_bands(solver, kpath; bands=1:4)
    gap = find_bandgap(bands, 1, 2)
    println("Resolution $res: gap = $(gap.gap)")
end
```

## Group Velocity

```julia
k = [0.5, 0.0]  # Point in BZ
vg = group_velocity(solver, k; bands=1:4)
# vg[i] = [∂ω/∂kx, ∂ω/∂ky] for band i
```

## K-path Options

### Simple k-paths

```julia
kpath = simple_kpath_square(a=1.0, npoints=30)    # Γ-X-M-Γ
kpath = simple_kpath_hexagonal(a=1.0, npoints=30) # Γ-M-K-Γ
```

### [Brillouin.jl](https://thchr.github.io/Brillouin.jl/stable/) integration

```julia
using Brillouin

# Automatic high-symmetry path
kpi = irrfbz_path(2, lat.vectors)
kpath = kpath_from_brillouin(kpi, npoints=100)
```

## Plotting

Requires [Plots.jl](https://docs.juliaplots.org/stable/):

```julia
using Plots
plot_bands(bands)
```

## Saving and Loading Results

PhoXonic.jl uses [JLD2.jl](https://juliaio.github.io/JLD2.jl/stable/) for saving and loading results.

### Bands

```julia
# Save
bands = compute_bands(solver, kpath)
save_bands("bands.jld2", bands)

# Load
bands_loaded = load_bands("bands.jld2")
```

### Eigenmodes

```julia
# Save with optional k-point and frequencies
ω, modes = solve(solver, k; bands=1:5)
save_modes("modes.jld2", modes; k=k, frequencies=ω)

# Load
data = load_modes("modes.jld2")
# data.modes, data.k, data.frequencies
```

### Material Arrays

```julia
# Photonic: saves ε, ε⁻¹, μ
save_epsilon("epsilon.jld2", solver)

# Phononic: saves C11, C12, C44, ρ
save_epsilon("elastic.jld2", solver_phononic)

# Load
eps = load_epsilon("epsilon.jld2")
# eps.ε, eps.resolution, eps.wave_type
```

### File Contents

| Function | Saved Data |
|----------|------------|
| `save_bands` | frequencies, kpoints, distances, labels, metadata |
| `save_modes` | modes, k (optional), frequencies (optional), metadata |
| `save_epsilon` | Material arrays depend on wave type (see above) |

## API Reference

- [Solver API](api-solver.md) - BandStructure, find_bandgap, group_velocity
- [Plotting API](api-plotting.md) - plot_bands (requires Plots.jl)
- [I/O API](api-io.md) - save/load functions (requires JLD2.jl)
