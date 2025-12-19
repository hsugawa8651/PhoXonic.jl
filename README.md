# PhoXonic.jl

[![CI](https://github.com/phoxonic-dev/PhoXonic.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/phoxonic-dev/PhoXonic.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/phoxonic-dev/PhoXonic.jl/actions/workflows/Documentation.yml/badge.svg)](https://phoxonic-dev.github.io/PhoXonic.jl/dev/)

Band structure calculation for photonic and phononic crystals using the plane wave expansion (PWE) method.

## Features

- **Photonic crystals**: TE/TM modes (2D), full vector EM (3D)
- **Phononic crystals**: SH/P-SV modes (2D), full elastic (3D)
- **Dimensions**: 1D, 2D, 3D
- **Solvers**: Dense, KrylovKit (iterative), LOBPCG
- **Analysis**: Band gaps, LDOS, Green's function
- **Lattices**: Square, hexagonal/triangular, cubic, FCC

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/phoxonic-dev/PhoXonic.jl")
```

## Quick Start

### 2D Photonic Crystal

```julia
using PhoXonic

# Triangular lattice photonic crystal
lat = hexagonal_lattice(1.0)

# Dielectric rods (ε=12) in air
air = Dielectric(1.0)
rod = Dielectric(12.0)
geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])

# Create TM solver
solver = Solver(TMWave(), geo, (64, 64); cutoff=7)

# Compute band structure along Γ → M → K → Γ
kpath = simple_kpath_hexagonal(npoints=30)
bands = compute_bands(solver, kpath; bands=1:8)

# Find band gaps
gaps = find_all_gaps(bands)
```

### 2D Phononic Crystal

```julia
# Steel cylinders in epoxy matrix
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
steel = IsotropicElastic(ρ=7800.0, λ=1.15e11, μ=8.28e10)

lat = square_lattice(1.0)
geo = Geometry(lat, epoxy, [(Circle([0.0, 0.0], 0.4), steel)])

# SH (out-of-plane) mode
solver = Solver(SHWave(), geo, (64, 64); cutoff=7)
bands = compute_bands(solver, simple_kpath_square(); bands=1:8)
```

### 1D Bragg Reflector

```julia
lat = lattice_1d(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Segment(0.0, 0.25), Dielectric(9.0))])
solver = Solver(Photonic1D(), geo, 128; cutoff=15)
```

### 3D FCC Photonic Crystal

```julia
lat = fcc_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Sphere([0.0, 0.0, 0.0], 0.25), Dielectric(12.0))])
solver = Solver(FullVectorEM(), geo, (12, 12, 12), KrylovKitMethod(shift=0.01); cutoff=3)
```

## Optional Dependencies

| Feature | Package | Description |
|---------|---------|-------------|
| `RSKGF`, `MatrixFreeGF` | [ReducedShiftedKrylov.jl](https://github.com/hsugawa8651/ReducedShiftedKrylov.jl) | Efficient Green's function computation |
| Plotting | Plots.jl | Band structure visualization |

To use RSKGF or MatrixFreeGF methods:

```julia
using PhoXonic
using ReducedShiftedKrylov  # Enables extension automatically

ldos = compute_ldos(solver, position, ω_values, k_points, MatrixFreeGF(); η=0.01)
```

## Examples

Examples are organized by category:

### `examples/` (No extra dependencies)

| Range | Category | Description |
|-------|----------|-------------|
| 1xx | 2D Photonic | Triangular, square, honeycomb lattices |
| 2xx | 2D Phononic | Steel/epoxy, etc. |
| 3xx | 1D Structures | Bragg reflector, elastic superlattice |
| 4xx | 3D Structures | FCC, SC lattices |
| 501 | Defects | Point defects with DirectGF |
| 8xx | Utilities | Structure plotting |
| 9xx | Benchmarks | MPB comparison, Joannopoulos book figures |

### `examples_rsk/` (Requires ReducedShiftedKrylov.jl)

| File | Description |
|------|-------------|
| 131_high_contrast_silicon.jl | RHSInvMethod comparison |
| 502_defect_mode_matrixfree.jl | Matrix-free LDOS for defect modes |

Run an example:

```bash
julia --project=. examples/101_triangular_rods.jl
```

See [`examples/000_examples.md`](examples/000_examples.md) for the full list.

## Wave Types

| Wave Type | Dimension | Description |
|-----------|-----------|-------------|
| `Photonic1D()` | 1D | Scalar EM |
| `TMWave()` | 2D | TM polarization (E_z) |
| `TEWave()` | 2D | TE polarization (H_z) |
| `FullVectorEM()` | 3D | Full vector Maxwell |
| `Longitudinal1D()` | 1D | Longitudinal elastic |
| `SHWave()` | 2D | Shear horizontal (u_z) |
| `PSVWave()` | 2D | P-SV coupled (u_x, u_y) |
| `FullElastic()` | 3D | Full 3D elastic |

## Materials

```julia
# Dielectric
Dielectric(ε)

# Isotropic elastic
IsotropicElastic(ρ=..., λ=..., μ=...)
IsotropicElastic(ρ=..., C11=..., C44=...)
```

## Shapes

| 1D | 2D | 3D |
|----|----|----|
| `Segment(start, stop)` | `Circle(center, r)` | `Sphere(center, r)` |
| | `Rectangle(center, w, h)` | `Cylinder(center, r, h)` |
| | `Polygon(vertices)` | |

## References

- J. D. Joannopoulos et al., "Photonic Crystals: Molding the Flow of Light", Princeton University Press (2008)
- M. S. Kushwaha et al., Phys. Rev. Lett. 71, 2022 (1993) - Phononic crystals

## License

MIT License
