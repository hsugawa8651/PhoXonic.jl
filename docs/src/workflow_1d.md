# 1D Calculations

1D photonic and phononic crystal calculations for multilayer structures (Bragg reflectors, superlattices).

## 1D Wave Types

![1D Wave Types](assets/wave_polarizations_1d.svg)

| Wave Type | Field | Description |
|-----------|-------|-------------|
| [`Photonic1D`](@ref) | E | Electromagnetic wave in dielectric multilayer |
| [`Longitudinal1D`](@ref) | u | Longitudinal elastic wave (P-wave) in elastic multilayer |

## PWE vs TMM

PhoXonic provides two methods for 1D band structure:

| Method | Description | Best For |
|--------|-------------|----------|
| PWE (Plane Wave Expansion) | Fourier-based eigenvalue solver | Band structure, eigenmodes |
| [TMM (Transfer Matrix)](tmm.md) | Exact matrix propagation | Transmission spectra, finite structures |

Both methods give identical band structures for periodic systems. TMM additionally provides transmission/reflection spectra.

## 1D Photonic Crystal

### Bragg Reflector

See also: [`examples/301_bragg_reflector.jl`](https://github.com/hsugawa8651/PhoXonic.jl/blob/main/examples/301_bragg_reflector.jl)

A quarter-wave stack of two dielectric materials.

```julia
using PhoXonic

lat = lattice_1d(1.0)
mat1 = Dielectric(1.0)
mat2 = Dielectric(12.0)

# Quarter-wave stack: d1/d2 = n2/n1
d1 = 0.5 * sqrt(12) / (1 + sqrt(12))
geo = Geometry(lat, mat1, [(Segment(0.0, d1), mat2)])

solver = Solver(Photonic1D(), geo, 128; cutoff=20)

# Compute bands
kpath = simple_kpath_1d(; a=1.0, npoints=50)
bands = compute_bands(solver, kpath; bands=1:5)
```

### Joannopoulos Ch4 Benchmark

See also: [`examples/311_joannopoulos_ch4_fig2.jl`](https://github.com/hsugawa8651/PhoXonic.jl/blob/main/examples/311_joannopoulos_ch4_fig2.jl)

Reproduces Figure 2 from "Photonic Crystals: Molding the Flow of Light" Chapter 4.

## 1D Phononic Crystal

### Elastic Superlattice

See also: [`examples/302_elastic_superlattice.jl`](https://github.com/hsugawa8651/PhoXonic.jl/blob/main/examples/302_elastic_superlattice.jl)

Alternating layers of two elastic materials (e.g., steel/epoxy).

```julia
using PhoXonic

lat = lattice_1d(0.01)  # 1 cm period
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
geo = Geometry(lat, epoxy, [(Segment(0.0, 0.5), steel)])

solver = Solver(Longitudinal1D(), geo, 128; cutoff=20)

# Compute bands
kpath = simple_kpath_1d(; a=0.01, npoints=50)
bands = compute_bands(solver, kpath; bands=1:5)
```

## Iterative Solvers

1D problems support all solver methods:

```julia
# Dense (default)
solver = Solver(Photonic1D(), geo, 128; cutoff=20)

# KrylovKit
solver = Solver(Photonic1D(), geo, 128, KrylovKitMethod(); cutoff=20)

# LOBPCG
solver = Solver(Photonic1D(), geo, 128, LOBPCGMethod(); cutoff=20)
```

For phononic problems with large eigenvalues, `KrylovKitMethod` applies automatic scaling.

## See Also

- [Transfer Matrix Method (TMM)](tmm.md) - Exact solution for transmission/reflection
- [Workflow (2D)](workflow.md) - 2D photonic and phononic crystals
- [Solver Methods](solver.md) - Dense, KrylovKit, LOBPCG comparison
