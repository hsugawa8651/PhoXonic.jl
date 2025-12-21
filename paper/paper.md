---
title: 'PhoXonic.jl: A Julia package for photonic and phononic crystal band structure calculations'
tags:
  - Julia
  - photonic crystals
  - phononic crystals
  - phoxonic crystals
  - plane wave expansion
  - transfer matrix method
  - band structure
authors:
  - name: Hiroharu Sugawara
    orcid: 0000-0002-0071-2396
    corresponding: true
    email: hsugawa@tmu.ac.jp
    affiliation: 1
affiliations:
  - name: Graduate School of Systems Design, Tokyo Metropolitan University, Japan
    index: 1
date: 19 December 2025
bibliography: paper.bib
---

# Summary

Photonic and phononic crystals are periodic structures that exhibit band gaps—frequency ranges where electromagnetic or elastic waves cannot propagate. These materials enable precise control over light and sound at the wavelength scale, with applications ranging from optical fibers and lasers to acoustic filters and vibration isolation. When a single structure simultaneously exhibits both photonic and phononic band gaps, it is called a phoxonic crystal, enabling coupled optomechanical interactions.

`PhoXonic.jl` is a Julia package for computing band structures of photonic, phononic, and phoxonic crystals. It provides a unified interface for simulating electromagnetic and elastic wave propagation in periodic media across one, two, and three dimensions. The package implements both the plane wave expansion (PWE) method for arbitrary periodic geometries in 1D, 2D, and 3D, and the transfer matrix method (TMM) for exact solutions of 1D multilayer structures.

# Statement of Need

Several tools exist for photonic crystal simulations, including MIT Photonic Bands (MPB) [@johnson2001block] and Peacock.jl [@palmer2020peacock]. However, these packages focus exclusively on electromagnetic waves. For phononic crystals, researchers typically rely on commercial finite element method (FEM) software such as COMSOL Multiphysics, or custom implementations. No open-source package provides a unified framework for both photonic and phononic crystals within a single, consistent API using the plane wave expansion method.

`PhoXonic.jl` addresses this gap by offering:

- **Unified API**: The same workflow for photonic and phononic crystals—define materials, geometry, and wave type, then compute band structures.
- **Full dimension support**: 1D, 2D, and 3D simulations, unlike Peacock.jl which is limited to 2D.
- **Multiple solver backends**: Dense eigensolvers for small systems, Krylov methods for large-scale problems, and LOBPCG with warm-start acceleration for efficient band structure sweeps.
- **Green's function methods**: Density of states (DOS) and local density of states (LDOS) calculations for defect mode analysis.
- **Supercell support**: Point and line defect simulations via supercell construction.
- **Void region modeling**: The `ElasticVoid` type enables phononic crystal simulations with vacuum or air inclusions using the Tanaka limit approach [@tanaka2000band].
- **Transfer matrix method**: Exact solutions for 1D multilayer structures, including transmission/reflection spectra, oblique incidence with TE/TM polarization, and support for lossy materials.

The package is implemented in pure Julia, requiring no external compiled dependencies, which simplifies installation and enables interactive development in Jupyter notebooks and the Julia REPL.

# Comparison with Existing Software

| Feature | MPB | Peacock.jl | PhoXonic.jl |
|---------|:---:|:----------:|:-----------:|
| Photonic crystals | Yes | Yes | Yes |
| Phononic crystals | No | No | Yes |
| 1D structures | Yes | No | Yes |
| 2D structures | Yes | Yes | Yes |
| 3D structures | Yes | No | Yes |
| Pure Julia | No | Yes | Yes |
| LDOS/DOS | Limited | No | Yes |
| Supercell defects | Yes | No | Yes |
| Void inclusions | N/A | N/A | Yes |
| Transfer matrix (1D) | No | No | Yes |

# Example Usage

## 2D Photonic Crystal

```julia
using PhoXonic

# Triangular lattice with dielectric rods in air
lat = hexagonal_lattice(1.0)
air = Dielectric(1.0)
rod = Dielectric(12.0)
geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])

# TM mode solver
solver = Solver(TMWave(), geo, (64, 64); cutoff=7)
bands = compute_bands(solver, simple_kpath_hexagonal(); bands=1:8)
gaps = find_all_gaps(bands)
```

## 2D Phononic Crystal

```julia
# Steel cylinders in epoxy matrix
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
steel = IsotropicElastic(ρ=7800.0, λ=1.15e11, μ=8.28e10)

geo = Geometry(square_lattice(1.0), epoxy,
    [(Circle([0.0, 0.0], 0.4), steel)])

# SH (shear horizontal) mode
solver = Solver(SHWave(), geo, (64, 64); cutoff=7)
bands = compute_bands(solver, simple_kpath_square(); bands=1:8)
```

## 3D FCC Photonic Crystal

```julia
lat = fcc_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0),
    [(Sphere([0.0, 0.0, 0.0], 0.25), Dielectric(12.0))])

solver = Solver(FullVectorEM(), geo, (12, 12, 12),
    KrylovKitMethod(shift=0.01); cutoff=3)
bands = compute_bands(solver, simple_kpath_fcc(); bands=1:10)
```

# Validation

`PhoXonic.jl` has been validated against:

- MIT Photonic Bands (MPB) for 3D FCC structures
- Published results from Kushwaha et al. [@kushwaha1993acoustic] for phononic crystals
- Tanaka et al. [@tanaka2000band] for phononic crystals with void inclusions
- Maldovan & Thomas [@maldovan2006simultaneous] for phoxonic crystals (\autoref{fig:phoxonic})
- Textbook examples from Joannopoulos et al. [@joannopoulos2008photonic]

![Phoxonic crystal band structure for Si with air holes (r/a=0.46), reproducing Maldovan & Thomas [@maldovan2006simultaneous]. Left: unit cell structure. Right: photonic (TE, TM) and phononic (SH, PSV) dispersion relations. Frequencies are normalized by the speed of light c (photonic) and the transverse wave velocity cT in silicon (phononic). Shaded regions indicate the complete band gaps reported in the original paper.\label{fig:phoxonic}](212_maldovan2006_phoxonic.png)

The package includes over 1,400 unit tests and 32 example scripts demonstrating various use cases.

# Acknowledgements

We acknowledge contributions from the Julia community and the developers of KrylovKit.jl, IterativeSolvers.jl, and Brillouin.jl, which provide essential functionality for this package.

# References
