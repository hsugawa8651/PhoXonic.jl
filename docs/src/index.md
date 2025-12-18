# PhoXonic.jl

Band structure calculation for photonic and phononic crystals using the plane wave expansion (PWE) method.

## Overview

PhoXonic.jl computes eigenfrequencies and eigenmodes of periodic structures:

- **Photonic crystals**: TE/TM modes (2D), full vector (3D)
- **Phononic crystals**: SH/P-SV modes (2D), full elastic (3D)
- **Dimensions**: 1D, 2D, 3D

The name "PhoXonic" comes from **Phoxonic crystals** - structures that exhibit both photonic and phononic band gaps simultaneously.

## Plane Wave Expansion Method

The PWE method solves the eigenvalue problem in reciprocal space by expanding fields in plane waves:

```math
\psi(\mathbf{r}) = \sum_{\mathbf{G}} \psi_{\mathbf{G}} e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}
```

### Photonic Crystals

For electromagnetic waves in periodic dielectric structures, Maxwell's equations lead to:

**TE mode** (H_z):
```math
-\nabla \cdot \left( \frac{1}{\varepsilon} \nabla H_z \right) = \frac{\omega^2}{c^2} \mu H_z
```

**TM mode** (E_z):
```math
-\frac{1}{\mu} \nabla^2 E_z = \frac{\omega^2}{c^2} \varepsilon E_z
```

### Phononic Crystals

For elastic waves in periodic structures:

**SH mode** (u_z, out-of-plane):
```math
-\nabla \cdot (C_{44} \nabla u_z) = \omega^2 \rho u_z
```

**P-SV mode** (u_x, u_y, in-plane):
```math
-\nabla \cdot \boldsymbol{\sigma} = \omega^2 \rho \mathbf{u}
```

## Units

PhoXonic.jl uses **dimensionless units** following the convention of MPB and Peacock.jl:

### Photonic Crystals
- Length: in units of lattice constant ``a``
- Frequency: in units of ``c/a`` (or ``\omega a / 2\pi c``)
- Wave vector: in units of ``2\pi/a``

To convert to physical units:
```
f [Hz] = ω × c / a
λ [m] = a / ω
```

### Phononic Crystals
- Length: in units of lattice constant ``a``
- Frequency: ``\omega`` from eigenvalue (units depend on material constants)

For physical frequencies, material constants should be in consistent SI units (Pa, kg/m³).

## Features

- Multiple lattice types: square, hexagonal, cubic, FCC, BCC
- Subpixel averaging for smooth boundaries
- Band gap detection
- Group velocity computation
- Brillouin.jl integration for k-paths
- Multiple solver methods: Dense, KrylovKit (iterative), LOBPCG
- Matrix-free operators for large-scale 3D calculations
- RSCG for Green's function / DOS / LDOS computation

## Related Projects

- [Peacock.jl](https://github.com/sp94/Peacock.jl) - 2D photonic crystals in Julia (design inspiration). S. J. Palmer and V. Giannini, JOSS 5, 2678 (2020). [DOI:10.21105/joss.02678](https://doi.org/10.21105/joss.02678)
- [MPB](https://mpb.readthedocs.io/) - MIT Photonic-Bands

## References

- J. D. Joannopoulos et al., "Photonic Crystals: Molding the Flow of Light", Princeton University Press (2008). [Book](https://press.princeton.edu/books/hardcover/9780691124568/photonic-crystals)
- M. S. Kushwaha et al., Phys. Rev. Lett. 71, 2022 (1993). [DOI:10.1103/PhysRevLett.71.2022](https://doi.org/10.1103/PhysRevLett.71.2022) - Phononic crystals
- S. G. Johnson and J. D. Joannopoulos, "Block-iterative frequency-domain methods for Maxwell's equations", Opt. Express 8, 173 (2001). [DOI:10.1364/OE.8.000173](https://doi.org/10.1364/OE.8.000173) - MPB
