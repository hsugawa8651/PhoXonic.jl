# Elastic Void (Tanaka Limit)

This document explains how to handle void (vacuum/air) regions in phononic crystal calculations using the plane wave expansion (PWE) method.

## The Problem: Void in PWE

The plane wave expansion method for elastic waves requires material properties (density ρ, elastic constants C) at every point. For void regions:

- **Vacuum**: ρ = 0, C = 0 (no material)
- **Air**: Very low ρ and C compared to solids

This causes numerical problems because:

1. **Division by zero**: The eigenvalue problem involves terms like C/ρ
2. **Spurious flat bands**: Setting ρ → 0 naively produces unphysical zero-frequency modes
3. **Ill-conditioned matrices**: Large property contrasts cause numerical instability

## The Tanaka Limit Solution

Tanaka et al. (2000) proposed an elegant solution: instead of ρ → 0, use **ρ/C → 0**.

The key insight is that wave velocity v = √(C/ρ). By keeping C finite while making ρ small:

```
v = √(C/ρ) → ∞  as ρ → 0 with C fixed
```

This pushes spurious modes to very high frequencies (ω ~ v → ∞), effectively removing them from the frequency range of interest.

### Implementation in PhoXonic.jl

```julia
void = ElasticVoid(ρ_ratio=1e-7)
```

Internally, this creates a material with:
- μ = 1.0 (shear modulus, normalized)
- λ = 1.0 (Lamé parameter)
- ρ = ρ_ratio × μ = 1e-7

Resulting velocities:
- v_T = √(μ/ρ) = √(1/1e-7) ≈ 3162 (transverse)
- v_L = √((λ+2μ)/ρ) ≈ 5477 (longitudinal)

## Choosing ρ_ratio

The `ρ_ratio` parameter controls the balance between:

| Smaller ρ_ratio | Larger ρ_ratio |
|-----------------|----------------|
| Spurious bands pushed higher | Spurious bands closer to physical bands |
| Better separation of physical/spurious modes | Less separation |
| Potentially worse numerical conditioning | Better numerical stability |
| More accurate results | May contaminate results |

### Recommended Values

| ρ_ratio | v_T | Use case |
|---------|-----|----------|
| 1e-6 | ~1000 | Quick calculations, lower accuracy |
| **1e-7** (default) | ~3162 | **Recommended for most cases** |
| 1e-8 | ~10000 | High accuracy, check for numerical issues |

### Verification

The default `ρ_ratio=1e-7` was validated against:

- **Maldovan 2006** (Si/void, triangular lattice, r/a=0.46)
  - Paper: complete gap at ωa/2πcT = 0.830-0.963
  - PhoXonic: PSV gap at 0.809-0.959 (within 3%)

- **Tanaka 2000** (Al/void, square lattice, f=0.55)
  - Complete elastic gap around ωa/vT ~ 3.5-4.5

## Usage Example

```julia
using PhoXonic

# Silicon matrix with vacuum holes
si = IsotropicElastic(ρ=2330.0, λ=4.67e10, μ=6.69e10)
void = ElasticVoid()  # default ρ_ratio=1e-7

# Triangular lattice
lat = hexagonal_lattice(1.0)
geo = Geometry(lat, si, [(Circle([0.0, 0.0], 0.46), void)])

# Solve
solver = Solver(PSVWave(), geo, (64, 64); cutoff=7)
bands = compute_bands(solver, kpath; bands=1:10)
```

## Identifying Spurious Bands

Spurious bands from void regions typically appear as:
- Flat or nearly flat bands at high frequencies
- Bands with v_T ≈ √(1/ρ_ratio) × (matrix v_T)

If you see suspicious flat bands in your frequency range of interest, try decreasing `ρ_ratio`.

## References

1. Y. Tanaka, Y. Tomoyasu, S. Tamura, "Band structure of acoustic waves in phononic lattices: Two-dimensional composites with large acoustic mismatch," Phys. Rev. B **62**, 7387 (2000)

2. M. Maldovan, E.L. Thomas, "Simultaneous complete elastic and electromagnetic band gaps in periodic structures," Appl. Phys. B **83**, 595 (2006)
