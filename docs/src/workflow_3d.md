# 3D Calculations

3D photonic and phononic crystal calculations have special considerations
due to the vector nature of the fields and the presence of spurious modes.

## 3D Wave Types

![3D Wave Polarizations](assets/wave_polarizations_3d.svg)

| Wave Type | Field Components | Description |
|-----------|-----------------|-------------|
| `TransverseEM` | H_⊥ (2N basis) | **Recommended** for 3D photonic crystals |
| `FullVectorEM` | H_x, H_y, H_z | Full vector H-field (3N basis, includes spurious modes) |
| `FullElastic` | u_x, u_y, u_z | Full elastic (1 P + 2 S modes per k) |

## 3D Photonic Crystals

### TransverseEM vs FullVectorEM

For 3D photonic crystal calculations, `TransverseEM` is **strongly recommended** over `FullVectorEM`:

| Feature | TransverseEM | FullVectorEM |
|---------|--------------|--------------|
| Matrix size | 2N × 2N | 3N × 3N |
| Spurious modes | None (∇·H = 0 enforced) | N longitudinal modes at ω ≈ 0 |
| Shift required | No | Yes (to filter spurious modes) |
| Memory usage | ~44% less | Full |

`TransverseEM` expands the H-field in a transverse polarization basis where each plane wave
has two orthonormal polarization vectors e₁, e₂ perpendicular to (k+G). This automatically
satisfies the transversality constraint ∇·H = 0, eliminating spurious longitudinal modes.

### H-field Formulation

The 3D implementation uses the H-field formulation:

```
∇ × (ε⁻¹ ∇ × H) = (ω/c)² μ H
```

**TransverseEM (Recommended)**

`TransverseEM` expands the H-field in a transverse polarization basis:

```
H_G = h₁(G) e₁(k+G) + h₂(G) e₂(k+G)
```

where e₁ and e₂ are orthonormal vectors perpendicular to (k+G). This construction
automatically satisfies ∇·H = 0, eliminating all spurious longitudinal modes.

```julia
# TransverseEM: no spurious modes, no shift needed
solver = Solver(TransverseEM(), geo, (16, 16, 16), DenseMethod(); cutoff=5)
```

**FullVectorEM (Legacy)**

`FullVectorEM` uses Cartesian components (Hx, Hy, Hz), producing:
- **N longitudinal modes** with ω ≈ 0 (unphysical, violate ∇·H = 0)
- **2N transverse modes** with physical frequencies

When using `FullVectorEM`, use `shift > 0` to filter out spurious modes:

```julia
# DenseMethod: post-hoc filtering
solver = Solver(FullVectorEM(), geo, (12, 12, 12), DenseMethod(shift=0.01); cutoff=7)

# KrylovKitMethod: shift-and-invert
solver = Solver(FullVectorEM(), geo, (16, 16, 16), KrylovKitMethod(shift=0.01); cutoff=7)
```

See [Solver Methods](@ref) for details on shift-and-invert.

### Examples

**FCC Lattice with Spheres**

See also: [`examples/401_fcc_spheres.jl`](https://github.com/hsugawa8651/PhoXonic.jl/blob/main/examples/401_fcc_spheres.jl)

```julia
lat = fcc_lattice(1.0)
air = Dielectric(1.0)
dielectric = Dielectric(12.0)
geo = Geometry(lat, air, [(Sphere([0.0, 0.0, 0.0], 0.25), dielectric)])

# TransverseEM: recommended for 3D photonic crystals
solver = Solver(TransverseEM(), geo, (16, 16, 16), DenseMethod(); cutoff=5)
```

**Simple Cubic Lattice**

See also: [`examples/402_sc_spheres.jl`](https://github.com/hsugawa8651/PhoXonic.jl/blob/main/examples/402_sc_spheres.jl)

```julia
lat = cubic_lattice(1.0)
geo = Geometry(lat, air, [(Sphere([0.0, 0.0, 0.0], 0.3), dielectric)])

solver = Solver(TransverseEM(), geo, (16, 16, 16), DenseMethod(); cutoff=5)
kpath = simple_kpath_cubic(a=1.0, npoints=20)
bands = compute_bands(solver, kpath; bands=1:6)
```

**More Examples**

- [`examples/411_joannopoulos_ch6_fig3.jl`](https://github.com/hsugawa8651/PhoXonic.jl/blob/main/examples/411_joannopoulos_ch6_fig3.jl) - FCC lattice benchmark
- [`examples/412_joannopoulos_ch6_fig8.jl`](https://github.com/hsugawa8651/PhoXonic.jl/blob/main/examples/412_joannopoulos_ch6_fig8.jl) - Diamond lattice
- [`examples/413_mpb_diamond.jl`](https://github.com/hsugawa8651/PhoXonic.jl/blob/main/examples/413_mpb_diamond.jl) - MPB comparison

### MPB Benchmark

PhoXonic 3D calculations have been validated against [MIT Photonic Bands (MPB)](https://mpb.readthedocs.io/).

See also: [`examples/413_mpb_diamond.jl`](https://github.com/hsugawa8651/PhoXonic.jl/blob/main/examples/413_mpb_diamond.jl)

**Simple Cubic (SC) Lattice**

For SC lattice with ε=12 sphere (r=0.3):

| cutoff | Plane waves | Error vs MPB |
|--------|------------|--------------|
| 3 | 123 | 24% |
| 5 | 515 | 7% |
| **7** | **1419** | **< 1%** |

**Recommended settings for SC:**
```julia
resolution = (16, 16, 16)
cutoff = 5
solver = Solver(TransverseEM(), geo, resolution, DenseMethod(); cutoff=cutoff)
```

**FCC Lattice**

For FCC lattice with ε=12 sphere (r=0.25):

**Recommended settings for FCC:**
```julia
resolution = (16, 16, 16)
cutoff = 5
solver = Solver(TransverseEM(), geo, resolution, DenseMethod(); cutoff=cutoff)
```

**Important Notes on MPB Comparison**

1. **Lattice normalization**: MPB internally normalizes primitive lattice vectors to unit length.
   This affects the fill fraction interpretation:
   - PhoXonic FCC primitive vector: |a| = 1/√2 ≈ 0.707
   - MPB FCC primitive vector: |a| = 1.0 (normalized)

2. **Frequency units**:
   - MPB: f = ω·a/(2πc)
   - PhoXonic: ω (angular frequency with c=1)
   - Conversion: ω_PhoXonic = 2π × f_MPB

3. **Fill fraction**: With the same `r=0.25`:
   - PhoXonic: ~26% fill
   - MPB: ~9% fill

   For exact comparison, match fill fractions or mean ε values.

### Cutoff Convergence

For 3D photonic crystals with high dielectric contrast (ε > 10), sufficient plane waves
are required for accurate band structure calculation.

**Recommended cutoff values:**

| Dielectric contrast | Minimum cutoff |
|--------------------|----------------|
| ε < 4 | 3–5 |
| 4 ≤ ε < 10 | 5–7 |
| ε ≥ 10 | 7+ |

**Example convergence test:**

```julia
for cutoff in [3, 5, 7]
    solver = Solver(TransverseEM(), geo, (16,16,16), DenseMethod(); cutoff=cutoff)
    bands = compute_bands(solver, kpath; bands=1:6)
    println("cutoff=$cutoff: ω₁ = $(bands.frequencies[1,1])")
end
```

**Note:** Higher cutoff increases computation time significantly.

## 3D Phononic Crystals

### FullElastic (Experimental)

`FullElastic` is an experimental implementation for 3D phononic crystal calculations.
It computes the full 3D elastic wave equation with three displacement components
(u_x, u_y, u_z), producing 3 modes per k-point (1 longitudinal P-wave + 2 transverse S-waves).

```julia
lat = cubic_lattice(0.01)  # 1 cm period
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)
steel = IsotropicElastic(ρ=7800.0, λ=1.15e11, μ=8.28e10)
geo = Geometry(lat, epoxy, [(Sphere([0.0, 0.0, 0.0], 0.004), steel)])

solver = Solver(FullElastic(), geo, (16, 16, 16), KrylovKitMethod(shift=0.01); cutoff=5)
```

### Why is FullElastic Challenging?

3D phononic calculations are more difficult than 3D photonic calculations for several reasons:

1. **Large contrast in elastic moduli**: Real materials have high elastic moduli contrast
   (e.g., C_steel / C_epoxy ~ 50), which causes numerical challenges regardless of the unit system.

2. **No transverse basis**: Unlike `TransverseEM` for photonic crystals, there is no
   established transverse basis for elastic waves that eliminates spurious modes.
   The shift parameter must be carefully chosen.

3. **Shift parameter sensitivity**: Too small a shift fails to filter spurious modes;
   too large a shift may miss physical low-frequency modes. The optimal value depends
   on material properties and k-point.

4. **High memory and computation time**: 3N × 3N matrices with N > 1000 plane waves
   require significant resources.

5. **Limited benchmarks**: Unlike 3D photonic crystals (validated against MPB),
   3D phononic benchmarks are less standardized in the literature.

## Common Topics

### 3D K-paths

```julia
# FCC lattice
kpath = simple_kpath_fcc(a=1.0, npoints=20)

# Simple cubic
kpath = simple_kpath_cubic(a=1.0, npoints=20)

# BCC lattice
kpath = simple_kpath_bcc(a=1.0, npoints=20)
```

See also: [`simple_kpath_fcc`](@ref), [`simple_kpath_cubic`](@ref), [`simple_kpath_bcc`](@ref)

### Memory Considerations

3D calculations require significantly more memory than 2D:

| Resolution | Plane waves (cutoff=7) | Dense matrix memory |
|------------|----------------------|---------------------|
| 12×12×12 | ~500 | ~6 GB |
| 16×16×16 | ~1400 | ~47 GB |
| 20×20×20 | ~2700 | ~175 GB |

For large 3D systems, use matrix-free methods. See [Matrix-Free Methods](@ref).

## API Reference

- [Solver API](api-solver.md) - TransverseEM, FullVectorEM, FullElastic, Solver methods
- [Advanced API](api-advanced.md) - Matrix-free operators
