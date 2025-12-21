# Transfer Matrix Method (1D)

The Transfer Matrix Method (TMM) provides exact solutions for wave propagation in 1D multilayer structures. PhoXonic.jl supports TMM for both photonic and phononic crystals.

## Overview

TMM uses an API similar to PWE, making it easy to compare results between the two methods. Both use the same wave types (`Photonic1D`, `Longitudinal1D`) and return compatible `BandStructure` objects, enabling direct validation of 1D calculations.

TMM is particularly useful for:
- Computing transmission and reflection spectra of multilayer structures
- Analyzing Bragg mirrors and Fabry-Pérot cavities
- Band structure calculation for 1D periodic systems
- Handling oblique incidence and polarization effects
- Modeling lossy materials (complex permittivity)

## Quick Start

### Photonic Bragg Mirror

```julia
using PhoXonic

# Define materials
mat_hi = Dielectric(2.5^2)  # n = 2.5
mat_lo = Dielectric(1.45^2) # n = 1.45

# Design wavelength
λ0 = 1.55

# Quarter-wave layers
d_hi = λ0 / (4 * 2.5)
d_lo = λ0 / (4 * 1.45)

# Create multilayer
unit_cell = [Layer(mat_hi, d_hi), Layer(mat_lo, d_lo)]
ml = periodic_multilayer(unit_cell, 10)

# Create solver
solver = TMMSolver(Photonic1D(), ml)

# Compute transmission spectrum
λ_values = range(0.8*λ0, 1.2*λ0, length=201)
R, T = tmm_spectrum(solver, collect(λ_values))
```

### Phononic Superlattice

```julia
using PhoXonic

# Define elastic materials
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)

# Create multilayer
unit_cell = [Layer(steel, 0.005), Layer(epoxy, 0.005)]
ml = periodic_multilayer(unit_cell, 10)

# Create solver for longitudinal waves
solver = TMMSolver(Longitudinal1D(), ml)

# Compute band structure
bands = tmm_bandstructure(solver; k_points=51, bands=1:5)
```

## Layer and Multilayer Types

### Layer

A single layer is defined by its material and thickness:

```julia
Layer(material, thickness)
```

Supported materials:
- `Dielectric(ε)` - lossless dielectric
- `LossyDielectric(ε)` - complex permittivity for absorption
- `IsotropicElastic(ρ, λ, μ)` - elastic material for phononic waves

### Multilayer

A complete structure including surrounding media:

```julia
# Manual construction
Multilayer(layers::Vector{Layer}, left_material, right_material)

# Periodic structure
periodic_multilayer(unit_cell::Vector{Layer}, n_periods)
```

## TMMSolver

The solver is created by specifying the wave type and structure:

```julia
# Photonic (electromagnetic) waves
solver = TMMSolver(Photonic1D(), multilayer)

# Phononic (longitudinal elastic) waves
solver = TMMSolver(Longitudinal1D(), multilayer)
```

## Transmission Spectrum

### Normal Incidence

```julia
# Single wavelength
result = tmm_spectrum(solver, λ)
# Returns: result.R (reflectivity), result.T (transmissivity)

# Multiple wavelengths
R, T = tmm_spectrum(solver, λ_values)
# Returns: vectors of R and T values
```

### Oblique Incidence (Photonic only)

```julia
# TE polarization (s-polarized)
result = tmm_spectrum(solver, λ; angle=π/6, polarization=:TE)

# TM polarization (p-polarized)
result = tmm_spectrum(solver, λ; angle=π/6, polarization=:TM)
```

Features:
- Handles total internal reflection automatically
- Brewster angle effects for TM polarization
- Energy conservation: R + T = 1 (for lossless materials)

## Band Structure

Compute the dispersion relation for periodic structures:

```julia
bands = tmm_bandstructure(solver; k_points=51, bands=1:5)
```

The result is a `BandStructure` object compatible with PWE results:
- `bands.k_points` - wave vectors along the path
- `bands.frequencies` - eigenfrequencies (k_points × n_bands matrix)
- `bands.k_labels` - high-symmetry point labels

### Bloch's Theorem

For a periodic structure with period `a`, the eigenfrequencies satisfy:

```
cos(ka) = Tr(M) / 2
```

where `M` is the transfer matrix for one period. Bandgaps occur where |Tr(M)/2| > 1.

## Lossy Materials

TMM supports complex permittivity for modeling absorption:

```julia
# Complex permittivity: ε = ε' + iε''
mat_lossy = LossyDielectric(complex(4.0, 0.1))

# Use in multilayer structure
ml = Multilayer([Layer(mat_lossy, d)], mat_air, mat_air)
solver = TMMSolver(Photonic1D(), ml)

# For lossy materials: R + T + A = 1
result = tmm_spectrum(solver, λ)
A = 1 - result.R - result.T  # Absorption
```

Note: `LossyDielectric` is only supported for TMM, not PWE.

## Examples

### Bragg Mirror

Source: `examples/601_tmm_bragg_mirror.jl`

Transmission spectrum of a dielectric Bragg mirror with varying number of layer pairs.

```julia
unit_cell = [Layer(mat_hi, d_hi), Layer(mat_lo, d_lo)]
ml = periodic_multilayer(unit_cell, 20)
solver = TMMSolver(Photonic1D(), ml)

λ_values = range(0.8*λ0, 1.2*λ0, length=201)
R, T = tmm_spectrum(solver, collect(λ_values))
```

![Bragg mirror spectrum](../examples/601_bragg_spectrum.png)

---

### Fabry-Pérot Cavity

Source: `examples/602_tmm_fabry_perot.jl`

Resonance behavior of a Fabry-Pérot cavity with Bragg mirrors.

```julia
# Symmetric structure: (HL)^N HH (LH)^N
left_layers = vcat([l for _ in 1:n_pairs for l in left_mirror_cell]...)
defect_layer = Layer(mat_hi, 2*d_hi)
right_layers = vcat([l for _ in 1:n_pairs for l in right_mirror_cell]...)

all_layers = vcat(left_layers, [defect_layer], right_layers)
ml = Multilayer(all_layers, mat_lo, mat_lo)
```

![Fabry-Pérot finesse](../examples/602_fabry_perot_finesse.png)

---

### Phononic Superlattice

Source: `examples/603_tmm_phononic.jl`

Steel/Epoxy multilayer for acoustic wave bandgaps.

```julia
steel = IsotropicElastic(ρ=7800.0, λ=115e9, μ=82e9)
epoxy = IsotropicElastic(ρ=1180.0, λ=4.43e9, μ=1.59e9)

unit_cell = [Layer(steel, d_steel), Layer(epoxy, d_epoxy)]
ml = periodic_multilayer(unit_cell, 10)
solver = TMMSolver(Longitudinal1D(), ml)

bands = tmm_bandstructure(solver; k_points=51, bands=1:4)
```

![Phononic band structure](../examples/603_phononic_bands.png)

---

### PWE vs TMM Comparison

Source: `examples/604_tmm_vs_pwe.jl`

Verifies consistency between PWE and TMM for 1D structures.

```julia
# TMM
solver_tmm = TMMSolver(Photonic1D(), ml)
bands_tmm = tmm_bandstructure(solver_tmm; k_points=21, bands=1:5)

# PWE
solver_pwe = Solver(Photonic1D(), geo, 256; cutoff=20)
```

![PWE vs TMM comparison](../examples/604_photonic_comparison.png)

## Comparison with PWE

| Feature | TMM | PWE |
|---------|-----|-----|
| Dimensionality | 1D only | 1D, 2D, 3D |
| Accuracy (1D) | Exact | Fourier truncation |
| Transmission/Reflection | Direct | Not available |
| Oblique incidence | Yes | No |
| Lossy materials | Yes | No |
| Arbitrary geometry | No (layers only) | Yes |
| Band structure | Via Bloch condition | Direct eigenvalue |

For 1D periodic structures, TMM and PWE should give identical band structures within numerical precision.

## API Reference

### Types

- `Layer{M}` - single layer with material type `M`
- `Multilayer{M}` - complete multilayer structure
- `TMMSolver{W,M}` - solver for wave type `W` and material type `M`

### Functions

- `tmm_spectrum(solver, λ)` - compute R, T at wavelength λ
- `tmm_spectrum(solver, λ_values)` - compute R, T spectrum
- `tmm_bandstructure(solver; ...)` - compute band structure
- `periodic_multilayer(unit_cell, n)` - create periodic structure
- `refractive_index(material)` - get refractive index
- `acoustic_impedance(material)` - get acoustic impedance (phononic)
