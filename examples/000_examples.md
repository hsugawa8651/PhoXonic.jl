# PhoXonic.jl Examples

Last-Modified: 2025-12-24

## Numbering System

| Range | Category |
|-------|----------|
| 000 | Documentation |
| 1xx | 2D Photonic |
| 2xx | 2D Phononic |
| 3xx | 1D Structures (PWE) |
| 4xx | 3D Structures |
| 5xx | Defects / Supercells |
| 6xx | Transfer Matrix Method (TMM) |
| 8xx | Utilities |
| 9xx | Benchmarks / References |

## File List

### 1xx: 2D Photonic

| File | Description |
|------|-------------|
| 101_triangular_rods.jl | Triangular lattice rods (TM band calculation) |
| 102_triangular_rods_plot.jl | Triangular lattice plotting |
| 103_square_rods.jl | Square lattice rods (TM/TE) |
| 104_honeycomb_rods.jl | Honeycomb lattice |
| 111_triangular_holes.jl | Triangular lattice air holes |
| 121_subpixel_comparison.jl | Subpixel discretization comparison |

### 2xx: 2D Phononic

| File | Description |
|------|-------------|
| 201_phononic_steel_epoxy.jl | Steel/Epoxy phononic crystal (solver comparison) |
| 202_phononic_pb_epoxy_benchmark.jl | Pb/Epoxy benchmark (Aravantinos-Zafiris 2014) |
| 203_vasseur2001_benchmark.jl | Vasseur 2001 benchmark |
| 205_kushwaha1993_benchmark.jl | Kushwaha 1993 benchmark |
| 207_solver_simple.jl | Simple solver usage |
| 208_solver_comparison.jl | Solver comparison |
| 209_warmstart_benchmark.jl | LOBPCG warm start benchmark |
| 210_vacuum_phononic.jl | Vacuum phononic crystal (ElasticVoid) |
| 211_kushwaha1993_fig1.jl | Kushwaha 1993 Fig.1: Ni/Al square lattice |
| 212_kushwaha1993_fig2.jl | Kushwaha 1993 Fig.2: Ni/Al square lattice (different f) |
| 213_tanaka2000_vacuum_al.jl | Tanaka 2000: Al/vacuum (Tanaka limit validation) |
| 214_maldovan2006_phoxonic.jl | Maldovan 2006: Phoxonic crystal (Si/air) |

### 3xx: 1D Structures

| File | Description |
|------|-------------|
| 301_bragg_reflector.jl | Bragg reflector (1D photonic) |
| 302_elastic_superlattice.jl | Elastic superlattice (1D phononic) |
| 311_joannopoulos_ch4_fig2.jl | Joannopoulos textbook Ch.4 Fig.2 |

### 4xx: 3D Structures

3D photonic crystal examples use `TransverseEM` wave type, which is recommended for 3D calculations.
`TransverseEM` uses a 2N×2N transverse basis (eliminating spurious longitudinal modes) and
automatically satisfies ∇·H = 0.

| File | Description |
|------|-------------|
| 401_fcc_spheres.jl | FCC lattice spheres (3D photonic, TransverseEM) |
| 402_sc_spheres.jl | SC lattice spheres (3D photonic, TransverseEM) |
| 403_fcc_phononic.jl | FCC phononic (not included in v1) |
| 411_joannopoulos_ch6_fig3.jl | Diamond air spheres (Joannopoulos Ch.6 Fig.3) |
| 412_joannopoulos_ch6_fig8.jl | Inverse opal (Joannopoulos Ch.6 Fig.8) |
| 413_mpb_diamond.jl | MPB diamond tutorial (dielectric spheres) |

### 5xx: Defects / Supercells

| File | Description |
|------|-------------|
| 501_defect_mode.jl | Defect mode / LDOS (DirectGF) |
| 511_supercell_study.jl | Supercell study |

### 6xx: Transfer Matrix Method (TMM)

| File | Description | Plots |
|------|-------------|-------|
| 601_tmm_bragg_mirror.jl | Bragg mirror transmission spectrum | 601_bragg_*.png |
| 602_tmm_fabry_perot.jl | Fabry-Pérot cavity resonance | 602_fabry_perot_*.png |
| 603_tmm_phononic.jl | Phononic TMM (Steel/Epoxy) | 603_phononic_*.png |
| 604_tmm_vs_pwe.jl | PWE vs TMM comparison | 604_*.png |

### 8xx: Utilities

| File | Description |
|------|-------------|
| 801_plot_structures.jl | Structure plot generation |

### 9xx: Benchmarks / References

| File | Description |
|------|-------------|
| 901_mpb_benchmark.jl | MPB benchmark |
| 911_joannopoulos_ch5_fig2.jl | Joannopoulos Ch.5 Fig.2 (rods) |
| 912_joannopoulos_ch5_fig10.jl | Joannopoulos Ch.5 Fig.10 (holes) |

## Examples Requiring ReducedShiftedKrylov.jl

The following examples require `using ReducedShiftedKrylov` and are located in a **separate folder** `examples_rsk/` (not in `examples/`):

| File (in `examples_rsk/`) | Description |
|---------------------------|-------------|
| 131_high_contrast_silicon.jl | High contrast Si/Air (RHSInvMethod comparison) |
| 502_defect_mode_matrixfree.jl | Defect mode (matrix-free method) |

> **Note**: These files were moved from `examples/` to `examples_rsk/` because they depend on the optional ReducedShiftedKrylov.jl package.

## How to Run

```bash
cd PhoXonic.jl
julia --project=. examples/101_triangular_rods.jl
```

Or from REPL:

```julia
using Pkg
Pkg.activate(".")
include("examples/101_triangular_rods.jl")
```

### Running RSK-dependent examples

```bash
cd PhoXonic.jl
julia --project=. -e '
using ReducedShiftedKrylov
include("examples_rsk/502_defect_mode_matrixfree.jl")
'
```

## Run All Examples

```bash
cd PhoXonic.jl
for f in examples/[1-9]*.jl; do
    echo "=== Running $f ==="
    julia --project=. "$f" || { echo "FAILED: $f"; exit 1; }
done
```

## Long-Running Examples

The following examples may take significant time to complete due to iterative solver comparisons or 3D calculations:

| File | Reason | Estimated Time |
|------|--------|----------------|
| 201_phononic_steel_epoxy.jl | KrylovKit solver comparison | 2-5 min |
| 208_solver_comparison.jl | Multiple solver comparisons | 2-5 min |
| 209_warmstart_benchmark.jl | LOBPCG benchmark | 2-5 min |
| 4xx (3D examples) | 3D band structure | 5-10 min |
| 501_defect_mode.jl | DirectGF calculation | 2-5 min |
