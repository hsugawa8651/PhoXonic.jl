# API Reference

## Dimensions

```@docs
Dim1
Dim2
Dim3
```

## Lattice

```@docs
Lattice
lattice_1d
square_lattice
hexagonal_lattice
cubic_lattice
fcc_lattice
reciprocal_vectors
```

## Materials

```@docs
Dielectric
IsotropicElastic
from_E_ν
```

## Shapes

```@docs
Circle
Rectangle
Polygon
Sphere
Cylinder
Segment
```

## Geometry

```@docs
Geometry
discretize
DiscretizationMethod
SimpleGrid
SubpixelAverage
```

## Plane Wave Basis

```@docs
PlaneWaveBasis
convolution_matrix
```

## Wave Types

### Abstract Types

```@docs
WaveType
PhotonicWave
PhononicWave
```

### 2D Wave Types

```@docs
TEWave
TMWave
SHWave
PSVWave
```

### 1D Wave Types

```@docs
Photonic1D
Longitudinal1D
```

### 3D Wave Types

```@docs
FullVectorEM
FullElastic
```

## Solver

```@docs
Solver
solve
solve_at_k
matrix_dimension
group_velocity
```

## Solver Method Types

### Abstract Types

```@docs
SolverMethod
IterativeMethod
RSCGMethod
```

### Concrete Methods

```@docs
DenseMethod
BasicRSCG
KrylovKitMethod
LOBPCGMethod
```

## Band Structure

```@docs
BandStructure
compute_bands
find_bandgap
find_all_gaps
frequencies
distances
labels
nbands
nkpoints
```

## K-path

### Simple K-paths

```@docs
SimpleKPath
simple_kpath_square
simple_kpath_hexagonal
simple_kpath_cubic
simple_kpath_fcc
simple_kpath_bcc
```

### Brillouin.jl Integration

```@docs
kpath_from_brillouin
kpath_square
kpath_hexagonal
kpath_cubic
kpath_fcc
kpath_bcc
```

## Matrix-Free Operators

```@docs
FFTContext
MatrixFreeWorkspace
MatrixFreeOperator
apply_lhs!
apply_rhs!
to_linear_map_lhs
to_linear_map_rhs
```

### Effective Hamiltonian

```@docs
EffectiveHamiltonian
MatrixFreeEffectiveHamiltonian
NegatedOperator
```

## Green's Function and DOS/LDOS

!!! note "Optional Dependency"
    `RSKGF` and `MatrixFreeGF` methods require `using ReducedShiftedKrylov`.
    See [Green's Function Methods](greens_function.md) for detailed usage.

### Unified API

```@docs
GFMethod
DirectGF
RSKGF
MatrixFreeGF
```

### RHS Inversion Methods

```@docs
RHSInvMethod
ApproximateRHSInv
CGRHSInv
```

### High-Level Functions

```@docs
compute_greens_function
compute_dos
compute_ldos
```

### Stochastic DOS

```@docs
compute_dos_stochastic
```

## Plotting

Requires `Plots.jl` to be loaded.

```@docs
plot_bands
plot_bands!
band_plot_data
```

## I/O

Requires `JLD2.jl` to be loaded.

```@docs
save_bands
load_bands
save_modes
load_modes
save_epsilon
load_epsilon
```

---

See [Workflow](workflow.md) for detailed usage examples and [Matrix-Free Methods](matrixfree.md) for large-scale calculations.
