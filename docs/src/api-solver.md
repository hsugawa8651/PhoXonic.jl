# Solver API

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

### [Brillouin.jl](https://thchr.github.io/Brillouin.jl/stable/) Integration

```@docs
kpath_from_brillouin
kpath_square
kpath_hexagonal
kpath_cubic
kpath_fcc
kpath_bcc
```
