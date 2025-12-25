# Core API

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

---

See also:
- [Solver API](api-solver.md) - Wave types, Solver, Band structure
- [Advanced API](api-advanced.md) - Matrix-free, Green's function
- [Plotting API](api-plotting.md) (requires Plots.jl)
- [I/O API](api-io.md) (requires JLD2.jl)

