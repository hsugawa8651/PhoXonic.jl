# Advanced API

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

## Band Tracking

The `track_bands` option in [`compute_bands`](@ref) uses eigenvector overlap to maintain band continuity across k-points.

See [Band Tracking](analysis.md#Band-Tracking) for usage examples.

## Field Visualization

See [Field Visualization](field_visualization.md) for usage examples.

### Core Functions

```@docs
reconstruct_field
get_epsilon_field
get_material_field
fix_phase
field_energy
```

### Plotting Functions

```@docs
plot_field
plot_epsilon
plot_field!
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

## Topological Invariants

See [Topological Invariants](topology.md) for usage examples.

### Zak Phase (1D)

```@docs
compute_zak_phase
ZakPhaseResult
```

### Wilson Loop (2D)

```@docs
compute_wilson_spectrum
WilsonSpectrumResult
winding_number
```
