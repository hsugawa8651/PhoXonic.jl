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
