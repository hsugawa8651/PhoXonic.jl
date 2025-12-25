# DOS / LDOS

Density of States (DOS) and Local Density of States (LDOS) computation using Green's function methods.

For large-scale calculations, `MatrixFreeGF()` uses FFT-based operators described in [Matrix-Free Methods](@ref).

## Setup

```julia
using PhoXonic

# For DirectGF (no additional dependencies required)
G = compute_greens_function(solver, k, ω_values, source, DirectGF(); η=0.01)

# For RSKGF and MatrixFreeGF, load ReducedShiftedKrylov.jl
using ReducedShiftedKrylov  # Automatically enables extension

G = compute_greens_function(solver, k, ω_values, source, RSKGF(); η=0.01)
G = compute_greens_function(solver, k, ω_values, source, MatrixFreeGF(); η=0.01)
```

!!! note "Optional Dependency"
    `RSKGF` and `MatrixFreeGF` methods require the [ReducedShiftedKrylov.jl](https://github.com/hsugawa8651/ReducedShiftedKrylov.jl) package.
    When you load it with `using ReducedShiftedKrylov`, the extension is automatically enabled.

## Available Methods

| Method | Memory | Best For | Requires RSK |
|--------|--------|----------|:------------:|
| [`DirectGF()`](api-advanced.md#PhoXonic.DirectGF) | O(N²) | Small systems, high accuracy | No |
| [`RSKGF()`](api-advanced.md#PhoXonic.RSKGF) | O(N²) | Many frequencies, 2D | Yes |
| [`MatrixFreeGF()`](api-advanced.md#PhoXonic.MatrixFreeGF) | O(N) | Large systems, 2D/3D | Yes |

### DirectGF

Direct LU factorization method. Most accurate but requires O(N²) memory for the matrix.

```julia
method = DirectGF()
```

### RSKGF

Uses the Reduced Shifted Conjugate Gradient (RSCG) method from ReducedShiftedKrylov.jl.
Efficient for computing Green's functions at many frequency points simultaneously.

```julia
method = RSKGF(; atol=1e-10, rtol=1e-10, itmax=1000, verbose=false)
```

### MatrixFreeGF

Matrix-free implementation using FFT-based operators. O(N) memory, suitable for large-scale 2D and 3D calculations.

```julia
# Default: ApproximateRHSInv (fast, good for low contrast)
method = MatrixFreeGF()

# Explicit ApproximateRHSInv
method = MatrixFreeGF(rhs_inv_method=ApproximateRHSInv())

# CGRHSInv for high-contrast materials (more accurate)
method = MatrixFreeGF(rhs_inv_method=CGRHSInv())

# CGRHSInv with custom parameters
method = MatrixFreeGF(rhs_inv_method=CGRHSInv(atol=1e-12, rtol=1e-12, maxiter=200))
```

## RHSInvMethod Options

For [`MatrixFreeGF`](api-advanced.md#PhoXonic.MatrixFreeGF), you can choose how to handle the RHS (right-hand side) matrix inversion:

| Method | Description | Best For |
|--------|-------------|----------|
| [`ApproximateRHSInv()`](api-advanced.md#PhoXonic.ApproximateRHSInv) | Element-wise 1/ε in real space | Low contrast materials |
| [`CGRHSInv()`](api-advanced.md#PhoXonic.CGRHSInv) | Iterative CG solve | High contrast materials |

```julia
# CGRHSInv parameters
CGRHSInv(; atol=1e-10, rtol=1e-10, maxiter=100)
```

## Functions

### [compute_greens_function](api-advanced.md#PhoXonic.compute_greens_function)

Compute the Green's function G(ω) = (ω² - H)⁻¹ for multiple frequencies.

```julia
G_values = compute_greens_function(solver, k, ω_values, source, method; η=1e-3)
```

**Arguments:**
- `solver`: Solver instance
- `k`: Wave vector (e.g., `[0.0, 0.0]`)
- `ω_values`: Vector of frequencies
- `source`: Source vector in plane wave basis
- `method`: `DirectGF()`, `RSKGF()`, or `MatrixFreeGF()`
- `η`: Broadening parameter (imaginary part of ω²)

**Returns:** Vector of Green's function solutions, one for each frequency.

### [compute_ldos](api-advanced.md#PhoXonic.compute_ldos)

Compute the Local Density of States at a specific position.

```julia
ldos = compute_ldos(solver, position, ω_values, k_points, method; η=1e-3)
```

**Arguments:**
- `solver`: Solver instance
- `position`: Real-space position (e.g., `[0.5, 0.5]`)
- `ω_values`: Vector of frequencies
- `k_points`: Vector of k-points for BZ sampling
- `method`: `DirectGF()`, `RSKGF()`, or `MatrixFreeGF()`
- `η`: Broadening parameter

**Returns:** Vector of LDOS values for each frequency.

### [compute_dos](api-advanced.md#PhoXonic.compute_dos)

Compute the Density of States (stochastic trace method for RSKGF/MatrixFreeGF).

```julia
dos = compute_dos(solver, ω_values, k_points, method; η=1e-3, n_random=10)
```

**Arguments:**
- `solver`: Solver instance
- `ω_values`: Vector of frequencies
- `k_points`: Vector of k-points for BZ sampling
- `method`: `DirectGF()`, `RSKGF()`, or `MatrixFreeGF()`
- `η`: Broadening parameter
- `n_random`: Number of random vectors for stochastic trace (RSKGF/MatrixFreeGF only)

**Returns:** Vector of DOS values for each frequency.

## Examples

See also: [`examples/501_defect_mode.jl`](https://github.com/hsugawa8651/PhoXonic.jl/blob/main/examples/501_defect_mode.jl) - Defect mode finding with LDOS

### Basic LDOS Calculation

```julia
using PhoXonic

# Setup
lat = square_lattice(1.0)
air = Dielectric(1.0)
rod = Dielectric(4.0)
geo = Geometry(lat, air, [(Circle([0.0, 0.0], 0.2), rod)])
solver = Solver(TEWave(), geo, (32, 32); cutoff=5)

# LDOS with DirectGF (no extra dependencies)
position = [0.5, 0.5]
ω_values = collect(0.2:0.01:0.8)
k_points = [[0.0, 0.0]]

ldos = compute_ldos(solver, position, ω_values, k_points, DirectGF(); η=0.01)
```

### Matrix-Free LDOS (requires ReducedShiftedKrylov)

```julia
using PhoXonic
using ReducedShiftedKrylov  # Enable extension

# ... (same setup as above)

# Matrix-free LDOS for large systems
ldos_mf = compute_ldos(solver, position, ω_values, k_points, MatrixFreeGF(); η=0.01)

# With CGRHSInv for high-contrast materials
ldos_cg = compute_ldos(solver, position, ω_values, k_points,
                        MatrixFreeGF(rhs_inv_method=CGRHSInv()); η=0.01)
```

## Dimension Support

| Method | 1D | 2D | 3D |
|--------|:--:|:--:|:--:|
| [`DirectGF()`](api-advanced.md#PhoXonic.DirectGF) | Yes | Yes | No |
| [`RSKGF()`](api-advanced.md#PhoXonic.RSKGF) | No | Yes | No |
| [`MatrixFreeGF()`](api-advanced.md#PhoXonic.MatrixFreeGF) | No | Yes | Yes |

## See Also

- [Advanced API](api-advanced.md) - GFMethod, compute_dos, compute_ldos
- [Matrix-Free Methods](matrixfree.md) for operators and FFT context
