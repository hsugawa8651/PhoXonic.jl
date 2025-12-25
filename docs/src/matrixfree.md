# Matrix-Free Methods

PhoXonic.jl provides matrix-free implementations for large-scale calculations where
storing dense matrices is prohibitive.

## Overview

### Memory and Computational Complexity

| Operation | Dense | Matrix-Free |
|-----------|-------|-------------|
| Memory | O(N²) | O(N) |
| Matrix construction | O(N³) | O(N) |
| Matrix-vector product | O(N²) | O(N log N) |

For a 3D calculation with N = 10,000 plane waves and 3 components:
- **Dense**: ~14 GB for matrices
- **Matrix-free**: ~few MB

### Supported Wave Types

| Wave Type | Dimension | LHS FFTs | RHS FFTs |
|-----------|-----------|----------|----------|
| TEWave | 2D | 4 | 2 |
| TMWave | 2D | 4 | 2 |
| SHWave | 2D | 4 | 2 |
| PSVWave | 2D | 16 | 4 |
| Photonic1D | 1D | 2 | 2 |
| Longitudinal1D | 1D | 2 | 2 |
| FullVectorEM | 3D | 24 | 6 |
| FullElastic | 3D | 42 | 6 |

## Grid Resolution Constraint

!!! warning "Important"
    The solver's grid resolution must be large enough to contain all plane wave indices.

The required resolution in each dimension is:

```
resolution >= 2 * cutoff + 1
```

| Cutoff | Min Resolution | Recommended |
|--------|---------------|-------------|
| 7 (default) | 15 | 16 or 32 |
| 10 | 21 | 32 |
| 15 | 31 | 32 |

### Dense vs Matrix-Free Equivalence

When the resolution constraint is satisfied, matrix-free and dense methods produce
**mathematically identical results** (up to floating-point rounding errors ~10⁻¹⁵).

**Dense method** (convolution matrix):
```
M[i,j] = ε̃(Gᵢ - Gⱼ)
y = M * x                    # O(N²) matrix-vector product
```

**Matrix-free method** (FFT-based):
```
y = FFT{ ε(r) · iFFT{x} }    # O(N log N) via FFT
```

By the convolution theorem of discrete FFT, both compute the same cyclic convolution.

### Why Aliasing Occurs

When the grid is too small, **aliasing** corrupts the results:

```
cutoff=7 → plane wave indices: -7 to +7

With resolution 16×16:
  - Grid can hold indices -8 to +7
  - Problem: multiplying two plane waves extends the range
    Example: G=5 × G=5 → needs G=10, but grid can only hold up to G=7
  - High-frequency components wrap around (cyclic convolution artifact)
  - Result: ~4% error

With resolution 32×32:
  - Grid can hold indices -16 to +15
  - Sufficient margin for all products of plane waves
  - Result: ~10⁻¹⁶ error (machine precision)
```

### Example: Correct vs Incorrect Resolution

```julia
# Incorrect: resolution too small for cutoff=7
solver_bad = Solver(TEWave(), geo, (8, 8); cutoff=7)  # Error ~96%

# Correct: resolution >= 2*7+1 = 15
solver_good = Solver(TEWave(), geo, (16, 16); cutoff=7)  # Error ~10⁻¹⁶
solver_good = Solver(TEWave(), geo, (32, 32); cutoff=7)  # Also correct
```

## Usage with Iterative Eigensolvers

Matrix-free operators are used automatically with iterative solvers:

```julia
using PhoXonic

# Create geometry
lat = square_lattice(1.0)
mat_rod = Dielectric(8.9)
mat_air = Dielectric(1.0)
rod = Circle([0.5, 0.5], 0.2)
geo = Geometry(lat, mat_air, [(rod, mat_rod)])

# Use resolution >= 2*cutoff+1
solver = Solver(TEWave(), geo, (32, 32), KrylovKitMethod(); cutoff=7)

# Solve - uses matrix-free internally
k = [0.5, 0.0]
ω, modes = solve(solver, k; bands=1:5)
```

For 3D calculations with matrix-free methods, `FullVectorEM` with shift-and-invert is used:

```julia
lat = cubic_lattice(1.0)
geo = Geometry(lat, Dielectric(1.0), [(Sphere([0.5,0.5,0.5], 0.3), Dielectric(12.0))])

# Resolution must be >= 2*cutoff+1 in each dimension
solver = Solver(FullVectorEM(), geo, (16, 16, 16), KrylovKitMethod(shift=0.01); cutoff=7)
ω, modes = solve(solver, [0.1, 0.1, 0.1]; bands=1:6)
```

**Note:** For most 3D photonic crystal calculations, `TransverseEM` with `DenseMethod()` is
recommended as it eliminates spurious modes without requiring shift-and-invert. See
[3D Calculations](@ref) for details.

## Unified Green's Function API

PhoXonic provides a unified API for Green's function, DOS, and LDOS calculations.
The same `GFMethod` types work with all three functions.

For detailed documentation on DOS/LDOS functions, see [DOS / LDOS](@ref).

### Available Methods

| Method | Memory | Description |
|--------|--------|-------------|
| [`DirectGF()`](api-advanced.md#PhoXonic.DirectGF) | O(N²) | LU factorization, most accurate |
| [`RSKGF()`](api-advanced.md#PhoXonic.RSKGF) | O(N²) | [ReducedShiftedKrylov.jl](https://github.com/hsugawa8651/ReducedShiftedKrylov.jl) |
| [`MatrixFreeGF()`](api-advanced.md#PhoXonic.MatrixFreeGF) | O(N) | Matrix-free RSCG, best for large systems |

### Basic Usage

```julia
using PhoXonic

solver = Solver(TEWave(), geo, (32, 32); cutoff=7)
k = [0.5, 0.0]
ω_values = [0.3, 0.4, 0.5, 0.6]
k_points = [[0.0, 0.0], [0.5, 0.0], [0.5, 0.5]]

# Default (DirectGF) - no method argument needed
G = compute_greens_function(solver, k, ω_values, source; η=1e-2)
dos = compute_dos(solver, ω_values, k_points; η=1e-2)
ldos = compute_ldos(solver, [0.5, 0.5], ω_values, k_points; η=1e-2)

# With explicit method
G = compute_greens_function(solver, k, ω_values, source, DirectGF(); η=1e-2)
dos = compute_dos(solver, ω_values, k_points, MatrixFreeGF(); η=1e-2)
ldos = compute_ldos(solver, [0.5, 0.5], ω_values, k_points, RSKGF(); η=1e-2)
```

### MatrixFreeGF Options

`MatrixFreeGF` has an `rhs_inv_method` option for RHS⁻¹ application:

| Option | Description | Speed | Accuracy |
|--------|-------------|-------|----------|
| `ApproximateRHSInv()` | Element-wise 1/ε (default) | Fast | Approximate for inhomogeneous media |
| `CGRHSInv()` | Inner CG iteration | Slower | Exact |

```julia
# Fast approximate method (default)
ldos = compute_ldos(solver, pos, ω_values, k_points, MatrixFreeGF(); η=1e-2)

# Exact CG method
ldos = compute_ldos(solver, pos, ω_values, k_points, MatrixFreeGF(rhs_inv_method=CGRHSInv()); η=1e-2)

# CGRHSInv with custom parameters
ldos = compute_ldos(solver, pos, ω_values, k_points,
    MatrixFreeGF(rhs_inv_method=CGRHSInv(atol=1e-12, rtol=1e-12, maxiter=200)); η=1e-2)
```

!!! note "When to use CGRHSInv"
    - `ApproximateRHSInv()` is sufficient for most cases, especially for peak detection
    - Use `CGRHSInv()` when accurate absolute values are needed for high-contrast structures

!!! warning "RSCG Convergence"
    The RSCG solver may not fully converge for all problems.
    This affects `RSKGF` and `MatrixFreeGF` equally.

    Symptoms of incomplete convergence:
    - Relative errors of a few percent compared to direct solve

    Mitigation strategies:
    - Increase `itmax` parameter
    - Use smaller `η` (broadening) values
    - For critical applications, use `DirectGF()`

!!! tip "Practical Usage"
    **`MatrixFreeGF()` is recommended for:**
    - Peak/resonance frequency detection (accurate even without full convergence)
    - Large-scale calculations where Dense methods are memory-prohibitive
    - Qualitative spectral analysis

    **`DirectGF()` is recommended for:**
    - Accurate absolute LDOS/DOS values
    - Quantitative analysis requiring high precision
    - Smaller systems where O(N²) memory is acceptable

### Method Comparison

```julia
# Direct solve - O(N²) memory, O(N³) per frequency
G_direct = compute_greens_function(solver, k, ω_values, source, DirectGF(); η=1e-2)

# Dense RSK - O(N²) memory, efficient for many frequencies
G_rsk = compute_greens_function(solver, k, ω_values, source, RSKGF(); η=1e-2)

# Matrix-free - O(N) memory, O(N log N) per iteration
G_mf = compute_greens_function(solver, k, ω_values, source, MatrixFreeGF(); η=1e-2)

# Matrix-free with exact RHS⁻¹
G_mf_cg = compute_greens_function(solver, k, ω_values, source, MatrixFreeGF(rhs_inv_method=CGRHSInv()); η=1e-2)
```

## Low-Level API

### [FFTContext](api-advanced.md#PhoXonic.FFTContext) and [MatrixFreeWorkspace](api-advanced.md#PhoXonic.MatrixFreeWorkspace)

For optimal performance in loops, separate FFT plan creation from workspace allocation:

```julia
# FFTContext: holds FFT plans (create once, reuse)
ctx = FFTContext(solver)

# MatrixFreeWorkspace: holds work arrays (one per thread)
workspace = MatrixFreeWorkspace(ctx)

# Create operator with explicit context and workspace
op = MatrixFreeOperator(solver, k, ctx, workspace)
```

!!! warning "Thread Safety"
    **FFT plan creation is NOT thread-safe** ([FFTW.jl](https://juliamath.github.io/FFTW.jl/stable/) limitation).
    Create `FFTContext` from a single thread before entering parallel regions.

    **FFT plan execution IS thread-safe**.
    The same `FFTContext` can be shared across threads, but each thread
    must have its own `MatrixFreeWorkspace` to avoid data races.

### [MatrixFreeOperator](api-advanced.md#PhoXonic.MatrixFreeOperator)

```julia
# Simple usage (creates new FFT plans internally)
op = MatrixFreeOperator(solver, k)

# Optimized usage (reuse FFT plans)
ctx = FFTContext(solver)
workspace = MatrixFreeWorkspace(ctx)
op = MatrixFreeOperator(solver, k, ctx, workspace)

# Apply LHS: y = LHS * x
y = zeros(ComplexF64, N)
apply_lhs!(y, op, x)

# Apply RHS: y = RHS * x
apply_rhs!(y, op, x)
```

### [MatrixFreeEffectiveHamiltonian](api-advanced.md#PhoXonic.MatrixFreeEffectiveHamiltonian)

For RSCG integration:

```julia
# Create effective Hamiltonian H = RHS⁻¹ * LHS
op = MatrixFreeOperator(solver, k)
H = MatrixFreeEffectiveHamiltonian(op)

# Negate for RSCG: A = -H
A = NegatedOperator(H)

# Solve (σI + A)x = b for multiple shifts
shifts = [ω^2 + im*η for ω in ω_values]
x_solutions, stats = ReducedShiftedKrylov.rscg(A, b, shifts)
```

### LinearMap Interface

```julia
using LinearMaps  # https://julialinearalgebra.github.io/LinearMaps.jl/stable/

# Convert to LinearMap for use with other iterative solvers
A_lhs = to_linear_map_lhs(op)  # LHS as LinearMap
B_rhs = to_linear_map_rhs(op)  # RHS as LinearMap
```

## Performance Tips

1. **Choose resolution wisely**: Use `resolution = 2^n` for optimal FFT performance
   (e.g., 16, 32, 64 instead of 15, 31, 63).

2. **Reuse FFT plans**: Use `FFTContext` and `MatrixFreeWorkspace` to avoid
   recreating FFT plans in loops:

   ```julia
   ctx = FFTContext(solver)           # Create once (not thread-safe)
   workspace = MatrixFreeWorkspace(ctx)  # One per thread

   for k in k_points
       op = MatrixFreeOperator(solver, k, ctx, workspace)
       # ... use op ...
   end
   ```

3. **Thread safety**: FFT plan execution is thread-safe, but creation is not.
   For parallel k-point loops:

   ```julia
   ctx = FFTContext(solver)  # Create once in main thread
   Threads.@threads for k in k_points
       ws = MatrixFreeWorkspace(ctx)  # Each thread needs its own workspace
       op = MatrixFreeOperator(solver, k, ctx, ws)
       # ... use op ...
   end
   ```

4. **3D calculations**: Matrix-free is essential for 3D. Dense matrices for
   N = 10,000 require ~14 GB, while matrix-free uses only a few MB.

5. **RSCG convergence**: Matrix-free RSCG often converges better than dense RSCG
   due to better numerical conditioning.

## API Reference

- [Advanced API](api-advanced.md) - MatrixFreeOperator, FFTContext, EffectiveHamiltonian
