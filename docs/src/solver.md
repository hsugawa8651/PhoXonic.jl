# Solver Methods

PhoXonic.jl provides multiple solver methods for different problem sizes:

| Method | Matrix Storage | Algorithm | Best For |
|--------|---------------|-----------|----------|
| [`DenseMethod()`](api-solver.md#PhoXonic.DenseMethod) | Dense (N×N) | LAPACK `eigen` | Small/medium systems |
| [`DenseMethod(shift=σ)`](api-solver.md#PhoXonic.DenseMethod) | Dense (N×N) | LAPACK `eigen` + filter | 3D (skip spurious modes) |
| [`KrylovKitMethod()`](api-solver.md#PhoXonic.KrylovKitMethod) | Matrix-free | Arnoldi iteration | Large 2D systems |
| [`KrylovKitMethod(shift=σ)`](api-solver.md#PhoXonic.KrylovKitMethod) | Dense + LU | Shift-and-invert | 3D, targeted frequencies |
| [`LOBPCGMethod()`](api-solver.md#PhoXonic.LOBPCGMethod) | Dense (N×N) | Block CG | Symmetric problems |
| [`LOBPCGMethod(shift=σ)`](api-solver.md#PhoXonic.LOBPCGMethod) | Dense + LU | Shift-and-invert LOBPCG | 3D systems |

## Memory and Computational Cost

| Method | Memory | Cost per k-point |
|--------|--------|------------------|
| Dense | O(N²) | O(N³) eigendecomposition |
| Matrix-free | O(N) | O(N log N) × iterations |
| Shift-invert | O(N²) | O(N³) LU + O(N²) × iterations |
| LOBPCG | O(N²) | O(N² × nev) × iterations |

**N** = number of plane waves = `basis.num_pw` × `ncomponents(wave)`

For large-scale calculations with N > 10,000, see [Matrix-Free Methods](@ref).

## Recommendations

| Dimension | Wave Type | N (typical) | Recommended Method |
|-----------|-----------|-------------|-------------------|
| 1D | [`Photonic1D`](api-solver.md#PhoXonic.Photonic1D) | < 100 | [`DenseMethod()`](api-solver.md#PhoXonic.DenseMethod) |
| 2D | [`TEWave`](api-solver.md#PhoXonic.TEWave)/[`TMWave`](api-solver.md#PhoXonic.TMWave)/[`SHWave`](api-solver.md#PhoXonic.SHWave) | 100–1,000 | [`DenseMethod()`](api-solver.md#PhoXonic.DenseMethod) |
| 2D | TE/TM/SH | 1,000–10,000 | [`KrylovKitMethod()`](api-solver.md#PhoXonic.KrylovKitMethod) or [`LOBPCGMethod()`](api-solver.md#PhoXonic.LOBPCGMethod) |
| 2D | [`PSVWave`](api-solver.md#PhoXonic.PSVWave) | 200–2,000 | [`DenseMethod()`](api-solver.md#PhoXonic.DenseMethod) |
| 2D | PSVWave | > 2,000 | [`KrylovKitMethod()`](api-solver.md#PhoXonic.KrylovKitMethod) or [`LOBPCGMethod()`](api-solver.md#PhoXonic.LOBPCGMethod) |
| 3D | [`FullVectorEM`](api-solver.md#PhoXonic.FullVectorEM) | < 500 | [`DenseMethod(shift=0.01)`](api-solver.md#PhoXonic.DenseMethod) |
| 3D | FullVectorEM | 500–5,000 | [`KrylovKitMethod(shift=0.01)`](api-solver.md#PhoXonic.KrylovKitMethod) or [`LOBPCGMethod(shift=0.01)`](api-solver.md#PhoXonic.LOBPCGMethod) |
| 3D | [`FullElastic`](api-solver.md#PhoXonic.FullElastic) | < 500 | [`DenseMethod(shift=0.01)`](api-solver.md#PhoXonic.DenseMethod) |
| 3D | FullElastic | 500–5,000 | [`KrylovKitMethod(shift=0.01)`](api-solver.md#PhoXonic.KrylovKitMethod) or [`LOBPCGMethod(shift=0.01)`](api-solver.md#PhoXonic.LOBPCGMethod) |

## DenseMethod

The default solver using LAPACK's dense eigenvalue decomposition.

### Basic Usage

```julia
solver = Solver(TEWave(), geo, (64, 64), DenseMethod(); cutoff=7)
```

### DenseMethod for 3D

For small 3D systems, `DenseMethod` with `shift > 0` can be used:

```julia
solver = Solver(FullVectorEM(), geo, (12, 12, 12), DenseMethod(shift=0.01); cutoff=5)  # Use cutoff≥7 for high ε
```

**How it works:**

Unlike iterative methods which use shift-and-invert transformation, `DenseMethod(shift=σ)`
computes **all** eigenvalues using LAPACK, then filters out eigenvalues where ω² < σ.
This post-hoc filtering effectively removes the spurious longitudinal modes at ω ≈ 0.

**Important difference from iterative methods:**

| Aspect | DenseMethod | KrylovKitMethod/LOBPCGMethod |
|--------|-------------|------------------------------|
| Computation | All eigenvalues | Only requested number |
| Filtering | Post-hoc (ω² < σ removed) | Shift-and-invert (finds λ closest to σ) |
| Band ordering | Strictly ascending ω | Depends on σ proximity |
| Memory | O(N²) | O(N²) with shift, O(N) without |

This means Dense and iterative methods may return **different bands** for the same
`bands=1:n` request if there are intermediate eigenvalues. Both results are valid
eigenvalues, just in different order.

**When to use DenseMethod for 3D:**
- Small systems (N × 3 < 2000)
- Validation and debugging
- Need all eigenvalues in ascending order

## KrylovKitMethod

Iterative solver based on KrylovKit.jl using Arnoldi iteration.

### Basic Usage

```julia
solver = Solver(TEWave(), geo, (128, 128), KrylovKitMethod(); cutoff=15)
```

### Shift-and-Invert for 3D

For 3D calculations, use shift-and-invert to skip longitudinal modes:

```julia
solver = Solver(FullVectorEM(), geo, (16, 16, 16), KrylovKitMethod(shift=0.01); cutoff=7)
```

The shift parameter σ transforms the generalized eigenvalue problem:

| Original | Transformed |
|----------|-------------|
| A x = λ B x | (A - σB)⁻¹ B x = μ x |
| λ = eigenvalue | μ = 1/(λ - σ) |

Eigenvalues closest to σ become the largest |μ|, making them easiest for iterative solvers to find.

### Why is Shift Needed for 3D?

The 3D H-field formulation does not explicitly enforce ∇·H = 0, producing:
- **Longitudinal modes**: λ ≈ 0 (unphysical, violate transversality)
- **Transverse modes**: λ > 0 (physical electromagnetic waves)

Without shift, iterative solvers converge to the spurious λ ≈ 0 modes first.
With shift σ > 0, eigenvalues λ < σ are effectively filtered out.

### Choosing σ

| Use case | Recommended σ | Explanation |
|----------|---------------|-------------|
| Skip longitudinal modes | 0.001 – 0.1 | Must be > 0, but smaller than the smallest physical λ |
| Target specific frequency | ω²_target | Finds modes near the target frequency first |
| Low-frequency modes near k | \|k\|²/ε_max | Approximate lowest transverse eigenvalue |

**Guidelines for choosing σ:**

1. **Lower bound**: σ must be larger than the longitudinal mode eigenvalues (≈ 0).
   A small positive value like `σ = 0.01` typically works.

2. **Upper bound**: σ should be smaller than the smallest physical eigenvalue you want to find.
   For a homogeneous medium: λ_min ≈ |k|²/ε_max.

   Example: |k| = 0.5, ε = 4 → λ_min ≈ 0.0625. Use σ < 0.06.

3. **Safe default**: For most 3D photonic calculations, `σ = 0.01` works well when
   |k| ≥ 0.1 (in units of 2π/a).

4. **Near k = 0**: At the Γ point (k = 0), the lowest transverse modes also have
   λ → 0. Use a very small shift (σ ~ 0.001) or skip the Γ point.

5. **High-contrast materials**: For large ε, eigenvalues scale as 1/ε.
   Reduce σ proportionally: `σ ≈ 0.01 / ε_max`.

### Troubleshooting

| Symptom | Likely cause | Solution |
|---------|--------------|----------|
| All eigenvalues ≈ 0 | σ too large | Reduce σ |
| Missing low-frequency modes | σ too large | Reduce σ below λ_min |
| Converges to spurious modes | σ too small or = 0 | Increase σ > 0 |
| Slow convergence | σ far from target λ | Adjust σ closer to desired eigenvalues |

### Phononic Eigenvalue Scaling

Phononic eigenvalues ω² are typically O(10¹⁰) for common materials (steel, epoxy, etc.),
which can cause numerical instability in iterative solvers. PhoXonic.jl automatically
scales the eigenvalue problem for `SHWave`, `PSVWave`, `FullElastic`, and `Longitudinal1D`:

```
Scaled problem: (A/s) x = (λ/s) B x
where s = c² × |k|² and c² = C₁₁/ρ (longitudinal wave speed squared)
```

This normalizes eigenvalues to O(1), ensuring stable convergence. The scaling is:
- **Automatic**: No user configuration needed
- **Transparent**: Results are unscaled before return
- **Photonic-safe**: Scale factor = 1.0 for TE/TM/FullVectorEM

**Example:** Steel/epoxy at |k| = 100 rad/m:
- Without scaling: ω² ~ 10¹⁰ → KrylovKit may fail
- With scaling: ω²/s ~ 1 → stable convergence

### KrylovKitMethod Parameters

```julia
method = KrylovKitMethod(
    tol = 1e-8,       # Convergence tolerance
    maxiter = 300,    # Maximum iterations
    krylovdim = 30,   # Krylov subspace dimension
    verbosity = 0,    # 0=silent, 1=warnings, 2=info
    shift = 0.0       # Spectral shift (0 = no shift)
)
```

## LOBPCGMethod

LOBPCG (Locally Optimal Block Preconditioned Conjugate Gradient) is an alternative
iterative solver. A. V. Knyazev, SIAM J. Sci. Comput. 23, 517 (2001). [DOI:10.1137/S1064827500366124](https://doi.org/10.1137/S1064827500366124)

### When to Use LOBPCG

LOBPCG is particularly effective for:
- Symmetric generalized eigenvalue problems A x = λ B x
- Problems where B is positive definite (mass matrix)
- Computing multiple eigenvalues simultaneously (block method)
- Phononic calculations (no eigenvalue scaling required)

### Basic Usage

```julia
# 2D photonic
solver = Solver(TEWave(), geo, (64, 64), LOBPCGMethod())

# 2D phononic
solver = Solver(SHWave(), geo, (64, 64), LOBPCGMethod())
solver = Solver(PSVWave(), geo, (64, 64), LOBPCGMethod())
```

### Shift-and-Invert for 3D

For 3D calculations, LOBPCG requires `shift > 0` to skip spurious longitudinal modes:

```julia
solver = Solver(FullVectorEM(), geo, (16, 16, 16), LOBPCGMethod(shift=0.01); cutoff=7)

solver = Solver(FullElastic(), geo, (16, 16, 16), LOBPCGMethod(shift=0.01); cutoff=5)
```

**Why is shift needed for 3D?**

The 3D H-field formulation does not enforce the transversality constraint ∇·H = 0,
leading to spurious longitudinal modes at λ ≈ 0. Without shift:
- LOBPCG tries to find smallest eigenvalues, converging to these spurious modes
- The RHS matrix B becomes singular in this subspace → `PosDefException`

With `shift=σ`, the problem transforms to `(A - σB)⁻¹ B x = μ x`, which only
returns eigenvalues λ > σ, effectively filtering out the spurious modes.

### LOBPCGMethod Parameters

```julia
method = LOBPCGMethod(
    tol = 1e-4,              # Convergence tolerance
    maxiter = 200,           # Maximum iterations
    shift = 0.0,             # Spectral shift (0 = no shift, required for 3D)
    warm_start = true,       # Use previous eigenvectors as initial guess
    scale = true,            # Scale matrix A for better conditioning
    first_dense = true,      # Solve first k-point with Dense
    preconditioner = :diagonal  # Preconditioner (:none, :diagonal, or custom)
)
```

### Warm Start for Band Structure Calculations

When computing band structures with `compute_bands`, LOBPCG can use **warm start**
to significantly speed up calculations. This uses the eigenvectors from the previous
k-point as the initial guess for the next k-point.

**How it works:**

1. First k-point is solved with `DenseMethod` for accurate eigenvectors (when `first_dense=true`)
2. Matrix A is scaled by `max|A|` for better conditioning (when `scale=true`)
3. Subsequent k-points use previous eigenvectors as initial guess
4. Diagonal preconditioner `P = diag(A)⁻¹` accelerates convergence

**Performance:**

For large problems, warm start achieves dramatic speedups:

| Cutoff | Matrix Size | Dense | LOBPCG (warm start) | Speedup |
|--------|-------------|-------|---------------------|---------|
| 12 | 882×882 | 33 s | 4.5 s | **7.5x** |
| 20 | 2514×2514 | 1055 s | 27 s | **38x** |

**Usage:**

```julia
# Automatic warm start (default)
solver = Solver(PSVWave(), geo, (64, 64), LOBPCGMethod(); cutoff=20)
bands = compute_bands(solver, kpath; bands=1:20)

# Disable warm start (traditional behavior)
method = LOBPCGMethod(warm_start=false, scale=false, first_dense=false)
solver = Solver(PSVWave(), geo, (64, 64), method; cutoff=20)
```

**Manual control with solve_at_k:**

For fine-grained control, use `solve_at_k` directly:

```julia
solver = Solver(PSVWave(), geo, (64, 64); cutoff=20)

# Get matrix dimension for X0
dim = matrix_dimension(solver)  # 2514 for cutoff=20

# First k-point with Dense
freqs1, vecs1 = solve_at_k(solver, k_points[1], DenseMethod();
                           bands=1:20, return_eigenvectors=true)

# Subsequent k-points with warm start
for k in k_points[2:end]
    freqs, vecs = solve_at_k(solver, k, LOBPCGMethod();
                             bands=1:20, X0=vecs, return_eigenvectors=true)
    vecs1 = vecs  # Update for next iteration
end
```

### Limitations

- Warm start works best for photonic problems (TE/TM waves)
- Phononic problems with band crossings at Γ point may require additional band tracking
- For small problems (N < 1000), Dense is often faster due to BLAS optimization

### LOBPCG vs KrylovKit

| Feature | [KrylovKitMethod](api-solver.md#PhoXonic.KrylovKitMethod) | [LOBPCGMethod](api-solver.md#PhoXonic.LOBPCGMethod) |
|---------|-----------------|--------------|
| Backend | [KrylovKit.jl](https://jutho.github.io/KrylovKit.jl/stable/) | [IterativeSolvers.jl](https://iterativesolvers.julialinearalgebra.org/stable/) |
| Algorithm | Arnoldi (Krylov) | Block CG |
| Matrix-free | Yes (2D/3D) | No |
| Phononic scaling | Required | Not required |
| Block computation | No | Yes |
| 3D shift-invert | Yes | Yes |

For details on matrix-free methods and their resolution constraints, see [Matrix-Free Methods](@ref).

## API Reference

- [Solver API](api-solver.md) - DenseMethod, KrylovKitMethod, LOBPCGMethod
