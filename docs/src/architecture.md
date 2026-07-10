# Architecture

This page describes how the solver types fit together in v0.3.0, and which
analysis functions exist in which dimension. It documents what is in `main`;
the last section says what is not.

## The (dimension, wave type) lattice

A wave type belongs to exactly one dimension. Dimension and wave type are not
independent axes, so the valid solvers are not the free product of the two.
There are nine:

| Dimension | Photonic | Phononic |
|-----------|----------|----------|
| 1D | [`Photonic1D`](api-solver.md#PhoXonic.Photonic1D) | [`Longitudinal1D`](api-solver.md#PhoXonic.Longitudinal1D) |
| 2D | [`TEWave`](api-solver.md#PhoXonic.TEWave), [`TMWave`](api-solver.md#PhoXonic.TMWave) | [`SHWave`](api-solver.md#PhoXonic.SHWave), [`PSVWave`](api-solver.md#PhoXonic.PSVWave) |
| 3D | [`TransverseEM`](api-solver.md#PhoXonic.TransverseEM), [`FullVectorEM`](api-solver.md#PhoXonic.FullVectorEM) | [`FullElastic`](api-solver.md#PhoXonic.FullElastic) |

A combination outside this table is rejected when the solver is constructed, not
later when the operator is assembled:

```julia
julia> lat = lattice_1d(1.0);
julia> geo = Geometry(lat, Dielectric(1.0), [(Segment(0.2, 0.8), Dielectric(9.0))]);
julia> Solver(TEWave(), geo, (64,); cutoff=5)
ERROR: ArgumentError: Unsupported wave type / dimension combination: TEWave with Dim1
```

[`TEWave`](api-solver.md#PhoXonic.TEWave) is a two-dimensional polarization, so
there is nothing for it to solve on a one-dimensional lattice. The same holds for
every combination absent from the table.

## Dimension support of the analysis functions

Band structures are available in every dimension. The post-analysis functions are
not, because the quantities they compute do not all exist in every dimension, and
because some of them are not implemented yet.

| Function | 1D | 2D | 3D |
|----------|:--:|:--:|----|
| [`solve`](api-solver.md#PhoXonic.solve), [`compute_bands`](api-solver.md#PhoXonic.compute_bands) | ✓ | ✓ | ✓ |
| [`compute_greens_function`](api-advanced.md#PhoXonic.compute_greens_function) | ✓ | ✓ | ✓ |
| [`compute_dos`](api-advanced.md#PhoXonic.compute_dos) | ✓ | ✓ | — |
| [`compute_ldos`](api-advanced.md#PhoXonic.compute_ldos) | ✓ | ✓ | [`MatrixFreeGF()`](api-advanced.md#PhoXonic.MatrixFreeGF) only |
| [`compute_zak_phase`](api-advanced.md#PhoXonic.compute_zak_phase) | ✓ | — | — |
| [`compute_wilson_spectrum`](api-advanced.md#PhoXonic.compute_wilson_spectrum) | — | ✓ | — |

A dash means the quantity is unavailable, but for two different reasons. The Zak
phase is a one-dimensional invariant and the Wilson loop spectrum a
two-dimensional one, so outside their dimension there is nothing to compute;
calling them there is a `MethodError`, never a wrong number. A three-dimensional
[`compute_dos`](api-advanced.md#PhoXonic.compute_dos), by contrast, would be
meaningful and is simply not implemented.

The one conditional cell is
[`compute_ldos`](api-advanced.md#PhoXonic.compute_ldos) in 3D: it needs
[`MatrixFreeGF()`](api-advanced.md#PhoXonic.MatrixFreeGF), from the
`ReducedShiftedKrylov` extension, and a wave type whose right-hand side can be
inverted ([`FullVectorEM`](api-solver.md#PhoXonic.FullVectorEM) or
[`FullElastic`](api-solver.md#PhoXonic.FullElastic), not
[`TransverseEM`](api-solver.md#PhoXonic.TransverseEM)). See
[DOS / LDOS](greens_function.md#Dimension-Support) for which `GFMethod` each
function accepts in each dimension.

## Two stages

Every calculation in PhoXonic.jl passes through the same two stages.

**Stage 1, dispersion.** A [`Solver`](api-solver.md#PhoXonic.Solver) assembles the
eigenvalue problem for a wave type on a geometry, and
[`solve`](api-solver.md#PhoXonic.solve) or
[`compute_bands`](api-solver.md#PhoXonic.compute_bands) returns the frequencies
ω_n(k) and the modes.

**Stage 2, post-analysis.** Everything else consumes the output of stage 1, or
the operator behind it: the density of states, the local density of states,
Green's functions, the Zak phase, the Wilson loop spectrum, field reconstruction.

The two stages are why the dimension support above is uneven. Stage 1 works
wherever a `(dimension, wave type)` pair exists, that is, everywhere. Stage 2
depends on the quantity: some do not exist in a given dimension, and some are
simply not implemented there.

## Not in v0.3.0

Three things a reader may expect to find here, and where they actually stand as
of v0.3.0.

**A unified `solve(problem, algorithm)` entry point.** The SciML idiom of
describing a problem and then choosing an algorithm to solve it is a design, not
code in `main`. The exported [`solve`](api-solver.md#PhoXonic.solve) is the
eigenvalue problem at one wave vector, and is unrelated to it. Stage 2 is reached
through the individual `compute_*` functions.

**Berry curvature and valley analysis.** Designed, not in `main`.

**A three-dimensional density of states.** Not implemented. Green's functions in
3D are available; the density of states built on them is not.
