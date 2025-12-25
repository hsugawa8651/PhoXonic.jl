# Dependencies

PhoXonic.jl uses the following packages:

## Core Dependencies

| Package | Description |
|---------|-------------|
| [Brillouin.jl](https://github.com/thchr/Brillouin.jl) | Brillouin zone paths and symmetry |
| [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) | FFT for matrix-free operators (O(N log N) convolution) |
| [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) | Fixed-size arrays for vectors |
| [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) | Standard linear algebra (stdlib) |
| [LinearMaps.jl](https://github.com/JuliaLinearAlgebra/LinearMaps.jl) | Lazy linear operators |

## Iterative Solvers

| Package | Description |
|---------|-------------|
| [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl) | Krylov subspace methods (Arnoldi, Lanczos) |
| [Krylov.jl](https://github.com/JuliaSmoothOptimizers/Krylov.jl) | CG, MINRES, and other Krylov methods |
| [IterativeSolvers.jl](https://github.com/JuliaLinearAlgebra/IterativeSolvers.jl) | LOBPCG eigensolver |

## Optional Extensions

| Package | Extension | Description |
|---------|-----------|-------------|
| [Plots.jl](https://github.com/JuliaPlots/Plots.jl) | PhoXonicPlotsExt | Band structure plotting |
| [ReducedShiftedKrylov.jl](https://github.com/hsugawa8651/ReducedShiftedKrylov.jl) | PhoXonicReducedShiftedKrylovExt | Matrix-free Green's function |

## I/O

| Package | Description |
|---------|-------------|
| [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) | Save/load band structures |
