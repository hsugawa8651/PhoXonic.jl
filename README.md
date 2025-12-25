# PhoXonic.jl

[![CI](https://github.com/phoxonic-dev/PhoXonic.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/phoxonic-dev/PhoXonic.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/phoxonic-dev/PhoXonic.jl/actions/workflows/Documentation.yml/badge.svg)](https://phoxonic-dev.github.io/PhoXonic.jl/dev/)

Band structure calculation for photonic and phononic crystals using the plane wave expansion (PWE) method.

## Features

- **Photonic crystals**: 1D, TE/TM (2D), TransverseEM (3D)
- **Phononic crystals**: 1D, SH/P-SV (2D), FullElastic (3D, experimental)
- **Solvers**: Dense, KrylovKit, LOBPCG
- **Analysis**: Band gaps, group velocity, LDOS

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/phoxonic-dev/PhoXonic.jl")
```

## Documentation

See the [full documentation](https://phoxonic-dev.github.io/PhoXonic.jl/dev/) for usage, examples, and API reference.

## License

MIT License
