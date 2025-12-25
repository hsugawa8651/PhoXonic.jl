# PhoXonic.jl

[![CI](https://github.com/hsugawa8651/PhoXonic.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/hsugawa8651/PhoXonic.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/hsugawa8651/PhoXonic.jl/actions/workflows/Documentation.yml/badge.svg)](https://hsugawa8651.github.io/PhoXonic.jl/dev/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18055242.svg)](https://doi.org/10.5281/zenodo.18055242)

Band structure calculation for photonic and phononic crystals using the plane wave expansion (PWE) method.

![Phoxonic crystal band structure](paper/214_maldovan2006_phoxonic.png)

*Phoxonic crystal (Si with air holes): photonic (TE, TM) and phononic (SH, P-SV) band structures. See [example code](examples/214_maldovan2006_phoxonic.jl).*

## Features

- **Photonic crystals**: 1D, TE/TM (2D), TransverseEM (3D)
- **Phononic crystals**: 1D, SH/P-SV (2D), FullElastic (3D, experimental)
- **Solvers**: Dense, KrylovKit, LOBPCG
- **Analysis**: Band gaps, group velocity, LDOS

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/hsugawa8651/PhoXonic.jl")
```

## Quick Start

See the [Getting Started guide](https://hsugawa8651.github.io/PhoXonic.jl/dev/getting_started/) for usage examples.

## Documentation

See the [full documentation](https://hsugawa8651.github.io/PhoXonic.jl/dev/) for workflow guides, examples, and API reference.

## Contributing

Contributions are welcome! Please open an issue or pull request on [GitHub](https://github.com/hsugawa8651/PhoXonic.jl).

## License

[MIT License](LICENSE)
