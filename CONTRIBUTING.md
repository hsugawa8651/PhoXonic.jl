# Contributing to PhoXonic.jl

## Running tests

PhoXonic.jl uses a split test runner:

- `Pkg.test()` — runs all tests (default, ~3 min)
- `Pkg.test(; test_args=["core"])` — core tests only (~1.5 min)
- `Pkg.test(; test_args=["ext"])` — extension tests only (~1 min)

CI runs both groups in parallel via `test_group: [core, ext]` matrix.

## Extensions

| Extension | Trigger package | Test file |
|---|---|---|
| `PhoXonicPlotsExt` | `Plots`, `RecipesBase` | `test/plotsext_test.jl` |
| `PhoXonicPythonPlotExt` | `PythonPlot` | `test/test_publication.jl` |
| `PhoXonicRecipesBaseExt` | `RecipesBase` | (covered transitively by `Plots` test) |
| `PhoXonicReducedShiftedKrylovExt` | `ReducedShiftedKrylov` | `test/rsk_ext/runtests.jl` |

`test/ext_tests.jl` includes the test files listed above.

## Adding a new extension test

1. Add the trigger package to `test/Project.toml`.
2. Create `test/<name>ext_test.jl`.
3. Add `include("<name>ext_test.jl")` to `test/ext_tests.jl`.

## Code style

JuliaFormatter v2 (Blue style via `.JuliaFormatter.toml`) is enforced
by the `Format Check` workflow. Run `JuliaFormatter.format(".")`
before commit.

## Quality checks

`Aqua.jl` (in core test group) checks for ambiguities, undefined
exports, stale deps, and missing compat entries.

## Release

Maintainer workflow:

1. Bump version in `Project.toml`, `CITATION.cff`, `README.md`,
   and `docs/src/index.md`.
2. Open a PR titled `chore: bump version to vX.Y.Z`.
3. After merge, comment `@JuliaRegistrator register` on the merge
   commit (see [General registry](https://github.com/JuliaRegistries/General)).
4. TagBot creates the git tag and GitHub Release. Zenodo issues a
   per-version DOI.
