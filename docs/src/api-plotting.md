# Plotting

PhoXonic.jl provides two plotting paths via package extensions:

| Extension | Trigger | Functionality |
|-----------|---------|---------------|
| RecipesBaseExt | `using RecipesBase` (or any backend) | `plot(bs)` recipe for `BandStructure` |
| PlotsExt | `using Plots` | `plot_bands`, `plot_bands!` (full-featured) |
| PythonPlotExt | `using PythonPlot` | Planned |

---

## RecipesBase Recipe

The lightweight recipe is activated when any RecipesBase-compatible backend is loaded.

```julia
using PhoXonic, Plots   # or CairoMakie, etc.

solver = Solver(TEWave(), geo, (64, 64))
bs = compute_bands(solver, kpath)

plot(bs)
```

| Attribute | Default |
|-----------|---------|
| xlabel | `"Wave vector"` |
| ylabel | `"Frequency"` |
| legend | `false` |
| grid | `true` |
| linewidth | 2 |

Override any attribute with standard keyword arguments:

```julia
plot(bs, ylabel="ω a / 2πc", title="TE Bands", linewidth=1)
```

---

## Plots.jl Functions

These functions require `using Plots` and provide additional features such as
high-symmetry point tick labels, scatter mode for 3D, and band gap highlighting.

### `plot_bands`

```julia
using PhoXonic, Plots

bs = compute_bands(solver, kpath)
p = plot_bands(bs; title="TE Bands", show_gaps=true)
savefig(p, "bands.png")
```

**Scatter mode** (recommended for 3D to avoid line artifacts at Γ point):

```julia
p = plot_bands(bs; scatter=true, markersize=2)
```

### `plot_bands!`

Add band structure to an existing plot:

```julia
p = plot_bands(bs_te; color=:blue, title="TE + TM")
plot_bands!(p, bs_tm; color=:red)
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `color` | Symbol | `:blue` | Line/marker color |
| `linewidth` | Int | 2 | Line width |
| `scatter` | Bool | `false` | Use scatter instead of lines |
| `markersize` | Int | 3 | Marker size (scatter mode) |
| `title` | String | `"Band Structure"` | Plot title |
| `ylabel` | String | `"Frequency"` | Y-axis label |
| `show_gaps` | Bool | `false` | Highlight band gaps |
| `gap_color` | Symbol | `:yellow` | Gap highlight color |
| `gap_alpha` | Float64 | 0.2 | Gap highlight transparency |
| `normalize` | Float64 | 1.0 | Frequency normalization factor |

### `band_plot_data`

Extract raw plot data from a `BandStructure` (backend-independent):

```julia
data = band_plot_data(bs; normalize=1.0)
# data.distances, data.frequencies, data.label_positions, data.label_names
```

---

## API Reference

```@docs
plot_bands
plot_bands!
band_plot_data
```
