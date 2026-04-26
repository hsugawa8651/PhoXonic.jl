# Plotting

PhoXonic.jl provides two plotting paths via package extensions:

| Extension | Trigger | Functionality |
|-----------|---------|---------------|
| RecipesBaseExt | `using RecipesBase` (or any backend) | `plot(bs)` recipe for `BandStructure` |
| PlotsExt | `using Plots` | `plot_bands`, `plot_bands!` (full-featured) |
| PythonPlotExt | `using PythonPlot` | `savefig_publication` (publication-quality PDF/PNG) |

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

## PythonPlot — Publication-Quality Figures

Activated when `using PythonPlot` is loaded. Provides `savefig_publication`
for matplotlib-backed PDF/PNG output with precise axis-size control suitable
for journal submissions.

### `savefig_publication`

Single `BandStructure`:

```julia
using PhoXonic, PythonPlot

bs = compute_bands(solver, kpath)
savefig_publication(bs, "bands.pdf"; axis_width_cm=8.0, axis_height_cm=6.0,
                    title="TE Bands")
```

Multiple `BandStructure` as subplots:

```julia
savefig_publication([bs_te, bs_tm], "te_tm.pdf"; layout=(1, 2))
```

Multiple `BandStructure` as overlay (single axis):

```julia
savefig_publication([bs_te, bs_tm], "te_tm_overlay.pdf"; overlay=true,
                    title="TE vs TM")
```

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `axis_width_cm` | Real | `8.0` | Axis width in centimeters |
| `axis_height_cm` | Real | `6.0` | Axis height in centimeters |
| `ylims` | Tuple or `nothing` | `nothing` | Y-axis range (auto if `nothing`) |
| `title` | String | `""` | Plot title (suppressed if empty) |
| `show_gaps` | Bool | `false` | Highlight band gaps with `axhspan` |
| `normalize` | Real | `1.0` | Frequency normalization factor |
| `layout` | Tuple | `(1, n)` | Subplot grid `(nrows, ncols)` |
| `overlay` | Bool | `false` | Plot all on a single axis |

The output format is inferred from the file extension (`.pdf`, `.png`, `.svg`,
etc., as supported by matplotlib).

---

## API Reference

```@docs
plot_bands
plot_bands!
band_plot_data
savefig_publication
```
