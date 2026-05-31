# Plotting

PhoXonic.jl provides three plotting paths via package extensions:

| Extension | Trigger | Functionality |
|-----------|---------|---------------|
| RecipesBaseExt | `using RecipesBase` (or any backend) | `plot(bs)` recipe for `BandStructure` |
| PlotsExt | `using Plots` | `plot_bands`, `plot_bands!` (full-featured) |
| PythonPlotExt | `using PythonPlot` | 3-layer API `plot_on_axis!` / `figure_publication` / `savefig_publication` for publication-quality static PDF / PNG |

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

Activated when `using PythonPlot` is loaded. Provides matplotlib-backed PDF / PNG
output with precise axis-size control (in millimeters) suitable for journal
submissions. The API has **three layers** — L1 / L2 are unexported (call them as
`PhoXonic.…`); L3 `savefig_publication` is exported:

| Layer | Function | Returns | Use |
|-------|----------|---------|-----|
| L3 | `savefig_publication(bs, path; ...)` | `path` | save to PDF or PNG in one call (format from the file extension; delegates to L2) |
| L2 | `PhoXonic.figure_publication(bs; ...)` | `(fig, ax)` | a sized figure + axis to tweak before saving |
| L1 | `PhoXonic.plot_on_axis!(ax, bs; ...)` | `ax` | draw onto your own matplotlib axis (compose a subplot grid) |

L2 and L3 (single `BandStructure`) take `axis_width_mm` / `axis_height_mm`; the
vector `savefig_publication` overload is **L3-only** and additionally takes
`layout` / `overlay` (see [Multiple band structures](@ref) below).

**Headless / CI.** matplotlib needs a non-interactive backend when there is no
display. Either set the environment variable *before* `using PythonPlot`

```julia
ENV["MPLBACKEND"] = "Agg"
using PhoXonic, PythonPlot
```

or select it explicitly after import: `PythonPlot.matplotlib.use("Agg")`.

!!! note "PhoXonic specifics"
    - **Millimeters only.** Axis sizes are given in mm (`axis_width_mm`,
      `axis_height_mm`); there is no `cm` keyword.
    - **`xlabel` is empty by default.** The horizontal axis is labelled with
      high-symmetry-point xticks (Γ, X, M, …), so no axis title is drawn unless
      you pass a non-empty `xlabel`.
    - **No colorbar.** Band structures are line plots; the colorbar-related
      keywords found in other PhoXonic ecosystem packages do not apply here.

### `savefig_publication` (L3)

Single `BandStructure`:

```julia
using PhoXonic, PythonPlot

bs = compute_bands(solver, kpath)
savefig_publication(bs, "bands.pdf"; axis_width_mm=80.0, axis_height_mm=60.0,
                    title="TE Bands")
```

#### Multiple band structures

A `Vector` of `BandStructure` is supported **only at L3** (`savefig_publication`);
there is no `figure_publication` overload for vectors. Two modes are available:

*Grid* — one panel per band structure (`layout = (nrows, ncols)`):

```julia
savefig_publication([bs_te, bs_tm], "te_tm.pdf"; layout=(1, 2))
```

*Overlay* — all band structures on a single axis, cycling through preset line
styles (`overlay = true`):

```julia
savefig_publication([bs_te, bs_tm], "te_tm_overlay.pdf"; overlay=true,
                    title="TE vs TM")
```

Plot keywords (`show_gaps`, `normalize`, `ylabel`, …) are forwarded to each panel
in both modes; in overlay mode an explicit `color` / `linewidth` / `linestyle`
overrides the preset style for every curve. To compose vectors more freely (e.g.
mixed layouts), use [`plot_on_axis!`](@ref PhoXonic.plot_on_axis!) with your own
`subplots()`.

### Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `axis_width_mm` | Real | `80.0` | Axis width in millimeters |
| `axis_height_mm` | Real | `60.0` | Axis height in millimeters |
| `ylims` | Tuple or `nothing` | `nothing` | Y-axis range (auto if `nothing`) |
| `title` | String | `""` | Plot title (suppressed if empty) |
| `show_gaps` | Bool | `false` | Highlight band gaps with `axhspan` |
| `normalize` | Real | `1.0` | Frequency normalization factor |
| `layout` | Tuple | `(1, n)` | Subplot grid `(nrows, ncols)` |
| `overlay` | Bool | `false` | Plot all on a single axis |

The output format is inferred from the file extension (`.pdf`, `.png`, `.svg`,
etc., as supported by matplotlib).

### Composing subplots with `plot_on_axis!` (L1) / `figure_publication` (L2)

Beyond `savefig_publication` (L3) shown above, the two lower layers let you
compose and tweak figures.

**L1 — compose multiple band structures into one figure.** You create the
figure with matplotlib's `subplots`/`add_subplot` and own the `close`:

```julia
using PhoXonic, PythonPlot

fig = PythonPlot.figure(figsize = (10, 4))
PhoXonic.plot_on_axis!(fig.add_subplot(1, 2, 1), bs_te; title = "TE")
PhoXonic.plot_on_axis!(fig.add_subplot(1, 2, 2), bs_tm; title = "TM")
fig.savefig("te_tm_panel.pdf")
PythonPlot.close(fig)            # caller owns the figure
```

The plot keyword arguments are the same as for `savefig_publication`
(e.g. `title`, `color`, `linewidth`, `linestyle`, `show_gaps`, `normalize`,
`ylabel`). PhoXonic uses high-symmetry-point xticks for the horizontal axis, so
`xlabel` is empty by default; pass a non-empty `xlabel` to add one.
`plot_on_axis!` returns the same `ax` for chaining.

**L2 — adjust before saving:**

```julia
fig, ax = PhoXonic.figure_publication(bs; axis_width_mm = 100.0)
ax.axhline(0.5; linestyle = "--")     # add your own annotations
fig.savefig("annotated.pdf")
PythonPlot.close(fig)
```

!!! note
    Inside the extension, figures are created with `PythonPlot.figure(...)` (not
    `subplots()`), but on the **caller** side `subplots()` is fine — you just own
    the resulting figure and must `close` it yourself.

On a headless machine or in CI (no display), select a non-interactive backend
before plotting: `PythonPlot.matplotlib.use("Agg")`.

### Installing PythonPlot

The extension loads automatically once both `PhoXonic` and `PythonPlot` are
imported in the same session. Install it with:

```julia
using Pkg
Pkg.add("PythonPlot")
```

CondaPkg.jl will install matplotlib on first use.

---

## API Reference

```@docs
plot_bands
plot_bands!
band_plot_data
savefig_publication
PhoXonic.figure_publication
PhoXonic.plot_on_axis!
```
