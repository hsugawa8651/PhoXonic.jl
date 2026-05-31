using Documenter
using PhoXonic
using PythonPlot                  # loads PhoXonicPythonPlotExt so its docstrings render
PythonPlot.matplotlib.use("Agg")  # headless backend for CI

# The PythonPlot extension carries the per-method docstrings for plot_on_axis! /
# figure_publication; include it so @docs can render them.
const PXPythonPlotExt = Base.get_extension(PhoXonic, :PhoXonicPythonPlotExt)

makedocs(;
    sitename="PhoXonic.jl",
    modules=[PhoXonic, PXPythonPlotExt],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Workflow (2D)" => "workflow.md",
        "1D Calculations" => "workflow_1d.md",
        "3D Calculations" => "workflow_3d.md",
        "Solver Methods" => "solver.md",
        "Analysis" => "analysis.md",
        "Field Visualization" => "field_visualization.md",
        "DOS / LDOS" => "greens_function.md",
        "Topological Invariants" => "topology.md",
        "Examples" => "examples.md",
        "Matrix-Free Methods" => "matrixfree.md",
        "Elastic Void (Tanaka Limit)" => "elastic_void.md",
        "Transfer Matrix Method (1D)" => "tmm.md",
        "API Reference" => [
            "Core" => "api.md",
            "Solver" => "api-solver.md",
            "Advanced" => "api-advanced.md",
            "Plotting" => "api-plotting.md",
            "I/O" => "api-io.md",
        ],
        "Dependencies" => "dependencies.md",
    ],
    # Source links always point to public repo
    repo=Documenter.Remotes.GitHub("hsugawa8651", "PhoXonic.jl"),
    warnonly=[:missing_docs, :cross_references],
)

# Deploy to GitHub Pages
deploydocs(;
    repo="github.com/hsugawa8651/PhoXonic.jl.git", devbranch="main", push_preview=true
)
