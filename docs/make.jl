using Documenter
using PhoXonic

makedocs(
    sitename = "PhoXonic.jl",
    modules = [PhoXonic],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = [
        "Home" => "index.md",
        "Workflow" => "workflow.md",
        "3D Calculations" => "workflow_3d.md",
        "Solver Methods" => "solver.md",
        "Analysis" => "analysis.md",
        "DOS / LDOS" => "greens_function.md",
        "Examples" => "examples.md",
        "Matrix-Free Methods" => "matrixfree.md",
        "Elastic Void (Tanaka Limit)" => "elastic_void.md",
        "API Reference" => "api.md",
    ],
    remotes = nothing,
    warnonly = [:missing_docs],
)

# Deploy to GitHub Pages
deploydocs(
    repo = "github.com/hsugawa8651/PhoXonic.jl.git",
    devbranch = "main",
    push_preview = true,
)
