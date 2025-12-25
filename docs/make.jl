using Documenter
using PhoXonic

makedocs(;
    sitename="PhoXonic.jl",
    modules=[PhoXonic],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Workflow (2D)" => "workflow.md",
        "1D Calculations" => "workflow_1d.md",
        "3D Calculations" => "workflow_3d.md",
        "Solver Methods" => "solver.md",
        "Analysis" => "analysis.md",
        "DOS / LDOS" => "greens_function.md",
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
