# Fallback errors for extension-dependent functions

function savefig_publication(args...; kwargs...)
    throw(
        ArgumentError(
            "savefig_publication requires PythonPlot.jl. " * "Run `using PythonPlot` first."
        ),
    )
end

function plot_on_axis!(args...; kwargs...)
    throw(
        ArgumentError(
            "plot_on_axis! requires PythonPlot.jl. " * "Run `using PythonPlot` first."
        ),
    )
end

function figure_publication(args...; kwargs...)
    throw(
        ArgumentError(
            "figure_publication requires PythonPlot.jl. " * "Run `using PythonPlot` first."
        ),
    )
end
