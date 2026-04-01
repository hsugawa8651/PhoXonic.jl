# Fallback errors for extension-dependent functions

function savefig_publication(args...; kwargs...)
    throw(ArgumentError(
        "savefig_publication requires PythonPlot.jl. " *
        "Run `using PythonPlot` first."
    ))
end
