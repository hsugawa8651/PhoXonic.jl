@testset "PhoXonic.jl (ext)" begin
    # PythonPlot extension: savefig_publication
    include("test_publication.jl")

    # ReducedShiftedKrylov extension
    include("rsk_ext/runtests.jl")

    # Plots extension (RecipesBase ext is covered transitively via Plots)
    include("plotsext_test.jl")
end
