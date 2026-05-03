@testset "PhoXonic.jl (ext)" begin
    # PythonPlot extension: savefig_publication
    include("test_publication.jl")

    # ReducedShiftedKrylov extension
    include("rsk_ext/runtests.jl")
end
