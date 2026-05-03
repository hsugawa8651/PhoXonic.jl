using Aqua
using PhoXonic

@testset "Aqua quality" begin
    Aqua.test_ambiguities(PhoXonic)
    Aqua.test_unbound_args(PhoXonic)
    Aqua.test_undefined_exports(PhoXonic)
    Aqua.test_project_extras(PhoXonic)
    Aqua.test_stale_deps(PhoXonic)
    Aqua.test_deps_compat(PhoXonic)
end
