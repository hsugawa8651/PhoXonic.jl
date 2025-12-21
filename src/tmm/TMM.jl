# Transfer Matrix Method module for 1D structures
#
# Provides transmission/reflection spectrum calculations
# for photonic and phononic multilayer structures.

# Layer types
include("layers.jl")

# Transfer matrix calculations
include("matrices.jl")

# TMM solver
include("solver.jl")

# Band structure calculation
include("bands.jl")
