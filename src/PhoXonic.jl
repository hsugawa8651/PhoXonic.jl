# Last-Modified: 2025-12-15T21:30:00+09:00

"""
    PhoXonic.jl

A Julia package for computing band structures of photonic and phononic crystals
using the plane wave expansion (PWE) method.

Supports 1D, 2D, and 3D periodic structures with various wave types:
- Photonic: 1D, TE/TM (2D), TransverseEM (3D)
- Phononic: 1D, SH/P-SV (2D), FullElastic (3D, experimental)
"""
module PhoXonic

using LinearAlgebra
using StaticArrays
using FFTW
using FFTW: plan_fft, plan_ifft
using Brillouin
using LinearMaps
using Krylov
using KrylovKit
using IterativeSolvers

# Core types
include("dimensions.jl")
include("lattice.jl")
include("material.jl")
include("methods.jl")

# Geometry
include("shapes.jl")
include("geometry.jl")

# Plane wave basis
include("planewave.jl")
include("convmat.jl")

# Wave types
include("waves.jl")

# Solvers
include("solver.jl")
include("matrixfree.jl")
include("rscg.jl")

# High-level API
include("kpath.jl")
include("bandstructure.jl")
include("plotting.jl")
include("io.jl")
include("supercell.jl")

# Topological invariants
include("wilson.jl")

# Transfer Matrix Method
include("tmm/TMM.jl")

# Exports - Dimensions
export Dim1, Dim2, Dim3

# Exports - Lattice
export Lattice, lattice_1d, square_lattice, hexagonal_lattice, cubic_lattice, fcc_lattice
export reciprocal_vectors

# Exports - Materials
export Material, PhotonicMaterial, ElasticMaterial
export Dielectric, LossyDielectric, IsotropicElastic, ElasticVoid, from_E_ν
export permittivity, permeability, refractive_index
export density, shear_modulus, transverse_velocity, longitudinal_velocity

# Exports - Shapes
export Circle, Rectangle, Ellipse, Polygon, Sphere, Cylinder, Slab, Segment
export translate  # Shape translation for supercell construction

# Exports - Supercell
export create_supercell, line_defect_positions

# Exports - Geometry
export Geometry, get_material
export DiscretizationMethod, SimpleGrid, SubpixelAverage
export discretize

# Exports - Basis
export PlaneWaveBasis, convolution_matrix

# Exports - Waves
export WaveType, PhotonicWave, PhononicWave
export TEWave, TMWave, SHWave, PSVWave
export Photonic1D, Longitudinal1D
export FullVectorEM, TransverseEM, FullElastic

# Exports - Solver methods
export SolverMethod, IterativeMethod, RSCGMethod
export DenseMethod, BasicRSCG, KrylovKitMethod, LOBPCGMethod

# Exports - Solver
export AbstractSolver, Solver
export solve, solve_at_k, solve_at_k_with_vectors
export build_matrices, get_weight_matrix
export group_velocity, matrix_dimension

# Exports - Matrix-free operators
export FFTContext, MatrixFreeWorkspace
export MatrixFreeOperator, apply_lhs!, apply_rhs!
export to_linear_map_lhs, to_linear_map_rhs

# Exports - Green's function and DOS/LDOS
export compute_greens_function, compute_dos, compute_ldos
export compute_dos_stochastic

# Exports - EffectiveHamiltonian (for advanced users)
export EffectiveHamiltonian, NegatedOperator
export MatrixFreeEffectiveHamiltonian

# Exports - Unified Green's function API
export GFMethod, DirectGF, RSKGF, MatrixFreeGF
export RHSInvMethod, ApproximateRHSInv, CGRHSInv

# Exports - K-path
export SimpleKPath, simple_kpath_square, simple_kpath_hexagonal
export simple_kpath_cubic, simple_kpath_fcc, simple_kpath_bcc, kpath_fcc_joannopoulos
export kpath_from_brillouin, kpath_square, kpath_hexagonal
export kpath_cubic, kpath_fcc, kpath_bcc

# Exports - Band structure
export BandStructure, compute_bands
export find_bandgap, find_all_gaps
export frequencies, distances, labels, nbands, nkpoints

# Exports - Plotting (requires Plots.jl)
export plot_bands, plot_bands!, band_plot_data

# Exports - I/O (requires JLD2.jl)
export save_bands, load_bands
export save_modes, load_modes
export save_epsilon, load_epsilon

# Exports - TMM (Transfer Matrix Method)
export Layer, Multilayer
export thickness, material, layers, incident, substrate
export nlayers, total_thickness
export periodic_multilayer, bragg_mirror
export TMMSolver, TMMResult, tmm_spectrum, tmm_bandstructure, bloch_k

# Exports - Wilson loop / Topological invariants
export ZakPhaseResult, compute_zak_phase
export WilsonSpectrumResult, compute_wilson_spectrum, winding_number

end # module
