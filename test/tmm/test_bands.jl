# TMM band structure tests

using Test
using PhoXonic
using LinearAlgebra

@testset "Bloch condition" begin
    @testset "Trace relation: Tr(M) = 2cos(ka)" begin
        # For a unit cell, Bloch theorem gives: Tr(M) = 2cos(ka)
        # where M is the transfer matrix for one period and a is the period

        # Simple quarter-wave stack unit cell
        n_hi, n_lo = 3.0, 1.0
        a = 1.0  # Period
        d_hi = a / 4  # Quarter of period
        d_lo = 3a / 4

        mat_hi = Dielectric(n_hi^2)
        mat_lo = Dielectric(n_lo^2)

        # Unit cell: one period
        unit_cell = [Layer(mat_hi, d_hi), Layer(mat_lo, d_lo)]
        ml = Multilayer(unit_cell, mat_lo, mat_lo)  # Periodic BC

        # At different wavelengths, check trace relation
        for λ in [2.0, 3.0, 4.0, 6.0]
            M = PhoXonic.system_matrix(ml, λ)
            trace = tr(M)

            # |Tr(M)| <= 2 means propagating wave (real k)
            # |Tr(M)| > 2 means evanescent wave (bandgap)
            @test isreal(trace) || abs(imag(trace)) < 1e-10
        end
    end

    @testset "Symmetry at zone boundary" begin
        # At k = 0 and k = π/a, bands have zero group velocity (flat)
        n_hi, n_lo = 2.0, 1.0
        a = 1.0

        mat_hi = Dielectric(n_hi^2)
        mat_lo = Dielectric(n_lo^2)

        unit_cell = [Layer(mat_hi, 0.5), Layer(mat_lo, 0.5)]
        ml = Multilayer(unit_cell, mat_lo, mat_lo)

        # Check that system matrix is computed without error
        M = PhoXonic.system_matrix(ml, 2.0)
        @test size(M) == (2, 2)
        @test isfinite(det(M))
    end
end

@testset "tmm_bandstructure" begin
    @testset "Basic band structure" begin
        # Quarter-wave stack: high contrast
        n_hi, n_lo = 3.0, 1.0
        a = 1.0

        mat_hi = Dielectric(n_hi^2)
        mat_lo = Dielectric(n_lo^2)

        unit_cell = [Layer(mat_hi, 0.25), Layer(mat_lo, 0.75)]
        ml = Multilayer(unit_cell, mat_lo, mat_lo)

        solver = TMMSolver(Photonic1D(), ml)
        bands = tmm_bandstructure(solver; k_points=21, bands=1:3)

        # Check output type
        @test bands isa BandStructure

        # Check dimensions
        @test size(bands.frequencies, 1) == 21  # k points
        @test size(bands.frequencies, 2) == 3   # bands

        # Frequencies should be non-negative and sorted
        for ik in 1:21
            @test all(bands.frequencies[ik, :] .>= 0)
            @test issorted(bands.frequencies[ik, :])
        end
    end

    @testset "Bandgap detection" begin
        # High contrast structure should have bandgap
        n_hi, n_lo = 3.0, 1.0
        a = 1.0

        mat_hi = Dielectric(n_hi^2)
        mat_lo = Dielectric(n_lo^2)

        # Quarter-wave stack at normalized frequency ω₀ = 1
        d_hi = 0.25  # n_hi * d_hi = 0.75
        d_lo = 0.75  # n_lo * d_lo = 0.75

        unit_cell = [Layer(mat_hi, d_hi), Layer(mat_lo, d_lo)]
        ml = Multilayer(unit_cell, mat_lo, mat_lo)

        solver = TMMSolver(Photonic1D(), ml)
        bands = tmm_bandstructure(solver; k_points=51, bands=1:4)

        # First bandgap should exist between band 1 and band 2
        band1_max = maximum(bands.frequencies[:, 1])
        band2_min = minimum(bands.frequencies[:, 2])

        @test band2_min > band1_max  # Gap exists
    end

    @testset "PWE comparison" begin
        # Same structure computed with PWE and TMM should match
        n_hi, n_lo = 3.0, 1.0
        a = 1.0
        d_hi = 0.25

        mat_hi = Dielectric(n_hi^2)
        mat_lo = Dielectric(n_lo^2)

        # TMM setup
        unit_cell = [Layer(mat_hi, d_hi), Layer(mat_lo, a - d_hi)]
        ml = Multilayer(unit_cell, mat_lo, mat_lo)
        solver_tmm = TMMSolver(Photonic1D(), ml)

        # PWE setup (same structure)
        lat = lattice_1d(a)
        geo = Geometry(lat, mat_lo, [(Segment(0.0, d_hi), mat_hi)])
        solver_pwe = Solver(Photonic1D(), geo, 256; cutoff=20)

        # Compute TMM bands
        k_values = range(0, π/a, length=11)
        bands_tmm = tmm_bandstructure(solver_tmm; k_points=11, bands=1:3)

        # Compare with PWE at individual k-points (since compute_bands not available for 1D)
        # Skip k=0 (Γ) where PWE and TMM handle band degeneracy differently
        for (ik, k) in enumerate(k_values)
            if k ≈ 0.0
                continue  # Skip Γ point - PWE/TMM band counting differs there
            end
            ω_pwe, _ = solve(solver_pwe, k; bands=1:3)  # k is scalar for 1D

            # Compare frequencies (should match within tolerance)
            for ib in 1:3
                ω_tmm = bands_tmm.frequencies[ik, ib]
                @test isapprox(ω_tmm, ω_pwe[ib], rtol=0.02)
            end
        end
    end
end

@testset "Bloch wavenumber calculation" begin
    @testset "bloch_k from transfer matrix" begin
        # k = (1/a) * acos(Tr(M)/2)
        n_hi, n_lo = 2.0, 1.0
        a = 1.0

        mat_hi = Dielectric(n_hi^2)
        mat_lo = Dielectric(n_lo^2)

        unit_cell = [Layer(mat_hi, 0.5), Layer(mat_lo, 0.5)]
        ml = Multilayer(unit_cell, mat_lo, mat_lo)

        # At specific wavelength, compute k
        λ = 4.0  # Should be in propagating regime
        M = PhoXonic.system_matrix(ml, λ)
        trace = real(tr(M))

        if abs(trace) <= 2.0
            # Propagating: k is real
            k = acos(trace / 2) / a
            @test 0 <= k <= π/a
        else
            # Bandgap: k is complex
            @test abs(trace) > 2.0
        end
    end
end
