# Tests for band sorting functionality

@testset "Band Sorting" begin

    @testset "find_best_permutation" begin
        # Identity case
        M = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
        @test PhoXonic.find_best_permutation(M) == [1, 2, 3]

        # Permuted case
        M = [0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0]
        @test PhoXonic.find_best_permutation(M) == [2, 3, 1]

        # Soft overlap (realistic case)
        M = [0.9 0.1 0.0; 0.1 0.8 0.1; 0.0 0.1 0.9]
        @test PhoXonic.find_best_permutation(M) == [1, 2, 3]

        # Another permutation
        M = [0.1 0.9 0.0; 0.9 0.1 0.0; 0.0 0.0 1.0]
        @test PhoXonic.find_best_permutation(M) == [2, 1, 3]
    end

    @testset "compute_bands with track_bands (2D)" begin
        # 2D TM with potential band crossing
        lat = square_lattice(1.0)
        geo = Geometry(lat, Dielectric(1.0), [
            (Circle([0.5, 0.5], 0.2), Dielectric(9.0))
        ])
        solver = Solver(TMWave(), geo, (32, 32); cutoff=5)

        # Simple k-path
        kpath = simple_kpath_square(; npoints=11)

        # Without tracking
        bands_unsorted = compute_bands(solver, kpath; bands=1:5, track_bands=false)

        # With tracking
        bands_sorted = compute_bands(solver, kpath; bands=1:5, track_bands=true)

        # Same frequencies at endpoints (Γ point)
        @test bands_unsorted.frequencies[1, :] ≈ bands_sorted.frequencies[1, :]

        # Sorted should be smoother (lower total variation)
        function total_variation(freqs)
            sum(abs.(diff(freqs, dims=1)))
        end
        @test total_variation(bands_sorted.frequencies) <= total_variation(bands_unsorted.frequencies) + 0.01
    end

    @testset "track_bands preserves eigenvalues" begin
        # Ensure tracking doesn't change the set of eigenvalues, just their order
        lat = square_lattice(1.0)
        geo = Geometry(lat, Dielectric(1.0), [
            (Circle([0.5, 0.5], 0.3), Dielectric(4.0))
        ])
        solver = Solver(TEWave(), geo, (24, 24); cutoff=4)

        kpath = simple_kpath_square(; npoints=5)

        bands_unsorted = compute_bands(solver, kpath; bands=1:4, track_bands=false)
        bands_sorted = compute_bands(solver, kpath; bands=1:4, track_bands=true)

        # At each k-point, the set of frequencies should be the same (just reordered)
        for ik in 1:5
            @test sort(bands_unsorted.frequencies[ik, :]) ≈ sort(bands_sorted.frequencies[ik, :])
        end
    end
end
