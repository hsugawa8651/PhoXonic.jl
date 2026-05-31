using Test
using PhoXonic
using PythonPlot
PythonPlot.matplotlib.use("Agg")
using StaticArrays

@testset "savefig_publication (PythonPlot)" begin
    # Lightweight 2D BandStructure fixture (4 k-points × 3 bands)
    bs2 = PhoXonic.BandStructure{2}(
        [SVector(0.0, 0.0), SVector(0.5, 0.0), SVector(0.5, 0.5), SVector(0.0, 0.0)],
        [0.0, 1.0, 2.0, 3.0],
        Float64[
            0.10 0.30 0.55;
            0.20 0.40 0.65;
            0.15 0.35 0.60;
            0.05 0.25 0.50
        ],
        [(1, "Γ"), (2, "X"), (3, "M"), (4, "Γ")],
    )

    # Lightweight 3D BandStructure fixture
    bs3 = PhoXonic.BandStructure{3}(
        [
            SVector(0.0, 0.0, 0.0),
            SVector(0.5, 0.0, 0.0),
            SVector(0.5, 0.5, 0.0),
            SVector(0.5, 0.5, 0.5),
            SVector(0.0, 0.0, 0.0),
        ],
        [0.0, 1.0, 2.0, 3.0, 4.0],
        Float64[
            0.10 0.30 0.55;
            0.20 0.40 0.65;
            0.15 0.35 0.60;
            0.18 0.38 0.62;
            0.05 0.25 0.50
        ],
        [(1, "Γ"), (2, "X"), (3, "M"), (4, "R"), (5, "Γ")],
    )

    @testset "BandStructure{2} 単体" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bs2.pdf")
            ret = savefig_publication(bs2, path)
            @test isfile(path)
            @test ret == path
            @test filesize(path) > 0
        end
    end

    @testset "BandStructure{3} 単体" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bs3.pdf")
            savefig_publication(bs3, path)
            @test isfile(path)
        end
    end

    @testset "axis_width_mm / axis_height_mm カスタム" begin
        mktempdir() do tmp
            path = joinpath(tmp, "custom_size.pdf")
            savefig_publication(bs2, path; axis_width_mm=100.0, axis_height_mm=50.0)
            @test isfile(path)
        end
    end

    @testset "title kwarg" begin
        mktempdir() do tmp
            path = joinpath(tmp, "with_title.pdf")
            savefig_publication(bs2, path; title="Test Bands")
            @test isfile(path)
        end
    end

    @testset "ylims kwarg" begin
        mktempdir() do tmp
            path = joinpath(tmp, "ylims.pdf")
            savefig_publication(bs2, path; ylims=(0.0, 1.0))
            @test isfile(path)
        end
    end

    @testset "show_gaps" begin
        mktempdir() do tmp
            path = joinpath(tmp, "gaps.pdf")
            savefig_publication(bs2, path; show_gaps=true)
            @test isfile(path)
        end
    end

    @testset "subplot layout=(1,2)" begin
        mktempdir() do tmp
            path = joinpath(tmp, "subplot_1x2.pdf")
            savefig_publication([bs2, bs2], path; layout=(1, 2))
            @test isfile(path)
        end
    end

    @testset "subplot layout=(2,1)" begin
        mktempdir() do tmp
            path = joinpath(tmp, "subplot_2x1.pdf")
            savefig_publication([bs2, bs2], path; layout=(2, 1))
            @test isfile(path)
        end
    end

    @testset "overlay モード" begin
        mktempdir() do tmp
            path = joinpath(tmp, "overlay.pdf")
            savefig_publication([bs2, bs2], path; overlay=true)
            @test isfile(path)
        end
    end

    @testset "PNG 出力" begin
        mktempdir() do tmp
            path = joinpath(tmp, "bs2.png")
            savefig_publication(bs2, path)
            @test isfile(path)
        end
    end

    @testset "ylabel kwarg" begin
        mktempdir() do tmp
            path = joinpath(tmp, "ylabel.pdf")
            savefig_publication(bs2, path; ylabel="Frequency [THz]")
            @test isfile(path)
        end
    end

    @testset "L1 plot_on_axis!" begin
        # returns ax, no exception, composes into a user subplot grid
        mktempdir() do tmp
            fig = PythonPlot.figure(; figsize=(8, 4))
            ax1 = fig.add_subplot(1, 2, 1)
            ax2 = fig.add_subplot(1, 2, 2)
            @test PhoXonic.plot_on_axis!(ax1, bs2; title="TE") === ax1
            PhoXonic.plot_on_axis!(ax2, bs2)
            path = joinpath(tmp, "panel.pdf")
            fig.savefig(path)
            @test isfile(path)
            PythonPlot.close(fig)
        end
        # title: non-empty applies, empty leaves no title
        let fig = PythonPlot.figure()
            ax = fig.add_subplot()
            PhoXonic.plot_on_axis!(ax, bs2; title="MyBands")
            @test occursin("MyBands", string(ax.get_title()))
            PythonPlot.close(fig)
        end
        let fig = PythonPlot.figure()
            ax = fig.add_subplot()
            PhoXonic.plot_on_axis!(ax, bs2)
            @test isempty(string(ax.get_title()))
            PythonPlot.close(fig)
        end
        # xlabel absorbed by the wrapper (PhoXonic-specific noop resolved)
        let fig = PythonPlot.figure()
            ax = fig.add_subplot()
            PhoXonic.plot_on_axis!(ax, bs2; xlabel="Wave vector")
            @test occursin("Wave vector", string(ax.get_xlabel()))
            PythonPlot.close(fig)
        end
        # PhoXonic-specific kwargs flow through L1 to the helper
        let fig = PythonPlot.figure()
            ax = fig.add_subplot()
            @test PhoXonic.plot_on_axis!(ax, bs2; show_gaps=true, normalize=2.0) === ax
            PythonPlot.close(fig)
        end
    end

    @testset "L2 figure_publication" begin
        res = PhoXonic.figure_publication(bs2)
        @test res isa Tuple && length(res) == 2
        fig, ax = res
        @test fig isa PythonPlot.Figure
        PythonPlot.close(fig)
        # finalization kwarg (ylims) routed through L2
        fig2, _ = PhoXonic.figure_publication(bs2; ylims=(0.0, 1.0))
        PythonPlot.close(fig2)
    end

    @testset "L3 Vector overlay forwards kwargs" begin
        mktempdir() do tmp
            path = joinpath(tmp, "overlay_kw.pdf")
            # overlay branch now forwards plot kwargs to L1 (previously dropped)
            savefig_publication([bs2, bs2], path; overlay=true, show_gaps=true)
            @test isfile(path)
        end
    end

    @testset "L1/L2 fallback on bad-type args" begin
        # ext loaded, but an unsupported type falls through to the untyped fallback
        @test_throws ArgumentError PhoXonic.plot_on_axis!(nothing, nothing)
        @test_throws ArgumentError PhoXonic.figure_publication(nothing)
    end
end
