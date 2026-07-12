using ExactWKB
using Test

@testset "Makie extension" begin
    # before a backend is loaded, the stub gives a helpful hint
    dw = SchrodingerProblem([3//4, 0, -2, 0, 1])
    g = stokes_graph(dw; theta = 0.0)

    # load the headless backend, which activates ExactWKBMakieExt
    using CairoMakie

    @testset "single graph renders to PNG" begin
        fig = plot_stokes_graph(g)
        path = tempname() * ".png"
        CairoMakie.save(path, fig)
        @test isfile(path)
        @test filesize(path) > 0
        rm(path; force = true)
    end

    @testset "θ-family renders" begin
        fam = stokes_graph_family(dw, [-0.05, 0.0, 0.05])
        fig = plot_stokes_graph(fam)
        path = tempname() * ".png"
        CairoMakie.save(path, fig)
        @test isfile(path)
        @test filesize(path) > 0
        rm(path; force = true)
    end
end
