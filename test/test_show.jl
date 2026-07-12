using ExactWKB
using Test

# helpers to capture both show forms as strings
compact(x) = sprint(show, x)
plain(x) = sprint((io, v) -> show(io, MIME"text/plain"(), v), x)

@testset "show" begin
    dw = SchrodingerProblem([3//4, 0, -2, 0, 1])

    @testset "SchrodingerProblem" begin
        @test occursin("SchrodingerProblem", compact(dw))
        @test occursin("deg 4", compact(dw))
        p = plain(dw)
        @test occursin("Q(z)", p)
        @test occursin("E =", p)
    end

    @testset "TurningPoint" begin
        tp = first(turning_points(dw))
        @test occursin("TurningPoint", compact(tp))
        @test occursin("simple", compact(tp))
    end

    @testset "WKBExpansion" begin
        w = wkb_expansion(dw; order = 2)
        @test occursin("WKBExpansion", compact(w))
        @test occursin("order 2", compact(w))
        @test occursin("Riccati", plain(w))
    end

    @testset "VorosSymbol" begin
        ho = SchrodingerProblem([0, 0, 1]; energy = 1)   # harmonic Q = z² − 1
        tps = turning_points(ho)
        vs = voros_symbol(wkb_expansion(ho; order = 2), encircling_contour(tps[1], tps[2]))
        @test occursin("VorosSymbol", compact(vs))
        @test occursin("v₋₁", compact(vs))
    end

    @testset "StokesGraph and StokesLine" begin
        g = stokes_graph(dw; theta = 0.0)
        @test occursin("StokesGraph", compact(g))
        @test occursin("saddles", compact(g))
        @test occursin("signature", plain(g))
        l = first(ExactWKB.lines(g))
        @test occursin("StokesLine", compact(l))
        @test occursin("mass", compact(l))
    end

    @testset "Saddle" begin
        s = first(saddles(dw))
        @test occursin("Saddle", compact(s))
        @test occursin("|Z|", compact(s))
        @test occursin("central charge", plain(s))
    end
end
