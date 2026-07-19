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

    @testset "M4 bridge types" begin
        cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
        t = ideal_triangulation(stokes_graph(cubic; theta = 0.3))
        @test occursin("5-gon", compact(t))
        @test occursin("2 diagonals", compact(t))
        @test occursin("γ(1, 2)", plain(t))
        cb = charge_basis(cubic, t)
        @test occursin("ChargeBasis(2 cycles)", compact(cb))
        @test occursin("signed pairing", plain(cb))
        # the signed frame is displayed: ε per cycle, and the decay representative
        # spelled out wherever it differs from the physical charge (here ε = [1, -1])
        @test occursin("ε = ", plain(cb))
        @test occursin("decay rep.", plain(cb))
        sp = bps_spectrum(cubic; theta = 0.3)
        @test occursin("BPSSpectrum(2 states", compact(sp))
        @test occursin("Ω = 1", plain(sp))
        @test occursin("BPSState", compact(sp.states[1]))
    end
end
