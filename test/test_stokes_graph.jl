using ExactWKB
using ExactWKB: lines, points
using Test

@testset "stokes_graph" begin
    @testset "Airy: three rays to infinity at 2π/3" begin
        airy = SchrodingerProblem([0, 1])          # Q = z, one simple TP at 0
        g = stokes_graph(airy; theta = 0.0)
        @test length(turning_points(g)) == 1
        @test n_infinite_lines(g) == 3
        @test isempty(edges(g))
        @test isempty(finite_lines(g))
        @test topology_signature(g) == (1, 3, Tuple{Int,Int}[])

        # exit angles are 0, 2π/3, 4π/3 (mod 2π), one per ray
        exits = sort([mod2pi(angle(points(l)[end])) for l in lines(g)])
        @test isapprox(exits, [0.0, 2π / 3, 4π / 3]; atol = 1e-4)
        # every ray escapes with the same mass (symmetry)
        @test all(l -> l.endpoint === :infinity, lines(g))
    end

    @testset "θ-rotation moves rays by the seeding formula" begin
        airy = SchrodingerProblem([0, 1])
        for θ in (0.0, 0.3, 0.6, -0.4)
            g = stokes_graph(airy; theta = θ)
            # ray k has exit direction (2θ + 2πk)/3; k = 0 → 2θ/3
            a0 = angle(points(lines(g)[1])[end])
            @test isapprox(a0, 2θ / 3; atol = 1e-4)
        end
    end

    @testset "w-drift is killed by re-projection" begin
        # Airy θ = 0, ray k = 0 is exactly the +real axis; any imaginary part is drift.
        g = stokes_graph(SchrodingerProblem([0, 1]); theta = 0.0)
        l0 = lines(g)[1]
        drift = maximum(abs(imag(z)) for z in points(l0))
        @test drift < 1e-8
    end

    @testset "order-≥2 turning points are refused" begin
        @test_throws UnsupportedTurningPoint stokes_graph(SchrodingerProblem([0, 0, 1]))  # Q = z²
    end

    @testset "double well: inner saddle + stable signature" begin
        dw = SchrodingerProblem([3//4, 0, -2, 0, 1])   # (z²−1)² − 1/4
        g = stokes_graph(dw; theta = 0.0)
        @test length(turning_points(g)) == 4

        # canonical turning-point order is by position: inner pair are indices 2,3
        locs = [real(location(t)) for t in turning_points(g)]
        @test issorted(locs)

        # exactly one saddle: the inner pair, double-traced
        es = edges(g)
        @test length(es) == 1
        @test es[1] == (2, 3)
        @test count(is_finite_line, lines(g)) == 2      # traced from both endpoints
        @test n_infinite_lines(g) == 10                 # 4·3 − 2

        # the two half-saddle line-masses agree (½|Z| each direction)
        halfmasses = [l.mass for l in finite_lines(g)]
        @test isapprox(halfmasses[1], halfmasses[2]; rtol = 1e-4)

        # topology signature is stable under seed / tolerance halving
        s1 = topology_signature(g)
        s2 = topology_signature(stokes_graph(dw; seed_radius = 0.005, reltol = 5e-10))
        @test s1 == s2
        @test s1 == (4, 10, [(2, 3)])
    end

    @testset "off-critical phase breaks the saddle" begin
        # away from θ_c = 0 the inner Stokes line no longer connects the two turning
        # points, so no finite edge survives.
        dw = SchrodingerProblem([3//4, 0, -2, 0, 1])
        g = stokes_graph(dw; theta = 0.25)
        @test isempty(edges(g))
        @test n_infinite_lines(g) == 12
    end

    @testset "BigFloat opt-in traces at high precision" begin
        setprecision(BigFloat, 128) do
            g = stokes_graph(SchrodingerProblem([0, 1]); theta = 0.0, bigfloat = true)
            @test eltype(points(lines(g)[1])) == Complex{BigFloat}
            @test n_infinite_lines(g) == 3
        end
    end
end
