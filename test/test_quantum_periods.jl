using Resurgence: coefficients
using QuadGK: quadgk

@testset "quantum_periods" begin
    # family fixtures: potentials with energy 0, energies passed per call
    harm = SchrodingerProblem([0.0, 0.0, 1.0])            # V = z²
    quart = SchrodingerProblem([0.0, 0.0, 1.0, 0.0, 1.0]) # V = z² + z⁴
    dwell = SchrodingerProblem([1.0, 0.0, -2.0, 0.0, 1.0]) # V = (z² − 1)²

    @testset "cycle identification" begin
        # harmonic at E = 1: one well cycle between ±1
        cs = spectral_cycles(harm, 1.0)
        @test length(cs) == 1
        @test kind(cs[1]) === :well
        l, r = endpoints(cs[1])
        @test l ≈ -1 && r ≈ 1
        # quartic at E = 1: one well; the two imaginary turning points are ignored
        cs = spectral_cycles(quart, 1.0)
        @test length(cs) == 1 && kind(cs[1]) === :well
        # double well at E = 1/4: well, barrier, well in ascending position
        cs = spectral_cycles(dwell, 0.25)
        @test [kind(c) for c in cs] == [:well, :barrier, :well]
        @test endpoints(cs[1])[1] ≈ -sqrt(1.5) && endpoints(cs[3])[2] ≈ sqrt(1.5)
        @test endpoints(cs[2]) == (endpoints(cs[1])[2], endpoints(cs[3])[1])
        # cubic Q = z³ − z at E = 0: the forbidden interval (−1, 0) borders the
        # unbounded allowed region, not two wells - only the (0, 1) well survives
        cs = spectral_cycles(SchrodingerProblem([0.0, -1.0, 0.0, 1.0]), 0.0)
        @test length(cs) == 1 && kind(cs[1]) === :well
        @test endpoints(cs[1])[1] ≈ 0 atol = 1e-10
        # E below the potential minimum: no cycles
        @test isempty(spectral_cycles(harm, -1.0))
    end

    @testset "coalescence refusal" begin
        # E at the double-well bottom (double turning points at ±1)
        @test_throws CoalescentTurningPoints spectral_cycles(dwell, 0.0)
        # E at the barrier top (double turning point at 0)
        @test_throws CoalescentTurningPoints spectral_cycles(dwell, 1.0)
        # E just below the barrier top: two real turning points collide near 0
        @test_throws CoalescentTurningPoints spectral_cycles(dwell, 1.0 - 1e-12)
        # non-real potential is out of the M5 scope
        @test_throws Resurgence.InvalidArgument spectral_cycles(
            SchrodingerProblem([1.0im, 0.0, 1.0]), 1.0)
    end

    @testset "harmonic quantum period truncates (orientation pin)" begin
        # v₋₁ = −iπE and every quantum correction integrates to zero around the
        # cycle - the all-orders Bohr–Sommerfeld condition is exact (M5 ledger).
        for E in (0.5, 1.0, 2.5)
            vs = quantum_period(harm, :well, E; order = 6)
            @test classical_period(vs) ≈ -im * π * E rtol = 1e-8
            @test all(abs.(coefficients(quantum_series(vs))) .< 1e-6)
        end
    end

    @testset "quartic classical period vs quadrature" begin
        E = 1.0
        cs = spectral_cycles(quart, E)
        vs = quantum_period(quart, cs[1], E; order = 2)
        # v₋₁ = −i J(E), J = 2∫√(E − V) over the allowed interval
        a = real(endpoints(cs[1])[2])
        J, _ = quadgk(x -> 2 * sqrt(max(E - x^2 - x^4, 0.0)), -a, a; rtol = 1e-10)
        @test classical_period(vs) ≈ -im * J rtol = 1e-7
        # well orientation: v₋₁ on the negative imaginary axis
        @test real(classical_period(vs)) ≈ 0 atol = 1e-8 * J
        @test imag(classical_period(vs)) < 0
    end

    @testset "double-well barrier cycle: real negative Z" begin
        cs = spectral_cycles(dwell, 0.25)
        vs = quantum_period(dwell, :barrier, 0.25; order = 3)
        @test abs(imag(classical_period(vs))) < 1e-8
        @test real(classical_period(vs)) < 0
        # the two well periods agree by symmetry
        v1 = quantum_period(dwell, cs[1], 0.25; order = 3)
        v3 = quantum_period(dwell, cs[3], 0.25; order = 3)
        @test classical_period(v1) ≈ classical_period(v3) rtol = 1e-8
        @test coefficients(quantum_series(v1))[1] ≈
              coefficients(quantum_series(v3))[1] rtol = 1e-6
    end

    @testset "kind-form ambiguity" begin
        @test_throws Resurgence.InvalidArgument quantum_period(dwell, :well, 0.25)
        @test_throws Resurgence.InvalidArgument quantum_period(harm, :barrier, 1.0)
    end
end
