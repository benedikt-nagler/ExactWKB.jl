using Resurgence: FormalSeries, coefficients, power_offset, borel, pade, poles

@testset "voros" begin
    @testset "harmonic Bohr–Sommerfeld: E = 2ħ(n+½)" begin
        # Our Q = z² − E ⟹ −ħ²ψ″ + z²ψ = Eψ, exact spectrum E_n = 2ħ(n+½).
        # The all-orders WKB truncates: quantum corrections vanish and the
        # classical period alone quantizes as |v₋₁|/(2πħ) = n+½.
        ħ = 0.1
        for n in 0:3
            E = 2 * ħ * (n + 0.5)
            harm = SchrodingerProblem([0.0, 0.0, 1.0]; energy = E)
            w = wkb_expansion(harm; order = 4)
            vs = voros_symbol(w, encircling_contour(turning_points(harm)...))
            @test classical_period(vs) ≈ -im * π * E rtol = 1e-8
            # quantum corrections vanish (Float64-honest bound; BigFloat reaches ~1e-70)
            @test all(abs.(coefficients(quantum_series(vs))) .< 1e-6)
            @test abs(classical_period(vs)) / (2π * ħ) ≈ n + 0.5 rtol = 1e-8
        end
        # BigFloat reaches the theorem's true smallness: |v_m| < 1e-30
        setprecision(256) do
            E = big"1.0"
            harm = SchrodingerProblem([big(0.0), big(0.0), big(1.0)]; energy = E)
            w = wkb_expansion(harm; order = 4)
            vs = voros_symbol(w, encircling_contour(turning_points(harm)...))
            @test all(abs.(coefficients(quantum_series(vs))) .< 1e-30)
        end
    end

    @testset "series shape & full_series" begin
        harm = SchrodingerProblem([0.0, 0.0, 1.0]; energy = 1.0)
        w = wkb_expansion(harm; order = 4)
        vs = voros_symbol(w, encircling_contour(turning_points(harm)...))
        s = quantum_series(vs)
        @test power_offset(s) == 1 // 1
        @test all(iszero, coefficients(s)[2:2:end])          # even coeffs exact zero
        fs = full_series(vs)
        @test power_offset(fs) == -1 // 1
        c = coefficients(fs)
        @test c[1] == classical_period(vs)                   # ħ^{-1}
        @test iszero(c[2])                                   # ħ^0 (Maslov slot)
        @test c[3] == coefficients(s)[1]                     # ħ^1 = v₁
    end

    @testset "double-well handoff to Resurgence" begin
        # Inner pair (barrier) of the double well: Q > 0 between them, so the
        # classical period (instanton action) is real. Handoff smoke test.
        dw = SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0])
        tps = sort(turning_points(dw); by = t -> real(location(t)))
        w = wkb_expansion(dw; order = 3)
        vs = voros_symbol(w, encircling_contour(tps[2], tps[3]; margin = 0.25))
        @test abs(imag(classical_period(vs))) < 1e-8         # real Z (instanton)
        @test real(classical_period(vs)) < 0
        B = borel(quantum_series(vs))                        # borel ∘ series runs
        # equal-degree Padé is degenerate on an alternating-zero series (documented);
        # reduction handles it and poles run.
        r = pade(B; reduce = true)
        @test poles(r) isa AbstractVector
    end

    @testset "the derivative Voros symbol" begin
        # ∂_λ commutes with Borel summation, so the derivative series can be pushed
        # through the ordinary summation path: `_log_voros` of the derivative symbol is
        # ∂_λ log V. Oracle: a finite difference of the summed log at frozen contour.
        prob = SchrodingerProblem([0 // 1, 0 // 1, -4 // 1, 0 // 1, 1 // 1];
                                  energy = -2 // 1)
        E0, ħ = -2 // 1, 0.15
        c = spectral_cycles(prob, E0)[1].contour
        w = wkb_expansion(prob; order = 8)
        dvs = voros_symbol(wkb_derivative(w), c)

        # the classical slot is the period derivative, and the even slots still vanish
        @test classical_period(dvs) ≈ period_derivative(prob, c) rtol = 1e-20
        @test all(iszero, coefficients(quantum_series(dvs))[2:2:end])
        @test any(!iszero, coefficients(quantum_series(dvs))[1:2:end])

        h = 1e-6
        lp = ExactWKB._log_voros(
            voros_symbol(wkb_expansion(with_energy(prob, E0 + h); order = 8), c), ħ)
        lm = ExactWKB._log_voros(
            voros_symbol(wkb_expansion(with_energy(prob, E0 - h); order = 8), c), ħ)
        @test ExactWKB._log_voros(dvs, ħ) ≈ (lp - lm) / (2h) rtol = 1e-6
    end

    @testset "argument checks" begin
        harm = SchrodingerProblem([0.0, 0.0, 1.0]; energy = 1.0)
        w0 = wkb_expansion(harm; order = 0)
        @test_throws Resurgence.InvalidArgument voros_symbol(
            w0, encircling_contour(turning_points(harm)...))
        @test_throws Resurgence.InvalidArgument voros_symbol(
            wkb_derivative(w0), encircling_contour(turning_points(harm)...))
    end
end
