using Resurgence: FormalSeries, median_sum

@testset "quantization" begin
    harm = SchrodingerProblem([0.0, 0.0, 1.0])            # V = z², E_n = 2ħ(n+½)
    quart = SchrodingerProblem([0, 0, 1, 0, 1])           # V = z² + z⁴ (exact coeffs)
    dwell = SchrodingerProblem([1.0, 0.0, -2.0, 0.0, 1.0]) # V = (z²−1)²

    @testset "harmonic exact at every order (convention pin)" begin
        # the WKB series truncates: the condition is exact, not asymptotic
        for ħ in (0.1, 0.35), n in 0:3
            E = wkb_eigenvalue(harm, n, ħ; order = 6)
            @test E ≈ 2 * ħ * (n + 0.5) rtol = 1e-8
            @test quantization_condition(harm, E, ħ; n, order = 6) ≈ 0 atol = 1e-6
        end
        # BigFloat: full-precision exactness of the pinned condition
        setprecision(160) do
            harmx = SchrodingerProblem([0, 0, 1])          # exact coeffs → BigFloat
            ħ = big"0.2"
            E = wkb_eigenvalue(harmx, 1, ħ; order = 6, rtol = 1e-30,
                               quad_rtol = 1e-35)
            @test E ≈ 2 * ħ * (1 + big"0.5") rtol = 1e-25
        end
    end

    @testset "quartic ground state: WKB vs diagonalization vs Bender–Wu" begin
        # Digit target: the single-well condition is exact up to the quartic's
        # COMPLEX-turning-point (hidden instanton) sectors, O(e^{−S_c/ħ}) - measured
        # residual floor at the true eigenvalue: −1.0e-4 (ħ = 0.1), −9.4e-8
        # (ħ = 0.05), order-converged, so it is physics, not summation error.
        # Complex saddles are Tier B; at ħ = 0.05 the honest target is ~7 digits.
        ħ = 0.05
        E_wkb = setprecision(160) do
            wkb_eigenvalue(quart, 0, big(ħ); order = 14, rtol = 1e-10,
                           quad_rtol = 1e-18)
        end
        E_diag = setprecision(160) do
            diagonalization_eigenvalues(quart, ħ; nev = 1, N = 80)[1]
        end
        # (a) independent diagonalization oracle
        @test E_wkb ≈ E_diag rtol = 1e-6
        @test abs(Float64(E_wkb - E_diag)) > 0    # the floor is real: not exact
        # (b) the Resurgence :quartic Bender–Wu series: −ħ²ψ″ + (y²+y⁴)ψ = Ẽψ maps
        # to H = p²/2 + x²/2 + gx⁴ by y = √ħ·x, giving Ẽ(ħ) = 2ħ·E₀(g = ħ/2). The
        # true Borel singularity sits at −1/3, but the Padé sprinkles a tiny-residue
        # spurious pole on the positive axis - the median sum integrates around it.
        bw = Resurgence.coefficients(FormalSeries(:quartic, 40))
        tail = FormalSeries(Float64.(bw[2:end]), :g; power_offset = 1 // 1)
        E0 = Float64(bw[1]) + real(median_sum(tail, ħ / 2))
        @test Float64(E_wkb) ≈ 2 * ħ * E0 rtol = 1e-6
        # excited state against the oracle (Float64 pipeline, ħ = 0.1: the complex-
        # instanton floor is ~1e-5 relative there)
        E2 = wkb_eigenvalue(SchrodingerProblem([0.0, 0.0, 1.0, 0.0, 1.0]), 2, 0.1;
                            order = 12)
        E2_diag = setprecision(96) do
            diagonalization_eigenvalues(quart, 0.1; nev = 3, N = 80)[3]
        end
        @test E2 ≈ Float64(E2_diag) rtol = 1e-4
    end

    @testset "spectral determinant zero ↔ eigenvalue" begin
        ħ = 0.25
        n = 1
        E = wkb_eigenvalue(harm, n, ħ; order = 6)
        D0 = spectral_determinant(harm, E, ħ; order = 6)
        δ = 1e-2
        @test abs(D0) < 1e-4 * abs(spectral_determinant(harm, E + δ, ħ; order = 6))
        # D(E) = 1 + V_B ≈ i·(J′/ħ)(E − E*) near the zero: the sign change is in
        # the imaginary part (the real part is quadratic and does not cross)
        @test sign(imag(spectral_determinant(harm, E - δ, ħ; order = 6))) !=
              sign(imag(spectral_determinant(harm, E + δ, ħ; order = 6)))
    end

    @testset "layout and argument errors" begin
        # double well needs parity; single well refuses parity
        @test_throws QuantizationError quantization_condition(dwell, 0.25, 0.05; n = 0)
        @test_throws QuantizationError quantization_condition(harm, 1.0, 0.1; n = 0,
                                                              parity = 1)
        @test_throws Resurgence.InvalidArgument quantization_condition(
            dwell, 0.25, 0.05; n = 0, parity = 2)
        @test_throws Resurgence.InvalidArgument wkb_eigenvalue(harm, -1, 0.1)
        # E below the potential floor has no allowed region
        @test_throws QuantizationError quantization_condition(harm, -1.0, 0.1; n = 0)
    end
end
