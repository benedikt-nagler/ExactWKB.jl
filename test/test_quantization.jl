using Resurgence: FormalSeries, median_sum
import QuadGK

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

    # ── the uniform (Weber) condition through the barrier top ───────────────────
    #
    # V = (z²−1)² has its barrier top at V(0) = 1. Below it the layout is the familiar
    # two wells and a barrier; above it there is no barrier and no doublet, and the
    # default route has nothing to stand on. The cycles do not care, and neither does
    # the Weber connection formula.

    @testset "uniform_cycles continue through the barrier top" begin
        for E in (0.9, 0.999, 1.001, 1.1, 2.0)
            c = uniform_cycles(dwell, E)
            @test c.below == (E < 1)
            pE = with_energy(dwell, E)
            # ledger item 1: −Im v₋₁/(2ħ) is the REAL half-well phase ∫_a^0√(E−V)/ħ on
            # both sides, which is what fixes the +i inner turning point above the top
            vW = period_integral(pE, c.well.contour)
            a = real(c.well.left)
            Φ0 = first(QuadGK.quadgk(x -> sqrt(complex(E - (x^2 - 1)^2)), a, 0.0;
                                     rtol = 1e-12))
            @test -imag(vW) / 2 ≈ real(Φ0) rtol = 1e-8
            # ledger item 3: the barrier index changes sign at the top, continuously
            vB = real(period_integral(pE, c.barrier.contour))
            ε = (c.below ? 1 : -1) * abs(vB) / 2π
            @test (E < 1) == (ε > 0)
        end
        # ε → 0 from both sides, and the two straddling values are close
        εof(E) = begin
            c = uniform_cycles(dwell, E)
            (c.below ? 1 : -1) *
            abs(real(period_integral(with_energy(dwell, E), c.barrier.contour))) / 2π
        end
        @test εof(0.999) > 0 > εof(1.001)
        @test εof(0.999) + εof(1.001) ≈ 0 atol = 1e-4
        # ledger item 2, the other half: the barrier ellipse must not swallow the OUTER
        # turning points, which it does at low energy if the padding is measured to the
        # pair's midpoint rather than to its segment. Deep down, its classical period
        # is 2S_I and nothing else.
        for E in (0.1, 0.19, 0.45)
            c = uniform_cycles(dwell, E)
            b = sqrt(1 - sqrt(E))                   # inner turning points ±b
            SI = first(QuadGK.quadgk(x -> sqrt((x^2 - 1)^2 - E), -b, b; rtol = 1e-12))
            @test abs(real(period_integral(with_energy(dwell, E), c.barrier.contour))) ≈
                  2 * SI rtol = 1e-8
        end
        # refusals: not a four-turning-point double well, and the merged pair itself
        @test_throws QuantizationError uniform_cycles(harm, 1.0)
        @test_throws CoalescentTurningPoints uniform_cycles(dwell, 1.0)
    end

    @testset "uniform quantization above the barrier top (ledger item 5)" begin
        # Dense diagonalization is the oracle. Levels 0-1 are the ground doublet, below
        # the top; levels 2-5 sit ABOVE it, where the default route throws - and where
        # the even/odd members interleave as N = 2n, 2n+1 rather than pairing up.
        ħ = 0.25
        Ed = diagonalization_eigenvalues(dwell, ħ; nev = 6, N = 160, omega = 2.0,
                                         refine = false)
        @test Ed[2] < 1 < Ed[3]                     # the barrier top splits the list
        for (i, (n, p)) in enumerate([(0, 1), (0, -1), (1, 1), (1, -1), (2, 1), (2, -1)])
            E = wkb_eigenvalue(dwell, n, ħ; parity = p, uniform = true, order = 10,
                               rtol = 1e-10)
            @test E ≈ Ed[i] rtol = 3e-4
            @test quantization_condition(dwell, E, ħ; n, parity = p, uniform = true,
                                         order = 10) ≈ 0 atol = 1e-8
        end
        # the four over-barrier levels are exactly the ones the default route cannot
        # reach at all
        for n in 1:2, p in (1, -1)
            @test_throws QuantizationError wkb_eigenvalue(dwell, n, ħ; parity = p,
                                                          order = 8)
        end
    end

    @testset "uniform reproduces the deep condition, up to two instantons" begin
        # Deep below the top the two conditions differ only in the connection constant
        # - Airy ½√(V_A/(1+V_A)) against Weber ½arctan(√V_A) - which is an O(V_A^{3/2})
        # difference. At ħ = 0.25, V_A ≈ 6.6e-3, so the eigenvalues must agree to ~1e-4
        # and no better.
        ħ = 0.25
        for p in (1, -1)
            a = wkb_eigenvalue(dwell, 0, ħ; parity = p, order = 8, rtol = 1e-11)
            b = wkb_eigenvalue(dwell, 0, ħ; parity = p, uniform = true, order = 8,
                               rtol = 1e-11)
            @test a ≈ b rtol = 1e-3
            @test abs(a - b) > 1e-6 * abs(a)        # they are NOT the same condition
        end
    end

    @testset "uniform limits: the over-barrier single well" begin
        # arctan(e^{−πε}) → π/2 as ε → −∞, so the condition becomes 2φ = π(2n+1−p/2):
        # the plain Bohr–Sommerfeld spectrum of the whole well with the even levels at
        # N = 2n and the odd at N = 2n+1. Checked against the classical action of the
        # single allowed region [a, d] at an energy far above the top.
        ħ = 0.25
        for (n, p, N) in ((2, 1, 4), (2, -1, 5))
            E = wkb_eigenvalue(dwell, n, ħ; parity = p, uniform = true, order = 8,
                               rtol = 1e-10)
            a = sqrt(1 + sqrt(E))                   # V = (z²−1)² = E ⇒ z² = 1 ± √E
            J = 2 * first(QuadGK.quadgk(x -> sqrt(max(E - (x^2 - 1)^2, 0.0)), -a, a;
                                        rtol = 1e-10))
            @test J / (2ħ) ≈ π * (N + 0.5) rtol = 2e-2   # leading BS, hence loose
        end
    end

    @testset "the pinched window is refused, and named" begin
        # `φ` is a period of the cycle PINCHED by the vanishing one, so its ħ-series
        # diverges near the top - while `ε`, the vanishing cycle's own index, stays
        # perfectly regular there. That asymmetry is the finding; the refusal says so.
        err = try
            quantization_condition(dwell, 0.999, 0.25; n = 1, parity = 1,
                                   uniform = true, order = 8)
        catch e
            e
        end
        @test err isa QuantizationError
        @test occursin("pinched", sprint(showerror, err))
        # ε at the same energy is finite and order-STABLE (the measured asymmetry -
        # φ moves by nineteen orders of magnitude between order 6 and order 8 there)
        c = uniform_cycles(dwell, 0.999)
        εs = map((2, 4, 6)) do order
            w = wkb_expansion(with_energy(dwell, 0.999); order)
            lVA = ExactWKB._log_voros(voros_symbol(w, c.barrier.contour), 0.25;
                                      theta = 0, side = :median)
            abs(real(lVA)) / 2π
        end
        @test all(isfinite, εs)
        @test maximum(εs) / minimum(εs) - 1 < 0.05
    end

    @testset "layout and argument errors" begin
        # the uniform condition is the parity-factorized one
        @test_throws Resurgence.InvalidArgument quantization_condition(
            dwell, 0.45, 0.25; n = 0, uniform = true)
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
