using Resurgence: FormalSeries, coefficients, power_offset, variable, n_terms

# Substitute the truncated S = Σ_{m=-1}^{order} S_m(z0) ħ^m back into the Riccati
# equation S² + S′ = ħ⁻² Q at a numeric point; the residual must vanish for every ħ
# power that the truncation determines (ħ^{-2} … ħ^{order-1}).
function _riccati_residual(w, z0)
    z = Complex{BigFloat}(z0)
    Qz = w.prob(z); u = sqrt(Qz)
    N = ExactWKB.order(w) + 2
    Sc = [ExactWKB._eval_term(ExactWKB._s_term(w, m), z, Qz, u) for m in -1:ExactWKB.order(w)]
    Spc = [ExactWKB._eval_term(ExactWKB._deriv_term(ExactWKB._s_term(w, m), w.Q, w.Qp, w.br),
                               z, Qz, u) for m in -1:ExactWKB.order(w)]
    S = FormalSeries(Sc, :ħ; power_offset = -1 // 1)
    Sp = FormalSeries(Spc, :ħ; power_offset = -1 // 1)
    resid = S * S + Sp                        # offset -2; ħ^{-2} coeff should be Q
    rc = coefficients(resid)
    rc[1] -= Qz                               # subtract Q(z0) from the ħ^{-2} term
    # check the safely-determined prefix (powers ħ^{-2} … ħ^{order-2})
    maximum(abs, rc[1:(ExactWKB.order(w))])
end

# The same substitution for the differentiated Riccati equation: with ∂Q = D,
# S² + S′ = ħ⁻²Q differentiates to 2S·∂S + (∂S)′ = ħ⁻²·D. Residual over the same
# safely-determined prefix.
function _driccati_residual(w, dw, z0)
    z = Complex{BigFloat}(z0)
    Qz = w.prob(z); u = sqrt(Qz)
    ev(t) = ExactWKB._eval_term(t, z, Qz, u)
    ms = -1:ExactWKB.order(w)
    S = FormalSeries([ev(ExactWKB._s_term(w, m)) for m in ms], :ħ; power_offset = -1 // 1)
    dS = FormalSeries([ev(ExactWKB._s_term(dw, m)) for m in ms], :ħ; power_offset = -1 // 1)
    dSp = FormalSeries(
        [ev(ExactWKB._deriv_term(ExactWKB._s_term(dw, m), w.Q, w.Qp, w.br)) for m in ms],
        :ħ; power_offset = -1 // 1)
    resid = 2 * (S * dS) + dSp                # offset -2; ħ^{-2} coeff should be D(z0)
    rc = coefficients(resid)
    rc[1] -= ExactWKB._poly_eval(dw.D, z)
    maximum(abs, rc[1:(ExactWKB.order(w))])
end

@testset "wkb_recursion" begin
    airy = SchrodingerProblem([0, 1])   # Q = z

    @testset "hand-derived Airy S_m (exact ring, ==)" begin
        w = wkb_expansion(airy; order = 4, arithmetic = :exact)
        # S_m = num·√z^ε / z^k
        expected = Dict(-1 => ([1 // 1], 0, 1),
                        0 => ([-1 // 4], 1, 0),
                        1 => ([-5 // 32], 3, 1),
                        2 => ([-15 // 64], 4, 0),
                        3 => ([-1105 // 2048], 6, 1),
                        4 => ([-1695 // 1024], 7, 0))
        for (m, (num, qp, sp)) in expected
            t = ExactWKB._s_term(w, m)
            @test Rational{BigInt}.(ExactWKB._numerator_coeffs(t)) == Rational{BigInt}.(num)
            @test t.q_pow == qp
            @test t.sqrt_pow == sp
            @test t.sqrt_pow == mod(m, 2)          # parity of sqrt_pow
        end
    end

    @testset "Riccati residual vanishes to order" begin
        w = wkb_expansion(airy; order = 8, arithmetic = :exact)
        @test _riccati_residual(w, 3 // 7 + 2im) < 1e-60
        # a nontrivial potential
        w2 = wkb_expansion(SchrodingerProblem([3 // 4, 0, -2, 0, 1]); order = 6)
        @test _riccati_residual(w2, 1 // 3 + 1im) < 1e-40
    end

    @testset "even/odd reduction oracle" begin
        w = wkb_expansion(airy; order = 8, arithmetic = :exact)
        @test even_odd_residual(w) < 1e-60
        w2 = wkb_expansion(SchrodingerProblem([0, 0, 1]; energy = 1); order = 8)
        @test even_odd_residual(w2) < 1e-30
    end

    @testset "exact vs BigFloat ring agreement" begin
        setprecision(256) do
            we = wkb_expansion(airy; order = 6, arithmetic = :exact)
            wb = wkb_expansion(airy; order = 6, arithmetic = :bigfloat)
            z = Complex{BigFloat}(5 // 4, 1 // 3)
            for m in -1:6
                Qz = airy(z); u = sqrt(Qz)
                a = ExactWKB._eval_term(ExactWKB._s_term(we, m), z, Qz, u)
                b = ExactWKB._eval_term(ExactWKB._s_term(wb, m), z, Qz, u)
                @test abs(a - b) < 1e-40
            end
        end
    end

    @testset "evaluate_s_odd → FormalSeries" begin
        w = wkb_expansion(airy; order = 6, arithmetic = :exact)
        S = evaluate_s_odd(w, 2.0 + 0.5im)
        @test variable(S) === :ħ
        @test power_offset(S) == -1 // 1
        c = coefficients(S)
        # even entries (m = 0, 2, 4, 6 → indices 2,4,6,8) are exactly zero
        @test all(iszero, c[2:2:end])
        @test !iszero(c[1])                              # S_{-1} = √Q ≠ 0
        # branch flips the sign of every odd term
        Sm = evaluate_s_odd(w, 2.0 + 0.5im; branch = -1)
        @test coefficients(Sm)[1] ≈ -c[1]
    end

    @testset "argument checks" begin
        @test_throws Resurgence.InvalidArgument wkb_expansion(airy; order = -1)
        @test_throws Resurgence.InvalidArgument wkb_expansion(
            SchrodingerProblem([0.0, 1.0]); arithmetic = :exact)
        @test_throws Resurgence.InvalidArgument wkb_expansion(airy; arithmetic = :nonsense)
    end

    @testset "parameter derivatives" begin
        # The oracle for the whole derivative tower at once: differentiate the Riccati
        # equation itself. S² + S′ = ħ⁻²Q  ⟹  2S·∂S + (∂S)′ = ħ⁻²·D. Nothing about the
        # numerics enters - if the differentiated recursion is right, the residual is
        # zero at every ħ power the truncation determines.
        quartic = SchrodingerProblem([0 // 1, 0 // 1, -4 // 1, 0 // 1, 1 // 1];
                                     energy = -2 // 1)
        for (prob, wrt, idx) in ((airy, :energy, 0), (airy, :coefficient, 1),
                                 (quartic, :energy, 0), (quartic, :coefficient, 0),
                                 (quartic, :coefficient, 3))
            w = wkb_expansion(prob; order = 6)
            dw = wkb_derivative(w; wrt, index = idx)
            @test _driccati_residual(w, dw, 0.4 + 1.1im) < 1e-20
        end

        @testset "∂/∂E = −∂/∂v₀ (Q = V − E)" begin
            w = wkb_expansion(quartic; order = 6)
            dE = wkb_derivative(w; wrt = :energy)
            d0 = wkb_derivative(w; wrt = :coefficient, index = 0)
            for m in -1:6
                a, b = ExactWKB._s_term(dE, m), ExactWKB._s_term(d0, m)
                @test a.q_pow == b.q_pow && a.sqrt_pow == b.sqrt_pow
                @test ExactWKB._numerator_coeffs(a) == -ExactWKB._numerator_coeffs(b)
            end
        end

        @testset "shape invariants" begin
            w = wkb_expansion(quartic; order = 6)
            dw = wkb_derivative(w; wrt = :energy)
            for m in -1:6
                s, ds = ExactWKB._s_term(w, m), ExactWKB._s_term(dw, m)
                @test ds.sqrt_pow == s.sqrt_pow      # ∂ preserves the √Q parity
                @test ds.sqrt_pow == mod(m, 2)
            end
            @test ExactWKB.order(dw) == ExactWKB.order(w)
        end

        @testset "argument checks" begin
            w = wkb_expansion(quartic; order = 2)
            @test_throws Resurgence.InvalidArgument wkb_derivative(w; wrt = :nonsense)
            @test_throws Resurgence.InvalidArgument wkb_derivative(
                w; wrt = :coefficient, index = -1)
            @test_throws Resurgence.InvalidArgument wkb_derivative(
                w; wrt = :coefficient, index = 5)   # V has degree 4
        end
    end
end
