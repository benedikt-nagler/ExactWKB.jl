# Tests for the pure-SU(2) Seiberg–Witten period module (`src/sw_curve.jl`).
# Oracles: the Mathieu characteristic values a_ν(q) (independent tridiagonal solver,
# checked against DLMF §28.6 small-q series and the a₀(1) constant), the weak-coupling
# u(a) identity, the leading NS quantum correction a₂ ≈ −Λ⁴/(4u^{5/2}), integer
# monodromy, and central-charge additivity.

using LinearAlgebra: SymTridiagonal, eigvals, det

# ── test-side Mathieu characteristic-value oracle (Hill tridiagonal eigenvalues) ──────
# Mathieu equation y″ + (a − 2q cos 2x) y = 0. The characteristic values a_r(q)/b_r(q)
# are eigenvalues of the Fourier-mode tridiagonal matrices (DLMF §28.4).

# even-index a_{2m}(q): diagonal (2k)² = 0,4,16,…, first off-diagonal q√2, rest q
function _mathieu_a_even(q, nmax)
    d = Float64[(2k)^2 for k in 0:nmax]
    e = Float64[k == 1 ? q * sqrt(2) : q for k in 1:nmax]
    eigvals(SymTridiagonal(d, e))
end

# odd-index a_{2m+1}(q): diagonal (2k+1)² = 1,9,25,…, with d₁ += q, off-diagonal q
function _mathieu_a_odd(q, nmax)
    d = Float64[(2k + 1)^2 for k in 0:nmax]; d[1] += q
    e = fill(float(q), nmax)
    eigvals(SymTridiagonal(d, e))
end

# odd-index b_{2m+1}(q): as a_odd but d₁ −= q
function _mathieu_b_odd(q, nmax)
    d = Float64[(2k + 1)^2 for k in 0:nmax]; d[1] -= q
    e = fill(float(q), nmax)
    eigvals(SymTridiagonal(d, e))
end

# even-index b_{2m}(q): diagonal (2k)² = 4,16,36,…, off-diagonal q
function _mathieu_b_even(q, nmax)
    d = Float64[(2k)^2 for k in 1:nmax]
    e = fill(float(q), nmax - 1)
    eigvals(SymTridiagonal(d, e))
end

mathieu_a(q, r; nmax = 24) = iseven(r) ? _mathieu_a_even(q, nmax)[r ÷ 2 + 1] :
                                         _mathieu_a_odd(q, nmax)[(r - 1) ÷ 2 + 1]
mathieu_b(q, r; nmax = 24) = iseven(r) ? _mathieu_b_even(q, nmax)[r ÷ 2] :
                                         _mathieu_b_odd(q, nmax)[(r - 1) ÷ 2 + 1]

@testset "sw_curve" begin
    @testset "Mathieu characteristic-value oracle (DLMF §28.6)" begin
        # small-q series expansions with exact rational coefficients
        for q in (0.05, 0.1, 0.2)
            a0 = -q^2 / 2 + 7q^4 / 128
            @test isapprox(mathieu_a(q, 0), a0; atol = 5e-5)
            a1 = 1 + q - q^2 / 8 - q^3 / 64
            @test isapprox(mathieu_a(q, 1), a1; atol = 5e-5)
            b1 = 1 - q - q^2 / 8 + q^3 / 64
            @test isapprox(mathieu_b(q, 1), b1; atol = 5e-5)
            a2 = 4 + 5q^2 / 12
            @test isapprox(mathieu_a(q, 2), a2; atol = 5e-4)
            b2 = 4 - q^2 / 12
            @test isapprox(mathieu_b(q, 2), b2; atol = 5e-4)
        end
        # the a₀(1) constant (DLMF §28 tabulated value)
        @test isapprox(mathieu_a(1.0, 0), -0.4551386; atol = 1e-6)
    end

    @testset "classical periods & weak-coupling u(a) identity" begin
        sw = SeibergWittenSU2()                     # Λ = 1
        @test dynamical_scale(sw) == 1
        for u in (25.0, 50.0, 100.0)
            p = sw_periods(sw, u)
            a = real(p.a)
            @test imag(p.a) ≈ 0 atol = 1e-9         # electric period is real
            # u(a) = a² + Λ⁴/(2a²) + 5Λ⁸/(32a⁶) + 9Λ¹²/(64a¹⁰) + …
            u_from_a = a^2 + 1 / (2a^2) + 5 / (32a^6) + 9 / (64a^10)
            @test isapprox(u, u_from_a; rtol = 1e-6)
            @test real(p.a_D) ≈ 0 atol = 1e-8       # dual period is imaginary
        end
        # Λ-scaling: singularities at u = ±2Λ²
        sw2 = SeibergWittenSU2(Λ = 2)
        @test sw_singularities(sw2).monopole == 8
        @test sw_singularities(sw2).dyon == -8
    end

    @testset "monopole massless & a_D growth" begin
        sw = SeibergWittenSU2()
        # a_D → 0 as u → 2Λ² = 2 (monopole becomes massless)
        aD_near = abs(sw_periods(sw, 2.001).a_D)
        aD_far  = abs(sw_periods(sw, 4.0).a_D)
        @test aD_near < aD_far
        @test aD_near < 0.05
        # a_D magnitude grows with u (asymptotic-freedom log)
        @test abs(sw_periods(sw, 100.0).a_D) > abs(sw_periods(sw, 10.0).a_D)
    end

    @testset "NS quantum period a₂ ≈ −Λ⁴/(4 u^{5/2})" begin
        sw = SeibergWittenSU2()
        # classical limit: order-0 quantum period equals the classical period
        p0 = quantum_sw_periods(sw, 25.0; order = 0)
        @test Resurgence.coefficients(p0.a)[1] ≈ complex(sw_periods(sw, 25.0).a)
        # leading large-u correction
        for u in (100.0, 400.0)
            p = quantum_sw_periods(sw, u; order = 1)
            a2 = real(Resurgence.coefficients(p.a)[3])    # coeff of ħ²
            @test a2 < 0
            @test isapprox(a2, -1 / (4 * u^2.5); rtol = 5e-3)
        end
    end

    @testset "quantum periods: operator form & magnetic ħ² correction" begin
        sw = SeibergWittenSU2()
        F = Float64

        # Oracle 1 (fully independent of PF): the operator-form electric correction
        # equals the direct Dunham-S₂ quadrature on the electric cycle.
        for u in (5.0, 12.0, 40.0, 150.0)
            op   = real(ExactWKB._electric_a2(sw, u, F))
            quad = ExactWKB._electric_a2_quad(sw, u, F)
            @test isapprox(op, quad; rtol = 1e-7)
        end

        # Oracle 2: the reduced first-order operator equals the raw third-order Dunham
        # operator  D₃ = (1/3)(u²−4Λ⁴)∂³ + u∂² + (1/4)∂  built from finite differences of
        # the closed-form derivatives - for BOTH cycles (certifies PF + the reduction).
        Λ = 1.0
        function d3(fp, u, h)   # fp(u) = closed-form Π′(u)
            Πp   = fp(u)
            Πpp  = (fp(u + h) - fp(u - h)) / (2h)
            Πppp = (fp(u + h) - 2fp(u) + fp(u - h)) / h^2
            (1 / 3) * (u^2 - 4Λ^4) * Πppp + u * Πpp + Πp / 4
        end
        for u in (7.0, 0.5 + 0.0im, -3.0 + 2.0im, 6.0 + 4.0im)
            ap(v)  = ExactWKB._da_du_closed(sw, v, F)
            aDp(v) = ExactWKB._daD_du_closed(sw, v, F)
            @test isapprox(d3(ap, u, 1e-3), ExactWKB._electric_a2(sw, u, F); rtol = 1e-4)
            @test isapprox(d3(aDp, u, 1e-3), ExactWKB._magnetic_aD2(sw, u, F); rtol = 1e-4)
        end

        # magnetic quantum period is now returned (3-term ħ-series) on the whole u-plane
        for u in (25.0, 3.0 + 0.0im, -2.0 + 3.0im, 5.0im)
            p = quantum_sw_periods(sw, u; order = 1)
            @test length(Resurgence.coefficients(p.a_D)) == 3
            aD2 = Resurgence.coefficients(p.a_D)[3]
            @test isfinite(abs(aD2)) && abs(aD2) > 0
            @test aD2 ≈ complex(ExactWKB._magnetic_aD2(sw, u, Float64))
        end

        # reality/branch consistency: on real u > 2Λ² the magnetic period is imaginary,
        # so its ħ² correction is imaginary too (Re ≈ 0)
        for u in (5.0, 30.0)
            aD2 = ExactWKB._magnetic_aD2(sw, u, F)
            @test isapprox(real(aD2), 0.0; atol = 1e-10)
        end

        # order-0 magnetic period equals the classical one
        p0 = quantum_sw_periods(sw, 9.0 + 2.0im; order = 0)
        @test length(Resurgence.coefficients(p0.a_D)) == 1
        @test Resurgence.coefficients(p0.a_D)[1] ≈ complex(sw_periods(sw, 9.0 + 2.0im).a_D)
    end

    @testset "central charge & monodromy" begin
        sw = SeibergWittenSU2()
        u = 30.0
        p = sw_periods(sw, u)
        @test central_charge(sw, u, (0, 1)) ≈ p.a          # electric Z = a
        @test central_charge(sw, u, (1, 0)) ≈ p.a_D        # magnetic Z = a_D
        @test central_charge(sw, u, (1, 1)) ≈
              central_charge(sw, u, (1, 0)) + central_charge(sw, u, (0, 1))
        # integer Sp(2,ℤ) monodromy, det 1, and M_∞ = M_dyon · M_monopole
        Minf = sw_monodromy(sw, :infinity)
        Mmon = sw_monodromy(sw, :monopole)
        Mdyon = sw_monodromy(sw, :dyon)
        for M in (Minf, Mmon, Mdyon)
            @test eltype(M) <: Integer
            @test det(float.(M)) ≈ 1
        end
        @test Minf == Mdyon * Mmon
    end

    @testset "Carlson elliptic layer: known values & Legendre relation" begin
        K(m) = ExactWKB._ellipK(complex(float(m)))
        E(m) = ExactWKB._ellipE(complex(float(m)))
        @test K(0.0) ≈ π / 2
        @test E(0.0) ≈ π / 2
        @test E(1.0) ≈ 1.0
        # DLMF values at parameter m = 1/2
        @test K(0.5) ≈ 1.8540746773013719 rtol = 1e-12
        @test E(0.5) ≈ 1.3506438810476755 rtol = 1e-12
        # Legendre relation E(m)K(1−m) + E(1−m)K(m) − K(m)K(1−m) = π/2, complex m
        for m in (0.3, 0.9, 0.4 + 0.3im, -0.5 + 2.0im, 1.7 - 0.6im)
            L = E(m) * K(1 - m) + E(1 - m) * K(m) - K(m) * K(1 - m)
            @test L ≈ π / 2 atol = 1e-10
        end
        @test_throws PeriodError K(1.0)                       # K diverges at m = 1
    end

    @testset "closed forms vs quadrature (the pinning oracle)" begin
        sw = SeibergWittenSU2()
        for u in (2.5, 4.0, 10.0, 100.0)
            p = sw_periods(sw, u)                             # closed form
            @test p.a ≈ ExactWKB._electric_a(sw, u, Float64) rtol = 1e-6
            @test p.a_D ≈ ExactWKB._magnetic_aD(sw, u, Float64) rtol = 1e-6
        end
        # closed-form derivatives vs finite differences of the closed forms
        for u0 in (5.0 + 0.0im, 1.0 + 2.0im, -3.0 + 1.5im)
            h = 1e-5
            da_fd = (sw_periods(sw, u0 + h).a - sw_periods(sw, u0 - h).a) / (2h)
            daD_fd = (sw_periods(sw, u0 + h).a_D - sw_periods(sw, u0 - h).a_D) / (2h)
            @test ExactWKB._da_du_closed(sw, u0, Float64) ≈ da_fd rtol = 1e-7
            @test ExactWKB._daD_du_closed(sw, u0, Float64) ≈ daD_fd rtol = 1e-7
        end
    end

    @testset "massless states at the singular moduli" begin
        sw = SeibergWittenSU2()
        # monopole exactly massless at u = 2Λ², where a = 4Λ/π (E(1) = 1)
        p = sw_periods(sw, 2.0)
        @test abs(p.a_D) < 1e-12
        @test p.a ≈ 4 / π rtol = 1e-10
        # a dyon (n_m = ±1, n_e = ∓2 up to sign) goes massless as u → −2Λ²
        for r in (1e-2, 1e-3)
            u = -2.0 + r * cis(π / 4)
            Zmin = min(abs(central_charge(sw, u, (-1, 2))), abs(central_charge(sw, u, (1, 2))))
            @test Zmin < 2 * sqrt(r)                          # Z ~ √(u+2Λ²) vanishing rate
        end
    end

    @testset "Picard–Fuchs continuation reproduces the closed forms" begin
        sw = SeibergWittenSU2()
        # single-point path = closed-form data
        base = continue_periods(sw, [3.0])
        @test base.a ≈ sw_periods(sw, 3.0).a
        @test base.a_D ≈ sw_periods(sw, 3.0).a_D
        # path staying in the closed forms' common analyticity domain ℂ ∖ (−∞, 2Λ²]
        target = -3.0 + 3.0im
        res = continue_periods(sw, [3.0 + 0im, 3.0 + 3.0im, target])
        p = sw_periods(sw, target)
        @test res.a ≈ p.a rtol = 1e-8
        @test res.a_D ≈ p.a_D rtol = 1e-8
        @test res.da ≈ ExactWKB._da_du_closed(sw, target, Float64) rtol = 1e-7
        @test res.da_D ≈ ExactWKB._daD_du_closed(sw, target, Float64) rtol = 1e-7
    end

    @testset "monodromy loops reproduce the exact integer matrices" begin
        sw = SeibergWittenSU2()
        base_u = 3.0 + 0im
        circle(c, r, K) = [c + r * cis(2π * k / K) for k in 0:K]
        # transport matrix on (a_D, a) from a closed loop at base_u: S′ = S·Mᵀ
        function loop_monodromy(path)
            b = continue_periods(sw, [base_u])
            S = [b.a_D b.a; b.da_D b.da]
            f = continue_periods(sw, path)
            S2 = [f.a_D f.a; f.da_D f.da]
            Mnum = transpose(S \ S2)
            Mint = round.(Int, real.(Mnum))
            @test maximum(abs.(Mnum - Mint)) < 1e-6
            Mint
        end
        # counterclockwise loop around ∞ = clockwise big circle enclosing both points
        Minf = loop_monodromy(circle(0.0 + 0im, 3.0, 48))
        # lasso around the monopole point u = 2Λ²
        approach_m = [base_u, 2.6 + 0im]
        Mmon = loop_monodromy(vcat(approach_m, circle(2.0 + 0im, 0.6, 48)[2:end],
                                   reverse(approach_m)))
        # lasso around the dyon point u = −2Λ² (detour through the upper half-plane)
        approach_d = [base_u, 2.6im, -1.4 + 0im]
        Mdyon = loop_monodromy(vcat(approach_d, circle(-2.0 + 0im, 0.6, 48)[2:end],
                                    reverse(approach_d)))
        @test Mmon == sw_monodromy(sw, :monopole)
        @test Mdyon == sw_monodromy(sw, :dyon)
        @test Minf == sw_monodromy(sw, :infinity)
    end

    @testset "wall of marginal stability" begin
        sw = SeibergWittenSU2()
        n = 16
        wall = ms_wall(sw; n = n)
        @test length(wall) == 2n + 2
        @test wall[1] ≈ 2.0                                   # monopole endpoint
        @test wall[n+2] ≈ -2.0                                # dyon endpoint
        upper = wall[2:n+1]
        # a_D/a is real on the wall, with |a_D/a| ≤ 2 (0 at the monopole, ∓2 at the dyon)
        for u in upper
            r = sw_periods(sw, u).a_D / sw_periods(sw, u).a
            @test abs(imag(r)) < 1e-6
            @test abs(real(r)) < 2 + 1e-3
        end
        # u → −u symmetry (monopole ↔ dyon exchange): upper[k] = −conj(upper[n+1−k])
        for k in 1:n
            @test upper[k] ≈ -conj(upper[n+1-k]) rtol = 1e-8
        end
        # the wall lies strictly inside |u| = 2Λ² off the real axis
        @test all(abs(u) < 2.0 for u in upper)
    end

    @testset "chambers" begin
        sw = SeibergWittenSU2()
        @test sw_chamber(sw, 1.0) == :strong                  # real axis: |u| < 2Λ²
        @test sw_chamber(sw, -1.0) == :strong
        @test sw_chamber(sw, 3.0) == :weak
        @test sw_chamber(sw, -3.0) == :weak
        @test sw_chamber(sw, 0.3im) == :strong
        @test sw_chamber(sw, 2.5im) == :weak
        @test sw_chamber(sw, -0.2 - 0.3im) == :strong
        # Λ-scaling: the wall scales as Λ²
        sw2 = SeibergWittenSU2(Λ = 2)
        @test sw_chamber(sw2, 4.0 * 0.3im) == :strong
        @test sw_chamber(sw2, 4.0 * 2.5im) == :weak
    end

    @testset "errors" begin
        sw = SeibergWittenSU2()
        @test_throws PeriodError sw_periods(sw, -2.0)         # dyon branch point
        @test_throws PeriodError quantum_sw_periods(sw, 25.0; order = 2)
        @test_throws PeriodError quantum_sw_periods(sw, 2.0; order = 1)   # u = +2Λ² pinch
        @test_throws PeriodError quantum_sw_periods(sw, -2.0; order = 1)  # u = −2Λ² pinch
        @test_throws PeriodError central_charge(sw, 2.0, (1, 0); ħ = 0.1, order = 1)
        @test_throws PeriodError sw_monodromy(sw, :bogus)
        @test_throws PeriodError continue_periods(sw, [3.0, 2.0])  # through the singularity
        @test_throws Resurgence.InvalidArgument continue_periods(sw, Float64[])
        @test_throws Resurgence.InvalidArgument ms_wall(sw; n = 1)
    end
end
