using Resurgence: FormalSeries, LogSeries, MultiTransseries, coefficients,
                  power_offset, log_degree, log_block, resonant_solve,
                  resonance_lattice, sector, weight
using Base.MathConstants: eulergamma

# the g⁰ coefficient of a LogSeries block, honouring its power_offset (resonant_solve
# generally leaves blocks on different power grids - g⁰ is not always index 1)
function _g0(block)
    idx = 1 - Int(power_offset(block))
    c = coefficients(block)
    (1 ≤ idx ≤ length(c)) ? c[idx] : zero(eltype(c))
end

@testset "double_well" begin
    # the ZJJ double well V = q²(1−q)²/2 (ZJJ-I eq. 2.10b; ħ = g/√2, E_ours = g·E_zjj)
    zjjdw = SchrodingerProblem([0//1, 0//1, 1//2, -1//1, 1//2])
    # our house double well V = (z²−1)² (z = 2q−1 rescaling of the above)
    dwell = SchrodingerProblem([1, 0, -2, 0, 1])

    @testset "B/A series shape (ledger item 1)" begin
        setprecision(128) do
            E = big(1) / 64
            B = perturbative_b(zjjdw, E; order = 5, quad_rtol = 1e-15)
            A = instanton_a(zjjdw, E; order = 5, quad_rtol = 1e-15)
            @test power_offset(B) == -1 // 1 && power_offset(A) == -1 // 1
            # leading coefficients real and positive; even ħ-slots exactly zero
            @test real(B[0]) > 0 && abs(imag(B[0])) < 1e-20
            @test real(A[0]) > 0 && abs(imag(A[0])) < 1e-20
            @test iszero(B[1]) && iszero(A[1])
            # classical instanton action: A₋₁ = 2S_I = 2∫√(V−E) ≈ 1/3·(1/√2)·√2 = 1/3
            # at E → 0; at E = 1/64 it is somewhat smaller
            @test 0 < real(A[0]) < 1 / 3
        end
        # the layer refuses a single well
        @test_throws QuantizationError perturbative_b(
            SchrodingerProblem([0.0, 0.0, 1.0]; energy = 1.0), 1.0)
    end

    @testset "Dunne–Ünsal relation, ZJJ potential: c = 6 (ledger item 2)" begin
        setprecision(192) do
            r = verify_zjj(zjjdw, big(1) / 64; order = 11, quad_rtol = 1e-20)
            @test real(r.c) ≈ 6 rtol = 1e-8
            @test abs(imag(r.c)) < 1e-8
            @test r.orders == [3, 5, 7, 9]
            @test all(r.residuals .< 1e-8)     # ħ³, ħ⁵, … - each an independent theorem
        end
    end

    @testset "Dunne–Ünsal relation, V = (z²−1)²: c = 3/2" begin
        setprecision(160) do
            r = verify_zjj(dwell, big(1) / 2; order = 9, quad_rtol = 1e-20)
            @test real(r.c) ≈ 3 // 2 rtol = 1e-6
            @test all(r.residuals .< 1e-8)
        end
    end

    @testset "parity doublet vs diagonalization (pins the ± mapping)" begin
        # Deep double well (barrier top V(0)=1, E₀≈0.19): the exact parity condition
        # converges with WKB order - order 10 already reaches ~1e-10 vs the dense
        # diagonalization oracle. This is the end-to-end pin of the ½ factor and the
        # even-below-odd ordering (ledger item 3).
        ħ = 0.1
        E_even, E_odd = setprecision(140) do
            (wkb_eigenvalue(dwell, 0, big(ħ); parity = +1, order = 10, rtol = 1e-11,
                            quad_rtol = 1e-16),
             wkb_eigenvalue(dwell, 0, big(ħ); parity = -1, order = 10, rtol = 1e-11,
                            quad_rtol = 1e-16))
        end
        Ed = setprecision(140) do
            diagonalization_eigenvalues(dwell, ħ; nev = 2, N = 140, omega = 2.0)
        end
        @test E_even < E_odd                       # even (symmetric) member is lower
        @test E_even ≈ Ed[1] rtol = 1e-8
        @test E_odd ≈ Ed[2] rtol = 1e-8
        # energy_splitting is exactly the gap of the two parity roots
        ΔE = setprecision(140) do
            energy_splitting(dwell, 0, big(ħ); order = 10, rtol = 1e-11,
                             quad_rtol = 1e-16)
        end
        @test ΔE ≈ Ed[2] - Ed[1] rtol = 1e-4
        @test ΔE > 0
    end

    @testset "one-instanton splitting law (ZJJ eqs. 2.19–2.20, ledger item 3)" begin
        # ΔE_zjj = 2ξ(1 − (71/12)g − (6299/288)g² + O(g³)), ξ = e^{−1/6g}/√(πg);
        # in our variables ΔE_ours = g·ΔE_zjj at ħ = g/√2 (mapping derived in the
        # ledger). g = 0.03 keeps the ground state below the barrier top V=1/32
        # (E₀ ≈ 0.015). This is the pinning oracle for the ½ in the parity condition:
        # the WKB order does not change the median-summed splitting, so the residual
        # is the truncated g-series (O(g³) unincluded) plus the O(ξ) two-instanton
        # correction - both shrink as g → 0.
        g = 0.03
        ħ = big(g) / sqrt(big(2))
        ΔE = setprecision(160) do
            energy_splitting(zjjdw, 0, ħ; order = 10, rtol = 1e-9, quad_rtol = 1e-16)
        end
        ξ = exp(-1 / (6 * big(g))) / sqrt(big(π) * g)
        lead = g * 2ξ
        full = lead * (1 - big(71) / 12 * g - big(6299) / 288 * g^2)
        @test ΔE > 0
        # the exponential rate and leading prefactor are right to ~20%
        @test Float64(ΔE / lead) ≈ 1 rtol = 0.25
        # the g and g² corrections bring it to ~2%, and improve on the bare leading
        @test Float64(ΔE / full) ≈ 1 rtol = 3e-2
        @test abs(Float64(ΔE / full) - 1) < abs(Float64(ΔE / lead) - 1)
    end

    @testset "M6b consumer: ZJJ two-instanton log sector (ledger item 4)" begin
        # The exact ZJJ condition 1/Γ(½−B) + (εi/√2π)(−2/g)^B e^{−A/2} = 0 expanded
        # to ξ² must produce E^{(2)} = ξ²·(χ + γ), χ = ln(−2/g): e_{0,210} = 1,
        # e_{0,200} = γ (ZJJ-I eq. 2.19c). The log is minted by the resonant channel
        # (g²∂_g)u = g - Resurgence.resonant_solve at detuning 0 - and the algebra
        # runs through LogSeries blocks.
        T = Complex{BigFloat}
        setprecision(192) do
            γ = T(eulergamma)
            # log g from the resonant channel (the log-creating primitive of M6b)
            logg = resonant_solve(FormalSeries(T[0, 1], :g), 0)
            @test logg isa LogSeries && log_degree(logg) == 1
            # χ = ln(−2/g) = (ln 2 + iπ) − log g   (principal branch for g > 0)
            χ = LogSeries([FormalSeries(T[log(big(2)) + im * big(π)], :g)], :g) +
                (-1) * logg
            # b = B − ½ solves −b(1 − γb + …) = −εi ξ̂ e^{bχ}: with u_k the ξ̂^k
            # coefficient of b/(εi), the one-instanton block is u₁ = 1, so the
            # two-instanton block is u₂ = u₁²·(χ + γ) = χ + γ (the ½-plane digamma
            # ψ(1) = −γ from the Γ(½−B) expansion supplies the γ).
            u2 = χ + LogSeries([FormalSeries(T[γ], :g)], :g)
            # E-shift: ΔE⁽²⁾ = (εiξ̂)²·u₂ = −ξ̂²u₂ = +ξ²·u₂ (ξ̂² = −ξ²), so u₂ IS
            # the two-instanton block in units of ξ²: coefficient of χ and constant:
            @test log_degree(u2) == 1
            e210 = -_g0(log_block(u2, 1))                      # χ = … − log g
            e200 = _g0(log_block(u2, 0)) - e210 * (log(big(2)) + im * big(π))
            @test e210 ≈ 1 atol = 1e-30
            @test e200 ≈ γ atol = 1e-30
            # package as the resonant A/−A two-parameter transseries: inexact
            # actions need the lattice passed explicitly (the recorded M6b lesson)
            a = BigFloat(1) / 6
            @test_throws Resurgence.InvalidArgument resonance_lattice((a, -a))
            @test resonance_lattice((a, -a); lattice = [(1, 1)]) == [(1, 1)]
            pert = FormalSeries(T[1 // 2, -1, -9 // 2], :g)     # ZJJ-I eq. 2.6
            one_inst = FormalSeries(T[-1 / sqrt(big(π))], :g;
                                    power_offset = -1 // 2)      # −εξ at ε = +1
            two_inst = (1 / big(π)) * LogSeries(
                [FormalSeries(coefficients(log_block(u2, 0)), :g;
                              power_offset = -1 // 1),
                 FormalSeries(coefficients(log_block(u2, 1)), :g;
                              power_offset = -1 // 1)], :g)      # ξ²·u₂, ξ² = e^{−⅓g}/πg
            mt = MultiTransseries((a, -a), Dict((0, 0) => pert, (1, 0) => one_inst,
                                                (2, 0) => two_inst))
            @test weight(mt, (2, 0)) ≈ 1 / big(3)
            @test sector(mt, (2, 0)) isa LogSeries
            @test log_degree(sector(mt, (2, 0))) == 1
        end
    end

    @testset "layout errors" begin
        harm = SchrodingerProblem([0.0, 0.0, 1.0])
        @test_throws QuantizationError energy_splitting(harm, 0, 0.1)
        @test_throws QuantizationError verify_zjj(harm, 1.0; order = 5)
    end
end
