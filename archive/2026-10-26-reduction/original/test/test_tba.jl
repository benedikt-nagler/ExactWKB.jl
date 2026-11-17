# TBA-layer oracles. The flagship identity: the TBA-side Voros values, computed from
# the BPS data alone (Z_γ, Ω, ⟨,⟩ - no WKB series), must agree with the Borel–Padé
# median sums of the WKB Voros symbols wherever both converge. Fixtures are the cubic
# (two-state chamber) and the symmetric quartic, as in test_ddp.jl / test_bps.jl.
# Series-side contours are rebuilt wide (margin 0.4, as in test_ddp.jl) because the
# chamber-frame basis contours are too tight for the order-12 even-period roundoff
# check; orientation is matched to the state by its central charge.

import Resurgence
using Resurgence: coefficients

# a Voros symbol carrying the state charge Z: try the chamber basis contours and all
# wide turning-point-pair contours; return (symbol, sign) with sign·∮√Q ≈ Z, or
# nothing if no candidate matches
function _state_symbol(sp, w, Z)
    prob = sp.basis.problem
    tps = simple_turning_points(prob)
    cands = collect(sp.basis.contours)
    for i in 1:(length(tps) - 1), j in (i + 1):length(tps)
        try
            push!(cands, encircling_contour(tps[i], tps[j]; margin = 0.4))
        catch
        end
    end
    for c in cands
        cp = try
            period_integral(prob, c)
        catch
            continue
        end
        s = abs(cp - Z) < 1e-6 * abs(Z) ? 1 :
            abs(cp + Z) < 1e-6 * abs(Z) ? -1 : 0
        s == 0 && continue
        v = try
            voros_symbol(w, c)
        catch
            continue
        end
        return v, s
    end
    nothing
end

# Borel–Padé log-Voros of the state with central charge Z, in the physical frame
function _borel_log_voros(sp, w, Z, ħ; side = :median)
    sym = _state_symbol(sp, w, Z)
    sym === nothing && return nothing
    v, s = sym
    s * ExactWKB._log_voros(v, ħ; side)
end

@testset "tba" begin
    cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
    sp = bps_spectrum(cubic)
    sol = solve_tba(sp)
    sys = sol.system
    ms = [mass(s) for s in sp.states]
    w12 = wkb_expansion(cubic; order = 12)

    @testset "system extraction" begin
        @test n_states(sys) == 2
        @test abs(sys.pairing[1, 2]) == 1
        @test sys.pairing[1, 1] == sys.pairing[2, 2] == 0
        @test sys.pairing[1, 2] == -sys.pairing[2, 1]
        @test all(abs(abs(r) - 1) < 1e-12 for r in sys.rays)
        # every state charge is carried by a turning-point cycle
        for a in 1:2
            @test _state_symbol(sp, w12, sys.Z[a]) !== nothing
        end
    end

    @testset "solver basics" begin
        @test n_iterations(sol) < 500
        @test residual(sol) < 2e-10
        # the state with real central charge has a real symbol on its own ray
        a = argmin([abs(imag(z)) for z in sys.Z])
        mid = length(sol.grid) ÷ 2
        @test abs(imag(sol.gplus[mid, a])) < 1e-8
    end

    @testset "IR plateau = constant Y-system (golden ratio)" begin
        # read the plateau well above the truncated IR boundary (its own error is
        # ~e^{-(θ-θmin)}) and well below the mass transition at θ ≈ -log m
        k = findfirst(θ -> θ > sol.grid[1] + 12, sol.grid)
        @test sol.grid[k] < -3
        φ = (1 + sqrt(5)) / 2
        for a in 1:2
            @test abs(exp(sol.gplus[k, a]) - φ) < 5e-3
        end
    end

    @testset "UV: semiflat + first quantum correction" begin
        for a in 1:2
            v, s = _state_symbol(sp, w12, sys.Z[a])
            v1 = s * coefficients(quantum_series(v))[1]
            ħc = ms[a] / 40 * -sys.rays[a]   # on the anti-ray: off every cubic pole
            tail = ExactWKB._log_x(sol, a, ħc, 0) - sys.Z[a] / ħc
            @test abs(tail - v1 * ħc) / abs(v1 * ħc) < 5e-2
        end
    end

    @testset "TBA vs Borel–Padé median (the flagship oracle)" begin
        for a in 1:2, x in (6.0, 8.0)
            ħ = ms[a] / x
            lb = _borel_log_voros(sp, w12, sys.Z[a], ħ)
            lt = ExactWKB._log_x(sol, a, ħ, 0)
            tail_b = lb - sys.Z[a] / ħ
            tail_t = lt - sys.Z[a] / ħ
            @test abs(tail_t - tail_b) / abs(tail_b) < 1e-3
            @test voros_value(sol, a, ħ) ≈ exp(lb) rtol = 1e-3
        end
    end

    @testset "lateral jump = DDP, Stokes constant an exact integer" begin
        # state i has its wall at θ_c = π/2; evaluate at real ħ on the θ_c = 0 ray
        i = argmax([abs(imag(z)) for z in sys.Z])
        j = 3 - i
        ħ = ms[j] / 6
        jump = ExactWKB._log_x(sol, i, ħ, 1) - ExactWKB._log_x(sol, i, ħ, -1)
        Xwall = voros_value(sol, j, ħ)
        κ = jump / log1p(Xwall)
        @test abs(κ - sys.pairing[i, j] * sys.omega[j]) < 1e-2
        # and the jump agrees with the measured DDP jump of the WKB symbols
        # (the test_ddp.jl cycle pair: γ₁ decaying on θ_c = 0, γ₂ the jumping one)
        tps = simple_turning_points(cubic)
        c1 = reverse(encircling_contour(tps[1], tps[2]; margin = 0.4))
        c2 = encircling_contour(tps[2], tps[3]; margin = 0.4)
        v1c = voros_symbol(w12, c1)
        v2c = voros_symbol(w12, c2)
        κp = intersection_pairing(cubic, c2, c1)
        r = verify_ddp(v2c, v1c, κp, ħ; theta = 0.0)
        s = abs(classical_period(v2c) - sys.Z[i]) < 1e-6 ? 1 : -1
        @test s * r.jump_measured ≈ jump rtol = 5e-2
    end

    @testset "chamber independence (wall-crossing invariance)" begin
        sp2 = bps_spectrum(cubic; theta = 0.3)
        sol2 = solve_tba(sp2)
        sys2 = sol2.system
        for a in 1:n_states(sys)
            a2 = argmin([abs(z - sys.Z[a]) for z in sys2.Z])
            @test sys2.Z[a2] ≈ sys.Z[a] rtol = 1e-6
            ħ = ms[a] / 6
            @test voros_value(sol2, a2, ħ) ≈ voros_value(sol, a, ħ) rtol = 1e-3
        end
    end

    @testset "harmonic: the one-state TBA is semiflat-exact" begin
        harm = SchrodingerProblem([-1.0, 0.0, 1.0])
        sph = bps_spectrum(harm)
        @test n_states(sph) == 1
        solh = solve_tba(sph)
        @test n_iterations(solh) ≤ 2
        wh = wkb_expansion(harm; order = 6)
        Zh = solh.system.Z[1]
        for ħ in (0.1, 0.3)
            lb = _borel_log_voros(sph, wh, Zh, ħ)
            @test ExactWKB._log_x(solh, 1, ħ, 0) ≈ lb rtol = 1e-8
        end
    end

    @testset "quartic flagship" begin
        quartic = SchrodingerProblem([-1.0, 0.0, 0.0, 0.0, 1.0])
        spq = bps_spectrum(quartic)
        solq = solve_tba(spq)
        wq = wkb_expansion(quartic; order = 12)
        msq = [mass(s) for s in spq.states]
        @test residual(solq) < 2e-10
        tested = 0
        for a in 1:n_states(solq.system)
            ħ = msq[a] / 7
            lb = _borel_log_voros(spq, wq, solq.system.Z[a], ħ)
            lb === nothing && continue
            lt = ExactWKB._log_x(solq, a, ħ, 0)
            tail_b = lb - solq.system.Z[a] / ħ
            tail_t = lt - solq.system.Z[a] / ħ
            @test abs(tail_t - tail_b) / abs(tail_b) < 3e-2
            tested += 1
        end
        @test tested ≥ 2
    end

    @testset "errors" begin
        airy = SchrodingerProblem([0.0, 1.0])
        spa = bps_spectrum(airy)
        @test n_states(spa) == 0
        @test_throws TBAError solve_tba(spa)
        @test_throws TBAError solve_tba(sp; maxiter = 1)
        @test_throws Resurgence.InvalidArgument solve_tba(sp; relax = 0.0)
        @test_throws Resurgence.InvalidArgument voros_value(sol, 1, 0.1; side = :up)
        @test_throws Resurgence.InvalidArgument voros_value(sol, 7, 0.1)
        # on a BPS ray but far outside the rapidity window
        @test_throws TBAError voros_value(sol, 1, 1e-12im)
    end
end
