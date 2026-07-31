# GMN hyperkähler-metric oracles. The layer takes lattice-level BPS data only, so the
# fixtures are the cubic flagship spectrum (two states, A₂ charge lattice) and pure SU(2)
# in both chambers. Rung 1 checks the Darboux coordinates 𝒳_γ(ζ): the semiflat limit, the
# twistor reality constraint, multiplicativity in γ, the Kontsevich–Soibelman ray jump,
# the solver against the independent one-pass evaluator, and the conformal limit R → 0
# against the existing TBA layer.

@testset "hyperkahler" begin
    cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
    sp = bps_spectrum(cubic)
    P = signed_pairing(sp.basis)

    @testset "torus construction" begin
        t = gmn_torus(sp; R = 0.7, theta = [0.3, -0.8])
        @test n_states(t) == n_states(sp)
        @test n_charges(t) == 2
        @test radius(t) ≈ 0.7
        @test torus_angles(t) == [0.3, -0.8]
        # state pairing is the Gram form evaluated on the charges
        cs = charges(sp)
        @test t.state_pairing[1, 2] == transpose(cs[1]) * P * cs[2]
        @test all(t.state_pairing[a, a] == 0 for a in 1:n_states(t))
        # θ_γ is linear in the charge
        @test t.angles[1] ≈ sum(cs[1] .* [0.3, -0.8])
        # the BPSSpectrum method fills in signed_pairing
        @test gmn_torus(sp; R = 0.7).pairing == P

        @test_throws Resurgence.InvalidArgument gmn_torus(sp; R = -1)
        @test_throws Resurgence.InvalidArgument gmn_torus(sp; theta = [0.1])
        @test_throws Resurgence.InvalidArgument gmn_torus(sp; sigma = 0)
        @test_throws Resurgence.InvalidArgument gmn_torus(sp.states, [0 1 0; -1 0 0; 0 0 0])
        @test_throws TBAError gmn_torus(BPSState{Float64}[], P)
        # basis central charges: inferred by default, refused when undetermined,
        # and checked against the spectrum when supplied
        cs = charges(sp)
        @test sum(cs[1] .* t.basisZ) ≈ sp.states[1].central_charge
        @test charge_central_charge(t, cs[2]) ≈ sp.states[2].central_charge
        @test_throws Resurgence.InvalidArgument gmn_torus([sp.states[1]], P)
        @test_throws Resurgence.InvalidArgument gmn_torus(sp; basis_Z = [1.0 + 0im])
        @test_throws Resurgence.InvalidArgument gmn_torus(sp; basis_Z = [1.0 + 0im, 2im])
    end

    @testset "vanishing pairing is exactly semiflat" begin
        # a single state corrects nothing on its own ray direction (⟨γ,γ⟩ = 0), so every
        # multiple of its charge stays semiflat - while a transverse charge does not
        bZ = gmn_torus(sp).basisZ
        t = gmn_torus([sp.states[1]], P; R = 1.3, theta = [0.4, 0.9], basis_Z = bZ)
        sol = solve_gmn(t)
        @test n_iterations(sol) == 1
        γ = t.charges[1]
        for ζ in (0.3 + 0.2im, -1.1 + 0.4im, 0.05im)
            @test xi_value(sol, γ, ζ) ≈ semiflat_xi(t, γ, ζ)
            @test xi_value(sol, -2 * γ, ζ) ≈ semiflat_xi(t, -2 * γ, ζ)
        end
        δ = [γ[2], -γ[1]]                       # ⟨δ, γ⟩ ≠ 0: instanton-corrected
        @test transpose(δ) * P * γ != 0
        @test !isapprox(xi_value(sol, δ, 0.3 + 0.2im), semiflat_xi(t, δ, 0.3 + 0.2im))
    end

    @testset "semiflat multiplicativity" begin
        t = gmn_torus(sp; R = 1.0, theta = [0.5, 1.2])
        ζ = 0.4 - 0.3im
        @test semiflat_xi(t, [1, 0], ζ) * semiflat_xi(t, [0, 1], ζ) ≈
              semiflat_xi(t, [1, 1], ζ)
    end

    t = gmn_torus(sp; R = 0.9, theta = [0.35, -0.6])
    sol = solve_gmn(t)

    @testset "solver" begin
        @test residual(sol) < 1e-8
        @test n_iterations(sol) < 500
        @test size(sol.gplus) == (length(sol.grid), n_states(t))
        # the window is symmetric (ledger item 2)
        @test sol.grid[1] ≈ -sol.grid[end]
        # warm start on the same grid converges immediately
        sol2 = solve_gmn(t; window = (sol.grid[1], sol.grid[end]),
                         n_points = length(sol.grid), seed = sol)
        @test n_iterations(sol2) ≤ 2
        @test maximum(abs, sol2.gplus - sol.gplus) < 1e-7
        @test_throws Resurgence.InvalidArgument solve_gmn(t; relax = 0)
        @test_throws Resurgence.InvalidArgument solve_gmn(t; seed = sol, n_points = 33)
        @test_throws TBAError solve_gmn(t; maxiter = 1, tol = 1e-14)
    end

    @testset "solver vs one-pass evaluation" begin
        # `_log_xi` re-does the integral with pole subtraction and a closed-form PV,
        # a completely different code path from the solver's precomputed kernel blocks
        for a in 1:n_states(t), k in (10, 40, 70)
            θ = sol.grid[k]
            ζ = exp(-θ) * t.rays[a]
            @test ExactWKB._log_xi(sol, t.charges[a], ζ, 0) ≈ sol.gplus[k, a] atol = 1e-6
            @test ExactWKB._log_xi(sol, t.charges[a], -ζ, 0) ≈ sol.gminus[k, a] atol = 1e-6
        end
    end

    @testset "twistor reality constraint" begin
        # 𝒳_γ(−1/ζ̄)* = 𝒳_{−γ}(ζ) - the real structure on the twistor space
        for ζ in (0.5 + 0.4im, -0.7 + 0.2im, 0.3 - 0.9im)
            for γ in ([1, 0], [0, 1], [1, 1])
                lhs = conj(xi_value(sol, γ, -1 / conj(ζ)))
                rhs = xi_value(sol, -γ, ζ)
                @test lhs ≈ rhs rtol = 1e-6
            end
        end
    end

    @testset "multiplicativity in the charge" begin
        ζ = 0.45 + 0.35im
        for (γ, γp) in (([1, 0], [0, 1]), ([1, 1], [1, -1]), ([2, 0], [-1, 3]))
            @test ExactWKB._log_xi(sol, γ, ζ, 0) + ExactWKB._log_xi(sol, γp, ζ, 0) ≈
                  ExactWKB._log_xi(sol, γ + γp, ζ, 0) rtol = 1e-8
        end
    end

    @testset "Kontsevich–Soibelman ray jump" begin
        # crossing ℓ_{γ_b}, 𝒳_γ picks up (1 − σ𝒳_{γ_b})^{Ω⟨γ,γ_b⟩}
        b = 1
        θ = sol.grid[fld(length(sol.grid), 2)]
        ζ = exp(-θ) * t.rays[b]
        γ = t.charges[2]
        c = transpose(γ) * t.pairing * t.charges[b]
        @test c != 0
        Xb = xi_value(sol, t.charges[b], ζ)
        jump = xi_value(sol, γ, ζ; side = :plus) / xi_value(sol, γ, ζ; side = :minus)
        @test jump ≈ (1 - t.sigma[b] * Xb)^(t.omega[b] * c) rtol = 1e-6
    end

    @testset "conformal limit R → 0 reproduces the TBA layer" begin
        # ζ = πRħ turns the semiflat term πRZ/ζ + πRZ̄ζ into Z/ħ + O(R²), so the
        # finite-R Darboux coordinates collapse onto the conformal-limit Voros symbols
        R = 1e-3
        tc = gmn_torus(sp; R = R)
        solc = solve_gmn(tc; tol = 1e-10)
        tba = solve_tba(sp)
        for ħ in (0.35, 0.5)
            for a in 1:n_states(sp)
                ζ = π * R * ħ * cis(0.4)
                lhs = ExactWKB._log_xi(solc, tc.charges[a], ζ, 0)
                rhs = log(voros_value(tba, a, ħ * cis(0.4)))
                @test lhs ≈ rhs rtol = 5e-3
            end
        end
    end

    # ── rung 2: the metric ────────────────────────────────────────────────────────────

    sw = SeibergWittenSU2()
    u0 = 1.0 + 0.3im                      # strong-coupling chamber
    θ0 = [0.4, 0.7]

    @testset "twistor symplectic form" begin
        f = (u, th) -> su2_torus(sw, u, th; R = 3.0, chamber = :strong)
        mp = metric_point(f, u0, θ0)
        d = hk_diagnostics(mp)
        # ϖ(ζ) is a Laurent polynomial of degrees −1, 0, 1 only: the twistor theorem.
        # The ζ^{±2} terms cancel because da_D = τ da (special geometry).
        @test d.residual < 1e-7
        # the hyperkähler triple is real, and the metric symmetric and positive definite
        @test d.imaginary < 1e-5
        @test d.asymmetry < 1e-5
        @test all(>(0), d.eigenvalues)

        # ϖ itself is antisymmetric at every ζ
        for ζ in (0.7 + 0.5im, -1.1 + 0.3im)
            ϖ = holomorphic_symplectic_form(mp, ζ)
            @test maximum(abs, ϖ + transpose(ϖ)) < 1e-9 * maximum(abs, ϖ)
        end

        @test_throws Resurgence.InvalidArgument symplectic_expansion(mp; zetas = [1.0])
        @test_throws Resurgence.InvalidArgument metric_point(f, u0, [0.1])
    end

    @testset "semiflat limit" begin
        # instanton corrections are O(exp(−2πR|Z|)); at R = 3 the lightest state has
        # 2πR|Z| ≈ 10, so hk_metric must sit on semiflat_metric to ~1e-4
        R = 3.0
        f = (u, th) -> su2_torus(sw, u, th; R, chamber = :strong)
        t0 = f(u0, θ0)
        dd = sw_period_derivatives(sw, u0)
        gsf = semiflat_metric(t0, [dd.da_D, dd.da])
        g = hk_metric(f, u0, θ0)
        @test maximum(abs, g - gsf) < 1e-4 * maximum(abs, gsf)
        # the correction is of the predicted instanton size, not zero
        @test maximum(abs, g - gsf) > 1e-9 * maximum(abs, gsf)

        # τ = da_D/da lives in the upper half plane and sets the fibre metric
        τ = dd.da_D / dd.da
        @test imag(τ) > 0
        @test gsf[1, 1] ≈ R * imag(τ) * abs2(dd.da)
        @test gsf[3, 3] * gsf[4, 4] - gsf[3, 4]^2 ≈ 1 / (4π^2 * R)^2

        @test_throws Resurgence.InvalidArgument semiflat_metric(t0, [1.0])
        # a period pair with Im τ ≤ 0 is not a Riemannian semiflat metric
        @test_throws TBAError semiflat_metric(t0, [1.0 + 0im, 1.0 + 0im])
    end

    @testset "instanton corrections at small radius" begin
        # at R = 0.4 the corrections are O(e^{−1.4}) - the metric must depart from
        # semiflat well beyond numerical noise, and stay Riemannian
        R = 0.4
        f = (u, th) -> su2_torus(sw, u, th; R, chamber = :strong)
        t0 = f(u0, θ0)
        dd = sw_period_derivatives(sw, u0)
        gsf = semiflat_metric(t0, [dd.da_D, dd.da])
        mp = metric_point(f, u0, θ0)
        g = hk_metric(mp)
        d = hk_diagnostics(mp)
        @test d.residual < 1e-6
        @test all(>(0), d.eigenvalues)
        @test maximum(abs, g - gsf) > 1e-3 * maximum(abs, gsf)
    end

    # ── rung 3: Ooguri-Vafa and SU(2) wall crossing ───────────────────────────────────

    @testset "Ooguri-Vafa" begin
        z = 0.25 + 0.1im
        θov = [0.6, 1.1]
        R = 0.8
        t = ooguri_vafa_torus(z, θov; R)
        @test n_states(t) == 1
        sol = solve_gmn(t)

        # the lone state's own coordinate is untouched; the dual one is corrected
        for ζ in (0.6 + 0.5im, -0.9 + 0.4im)
            @test xi_value(sol, [0, 1], ζ) ≈ semiflat_xi(t, [0, 1], ζ)
            @test !isapprox(xi_value(sol, [1, 0], ζ), semiflat_xi(t, [1, 0], ζ))
        end

        # adaptive quadrature of the same integral - a completely independent solver
        for ζ in (0.6 + 0.5im, -0.9 + 0.4im, 1.7 - 1.2im)
            for γ in ([1, 0], [1, 1], [2, -1])
                @test xi_value(sol, γ, ζ) ≈ ooguri_vafa_xi(t, γ, ζ) rtol = 1e-6
            end
        end
        @test_throws TBAError ooguri_vafa_xi(t, [1, 0], t.rays[1])
        @test_throws Resurgence.InvalidArgument ooguri_vafa_xi(gmn_torus(sp), [1, 0], 1.0)

        # the ζ → ∞ limit is the closed-form Bessel instanton sum, which pins the 2πR
        # in the semiflat exponent
        for γ in ([1, 0], [2, 1])
            ζ = 1.0e5 * cis(0.31)
            corr = instanton_correction(sol, γ, ζ)
            @test real(corr) ≈ ooguri_vafa_instantons(t, γ) rtol = 1e-5
            @test abs(imag(corr)) < 1e-5
        end
        # a charge with no pairing against the state gets no correction at all
        @test ooguri_vafa_instantons(t, [0, 3]) == 0

        # the one-instanton term alone already dominates: sin θ_s K₀(2πR|z|)/π
        one_inst = -1 * 1 * sin(θov[2]) * ExactWKB._besselk0(2π * R * abs(z)) / π
        @test ooguri_vafa_instantons(t, [1, 0]) ≈ one_inst rtol = 0.15
        @test ExactWKB._besselk0(1.0) ≈ 0.4210244382407083 rtol = 1e-8
        @test ExactWKB._besselk0(2.5) ≈ 0.0623475532003457 rtol = 1e-8
    end

    @testset "Ooguri-Vafa metric" begin
        # the local model is Riemannian and instanton-corrected; the corrections grow as
        # the state gets light (|z| → 0 at fixed R), where semiflat alone degenerates
        θov = [0.6, 1.1]
        R = 1.0
        f = (z, th) -> ooguri_vafa_torus(z, th; R)
        gaps = Float64[]
        for z in (0.30 + 0.05im, 0.10 + 0.02im)
            t = f(z, θov)
            mp = metric_point(f, z, θov)
            d = hk_diagnostics(mp)
            @test d.residual < 1e-6
            @test all(>(0), d.eigenvalues)
            gsf = semiflat_metric(t, ooguri_vafa_period_derivatives(z))
            g = hk_metric(mp)
            push!(gaps, maximum(abs, g - gsf) / maximum(abs, gsf))
        end
        @test gaps[2] > gaps[1]
    end

    @testset "SU(2) metric: chamber independence across the wall" begin
        # The Kontsevich–Soibelman identity between the two SU(2) chambers says the two
        # ray products agree, so the GMN equations - and the metric they build - do not
        # care which chamber's spectrum is used. Evaluated just outside the wall of
        # marginal stability, where the BPS rays are still distinct (exactly on the wall
        # every ray coincides and the kernel is singular).
        sw = SeibergWittenSU2()
        uw = ms_wall(sw; n = 64)[9]
        u = uw * 1.05
        R = 1.2
        θw = [0.5, 0.9]
        gs = hk_metric((x, th) -> su2_torus(sw, x, th; R, chamber = :strong), u, θw)
        gw = hk_metric((x, th) -> su2_torus(sw, x, th; R, chamber = :weak, tower = 4),
                       u, θw)
        @test maximum(abs, gs - gw) < 1e-6 * maximum(abs, gs)

        # and the agreement is not trivial: both are ~1.5% away from the semiflat metric
        t = su2_torus(sw, u, θw; R, chamber = :strong)
        dd = sw_period_derivatives(sw, u)
        gsf = semiflat_metric(t, [dd.da_D, dd.da])
        @test maximum(abs, gs - gsf) > 1e-3 * maximum(abs, gsf)
        @test maximum(abs, gw - gsf) > 1e-3 * maximum(abs, gsf)

        # the ten-state weak-chamber solution is a bona fide hyperkähler metric too
        mp = metric_point((x, th) -> su2_torus(sw, x, th; R, chamber = :weak, tower = 4),
                          u, θw)
        d = hk_diagnostics(mp)
        @test d.residual < 1e-7
        @test d.imaginary < 1e-5
        @test d.asymmetry < 1e-5
        @test all(>(0), d.eigenvalues)
    end

    @testset "chamber independence under the NS deformation" begin
        # Does the previous testset's result survive when the classical central charges
        # are replaced by the Nekrasov–Shatashvili quantum periods? There is no theorem
        # either way: KS wall-crossing is a statement about the classical Z, and the wall
        # itself moves under the deformation.
        #
        # The answer is yes, and structurally so: the GMN equations see the BPS data only
        # through (charges, Ω, Z), and the KS identity that makes the two chamber spectra
        # interchangeable lives at the charge/Ω level, which the deformation does not
        # touch. What the deformation actually does is move the wall (see the deformed-wall
        # testset in test_sw_curve.jl), i.e. relabel which chamber a given u belongs to.
        sw = SeibergWittenSU2()
        u = ms_wall(sw; n = 64)[9] * 1.05
        R = 1.2
        θw = [0.5, 0.9]
        strong(h, o) = (x, th) -> su2_torus(sw, x, th; R, chamber = :strong, ħ = h, order = o)
        weak(h, o) = (x, th) -> su2_torus(sw, x, th; R, chamber = :weak, tower = 4,
                                          ħ = h, order = o)

        # the gate: ħ = 0 at order 1 must reproduce the classical metric exactly. If this
        # slips, nothing below means anything.
        gs0 = hk_metric(strong(0, 0), u, θw)
        @test hk_metric(strong(0.0, 1), u, θw) == gs0

        # the headline: chamber independence holds at the same 1e-6 relative bar the
        # classical result meets, at every ħ we can trust the order-1 truncation at, and
        # at order 2 as well
        for (h, o) in ((0.1, 1), (0.3, 1), (0.6, 1), (0.2, 2))
            gs = hk_metric(strong(h, o), u, θw)
            gw = hk_metric(weak(h, o), u, θw)
            @test maximum(abs, gs - gw) < 1e-6 * maximum(abs, gs)
        end

        # and it is not vacuous: the deformed metric genuinely moves, by O(ħ²) - three
        # orders of magnitude more than the chamber discrepancy it leaves untouched
        move(h) = maximum(abs, hk_metric(strong(h, 1), u, θw) - gs0) / maximum(abs, gs0)
        m1, m2 = move(0.1), move(0.2)
        @test m1 > 1e-4
        @test isapprox(log(m2 / m1) / log(2), 2; atol = 0.05)

        # deformed special geometry still holds: ϖ(ζ) stays Laurent of degrees −1, 0, 1,
        # and the deformed metric is still a bona fide hyperkähler one
        mpd = metric_point(weak(0.3, 1), u, θw)
        @test symplectic_expansion(mpd).residual < 1e-7
        dd = hk_diagnostics(mpd)
        @test dd.residual < 1e-7
        @test dd.imaginary < 1e-5
        @test dd.asymmetry < 1e-5
        @test all(>(0), dd.eigenvalues)
    end
end
