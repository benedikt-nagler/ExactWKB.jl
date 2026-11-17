# Oracles for the implicit differentiation of the GMN fixed point. Three of them are exact
# and cost one solve: the single-state case is semiflat in closed form, the derivative is
# linear in the charge, and a common rotation of every central charge is a symmetry of the
# whole system. The quantitative ones are Ooguri-Vafa against an independent quadrature,
# and pure SU(2) against the central difference this layer replaces.

@testset "gmn_derivative" begin
    cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
    sp = bps_spectrum(cubic)
    P = signed_pairing(sp.basis)

    # branch-safe log 𝒳_γ(ζ) from the public API: the instanton correction never wraps
    function logxi(t, sol, γ, ζ)
        Zγ = charge_central_charge(t, γ)
        θγ = sum(γ .* torus_angles(t))
        instanton_correction(sol, γ, ζ) +
        π * radius(t) * (Zγ / ζ + conj(Zγ) * ζ) + im * θγ
    end

    @testset "argument checks" begin
        t = gmn_torus(sp; R = 0.9, theta = [0.35, -0.6])
        sol = solve_gmn(t)
        dZ = [1.0 + 0.0im, 0.3 - 0.2im]
        @test_throws Resurgence.InvalidArgument solve_gmn_derivative(sol, [1.0 + 0im])
        @test_throws Resurgence.InvalidArgument solve_gmn_derivative(sol, dZ; relax = 0)
        @test_throws TBAError solve_gmn_derivative(sol, dZ; maxiter = 1, tol = 1e-30)
        gd = solve_gmn_derivative(sol, dZ)
        @test n_parameters(gd) == 4
        @test n_charges(gd) == 2
        @test n_states(gd) == n_states(sp)
        @test residual(gd) < 1e-9
        @test n_iterations(gd) ≥ 1
        @test_throws Resurgence.InvalidArgument log_xi_derivative(gd, [1], 0.4 + 0.2im)
        @test_throws TBAError log_xi_derivative(gd, [1, 0], 0)
        # on a BPS ray it pairs with, the coordinate jumps and the derivative is one-sided
        b = findfirst(c -> !iszero(transpose([1, 0]) * t.pairing * c), t.charges)
        @test_throws TBAError log_xi_derivative(gd, [1, 0], 0.5 * t.rays[b])
        # a ray it does not pair with is harmless: 𝒳_{[1,0]} does not jump there
        s = findfirst(c -> iszero(transpose([1, 0]) * t.pairing * c), t.charges)
        @test log_xi_derivative(gd, [1, 0], 0.5 * t.rays[s]) isa Vector{ComplexF64}
    end

    @testset "the single-state case is semiflat in closed form" begin
        # ⟨γ,γ⟩ = 0 kills every correction to 𝒳_γ, so its derivative is the derivative of
        # πR(Z_γ/ζ + Z̄_γζ) + iθ_γ, term by term and with no discretization at all
        bZ = gmn_torus(sp).basisZ
        R = 1.3
        t = gmn_torus([sp.states[1]], P; R, theta = [0.4, 0.9], basis_Z = bZ)
        dZ = [0.7 + 0.2im, -0.4 + 1.1im]
        gd = solve_gmn_derivative(solve_gmn(t), dZ)
        γ = t.charges[1]
        for c in (1, -2), ζ in (0.3 + 0.2im, -1.1 + 0.4im, 0.05im)
            g = c * γ
            Z = sum(g .* dZ)
            closed = [π * R * (Z / ζ + conj(Z) * ζ),
                      π * R * (im * Z / ζ + conj(im * Z) * ζ),
                      im * g[1], im * g[2]]
            @test log_xi_derivative(gd, g, ζ) ≈ closed atol = 1e-14
        end
    end

    t = gmn_torus(sp; R = 0.9, theta = [0.35, -0.6])
    sol = solve_gmn(t)

    @testset "linearity in the charge" begin
        gd = solve_gmn_derivative(sol, [1.0 + 0.0im, 0.3 - 0.2im])
        ζ = 0.45 + 0.7im
        @test log_xi_derivative(gd, [1, 0], ζ) + log_xi_derivative(gd, [0, 1], ζ) ≈
              log_xi_derivative(gd, [1, 1], ζ)
        @test 3 * log_xi_derivative(gd, [1, -1], ζ) ≈ log_xi_derivative(gd, [3, -3], ζ)
    end

    @testset "a common rotation of Z is a symmetry" begin
        # Z → e^{iα}Z rotates every ray with the ζ-plane, so 𝒳^{(α)}_γ(ζ) = 𝒳_γ(e^{-iα}ζ)
        # and the generator is ∂_α log 𝒳_γ(ζ) = −iζ ∂_ζ log 𝒳_γ(ζ). Taking dZ/du = iZ puts
        # that rotation in the Re u direction. Nothing else pins the ray-motion term of
        # the kernel on its own, since the source term carries only differences λ_a − λ_b.
        gd = solve_gmn_derivative(sol, im * t.basisZ)
        # the ray data itself does not move: the sample points rotate with the rays
        @test maximum(abs, gd.dgplus[:, :, 1]) < 1e-12
        @test maximum(abs, gd.dgminus[:, :, 1]) < 1e-12
        for γ in ([1, 0], [1, 1]), ζ in (0.6 + 0.5im, -0.9 + 0.3im)
            h = 1e-6
            dζ = (logxi(t, sol, γ, ζ * (1 + h)) - logxi(t, sol, γ, ζ * (1 - h))) /
                 (2 * h * ζ)
            @test log_xi_derivative(gd, γ, ζ)[1] ≈ -im * ζ * dζ atol = 1e-7
        end
    end

    @testset "Ooguri-Vafa against an independent quadrature" begin
        # `ooguri_vafa_xi` is adaptive quadrature of the same integral, independent of the
        # solver in grid, weights, window and iteration. Differencing it in each of the
        # four real directions is therefore an oracle for the whole derivative layer.
        z = 0.3 + 0.05im
        θov = [0.6, 1.1]
        f = (x, th) -> ooguri_vafa_torus(x, th; R = 1.0)
        gd = solve_gmn_derivative(solve_gmn(f(z, θov)), ooguri_vafa_period_derivatives(z))
        ζ = 0.7 * cis(1.0)
        h = 1e-5
        for γ in ([1, 0], [1, 2])
            q(x, th) = log(ooguri_vafa_xi(f(x, th), γ, ζ))
            fd = [(q(z + h, θov) - q(z - h, θov)) / (2h),
                  (q(z + im * h, θov) - q(z - im * h, θov)) / (2h),
                  (q(z, θov + [h, 0]) - q(z, θov - [h, 0])) / (2h),
                  (q(z, θov + [0, h]) - q(z, θov - [0, h])) / (2h)]
            @test log_xi_derivative(gd, γ, ζ) ≈ fd atol = 1e-8
        end
    end

    # ── the metric layer: the exact route against the central difference it replaces ────

    sw = SeibergWittenSU2()
    uw = ms_wall(sw; n = 64)[9] * 1.05
    Rw = 1.2
    θw = [0.5, 0.9]
    strong = (x, th) -> su2_torus(sw, x, th; R = Rw, chamber = :strong)
    dzdu(x) = (d = sw_period_derivatives(sw, x); [d.da_D, d.da])

    @testset "SU(2) metric: exact against the finite difference" begin
        mpfd = metric_point(strong, uw, θw)
        mpex = metric_point(strong, uw, θw; dZdu = dzdu(uw))
        @test !is_exact(mpfd)
        @test is_exact(mpex)
        @test isempty(mpex.plus) && isempty(mpex.minus)
        gfd, gex = hk_metric(mpfd), hk_metric(mpex)
        @test maximum(abs, gfd - gex) < 1e-6 * maximum(abs, gex)

        # the central difference converges to the exact derivative at exactly its own
        # order: halving h divides the gap by four. This is the Richardson oracle.
        errs = [maximum(abs, hk_metric(strong, uw, θw; h) - gex) for h in (4e-3, 2e-3, 1e-3)]
        @test all(isapprox(errs[i] / errs[i + 1], 4; atol = 0.05) for i in 1:2)
        @test errs[3] < 1e-7 * maximum(abs, gex)
    end

    @testset "the twistor oracles stop being h-limited" begin
        # `symplectic_expansion`'s residual is the numerical test that the 𝒳_γ are Darboux
        # coordinates, and the reality of the hyperkähler triple is the test that ω₁,ω₂,ω₃
        # are a triple. Both were bounded by h, not by the theory: with exact derivatives
        # they land at the solver's own precision, four to six orders lower.
        dfd = hk_diagnostics(metric_point(strong, uw, θw))
        dex = hk_diagnostics(metric_point(strong, uw, θw; dZdu = dzdu(uw)))
        @test dex.residual < 1e-12
        @test dex.imaginary < 1e-10
        @test dex.asymmetry < 1e-10
        @test dex.residual < 1e-3 * dfd.residual
        @test dex.imaginary < 1e-3 * dfd.imaginary
        @test dex.asymmetry < 1e-3 * dfd.asymmetry
        @test all(>(0), dex.eigenvalues)
    end

    @testset "M6b headlines through the exact route" begin
        # chamber independence across the wall of marginal stability, the M6b headline,
        # reproduced with no finite difference anywhere in the metric
        ge = hk_metric(strong, uw, θw; dZdu = dzdu(uw))
        gw = hk_metric((x, th) -> su2_torus(sw, x, th; R = Rw, chamber = :weak, tower = 4),
                       uw, θw; dZdu = dzdu(uw))
        @test maximum(abs, ge - gw) < 1e-8 * maximum(abs, ge)

        # and the metric is still the instanton-corrected one, not the semiflat limit
        gsf = semiflat_metric(strong(uw, θw), dzdu(uw))
        @test maximum(abs, ge - gsf) > 1e-3 * maximum(abs, gsf)

        # Ooguri-Vafa: Riemannian, and the corrections grow as the state gets light. Here
        # the twistor residual is bounded by the rapidity grid rather than by h, so the
        # exact route improves it only where the grid is comfortable (the lighter state).
        θov = [0.6, 1.1]
        f = (x, th) -> ooguri_vafa_torus(x, th; R = 1.0)
        gaps = Float64[]
        for z in (0.30 + 0.05im, 0.10 + 0.02im)
            dz = ooguri_vafa_period_derivatives(z)
            mp = metric_point(f, z, θov; dZdu = dz)
            d = hk_diagnostics(mp)
            @test d.residual < 1e-6
            @test all(>(0), d.eigenvalues)
            gsf = semiflat_metric(f(z, θov), dz)
            push!(gaps, maximum(abs, hk_metric(mp) - gsf) / maximum(abs, gsf))
        end
        @test gaps[2] > gaps[1]
        zl = 0.10 + 0.02im
        @test hk_diagnostics(metric_point(f, zl, θov;
                                          dZdu = ooguri_vafa_period_derivatives(zl))).residual <
              1e-2 * hk_diagnostics(metric_point(f, zl, θov)).residual
    end

    @testset "second derivatives of the metric" begin
        # ∂g by one difference of exact metrics, against one difference of finite-difference
        # metrics. **Recorded finding**: the second route is *not* noise-limited, contrary
        # to the expectation this stage was scoped on. A central difference of a smooth
        # deterministic solver carries a smooth O(h²) error, so differencing it again adds
        # no 1/h amplification - the iterated difference is truncation-limited in both
        # steps and keeps converging as the inner h shrinks. What the exact route buys on
        # second derivatives is therefore a factor, not a capability: one h to control
        # instead of two, at a ninth of the solves. The capability it does buy is the
        # testset above.
        gex(x) = hk_metric(strong, x, θw; dZdu = dzdu(x))
        scale = maximum(abs, gex(uw))
        A(H) = (gex(uw + H) - gex(uw - H)) / (2H)
        A1, A2 = A(1e-3), A(2e-3)
        @test maximum(abs, A1) > 1e-2 * scale           # the derivative is not small
        @test maximum(abs, A1 - A2) < 1e-5 * scale      # and it is stable in H

        B = (hk_metric(strong, uw + 1e-3, θw) - hk_metric(strong, uw - 1e-3, θw)) / 2e-3
        @test maximum(abs, A1 - B) < 1e-5 * scale
    end
end
