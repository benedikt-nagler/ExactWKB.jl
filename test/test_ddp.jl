# DDP-layer oracles. The fixture is the cubic Q = z³ − z (type A₂): turning points
# −1, 0, 1; the (−1,0) cycle has real central charge (wall at θ_c = 0, the test wall),
# the (0,1) cycle imaginary Z (wall at θ_c = π/2); the composite (1,3) candidate is not
# a confirmed saddle - this chamber holds two BPS states, the minimal A₂ chamber.
#
# Airy (Q = z) cannot pin the jump - it has no closed even-enclosure cycle, so no
# Voros symbol; its normalization role is carried by Resurgence's Euler/Airy pins,
# inherited here through lateral_sum. The ExactWKB-side normalization is pinned by the
# truncating harmonic oscillator instead (first testset).
#
# Contour parameters: WKB order 12 with margin = 0.4 keeps the even-period roundoff
# check clean in Float64 (the theorem itself holds to 1e-63 in BigFloat).

using Resurgence: evaluate, borel, pade, poles
import ClusterAlgebras

@testset "ddp" begin
    # shared cubic fixture
    prob = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
    tps = simple_turning_points(prob)
    w = wkb_expansion(prob; order = 12)
    # γ₁: the (−1,0) cycle, reversed so Z₁ < 0 (decaying along its wall θ_c = 0);
    # γ₂: the (0,1) cycle, decaying along its own wall θ_c = π/2 as built.
    c1 = reverse(encircling_contour(tps[1], tps[2]; margin = 0.4))
    c2 = encircling_contour(tps[2], tps[3]; margin = 0.4)
    v1 = voros_symbol(w, c1)
    v2 = voros_symbol(w, c2)
    Zwall = abs(classical_period(v1))            # |Z_{γ₁}| ≈ 0.9585
    hs = [Zwall / x for x in (5.0, 6.5, 8.0)]

    @testset "voros_value normalization (harmonic, truncating WKB)" begin
        harm = SchrodingerProblem([0.0, 0.0, 1.0]; energy = 1.0)
        wh = wkb_expansion(harm; order = 4)
        vs = voros_symbol(wh, encircling_contour(turning_points(harm)...))
        for ħ in (0.1, 0.35), side in (:plus, :minus, :median)
            direct = exp(evaluate(full_series(vs), ħ))
            @test voros_value(vs, ħ; side) ≈ direct rtol = 1e-8
        end
    end

    @testset "Borel walls sit on the Z-lattice" begin
        # the Borel singularity of V_γ₂ on the θ_c = 0 ray sits at ζ = −Z_{γ₁}
        # (so that the jump carries e^{-ζ/ħ} = e^{Z_{γ₁}/ħ} = the decaying V_{γ₁})
        B2, r2 = ExactWKB._voros_approximant(v2)
        ps = poles(r2)
        ζ = -classical_period(v1)
        nearest = ps[argmin(abs.(ps .- ζ))]
        @test abs(nearest - ζ) < 0.02 * abs(ζ)
        @test abs(mod(angle(nearest), π)) < 1e-4 ||
              abs(mod(angle(nearest), π) - π) < 1e-4     # arg ≡ θ_c = 0 (mod π)
    end

    @testset "the DDP jump: sign, normalization, integer Stokes constant" begin
        κ = intersection_pairing(prob, c2, c1)
        @test abs(κ) == 1                                # primitive charge, Ω(γ') = 1
        for ħ in hs
            r = verify_ddp(v2, v1, κ, ħ; theta = 0.0)
            @test r.relative_residual < 3e-2
            @test round(Int, real(r.kappa_measured)) == κ
            @test abs(r.kappa_measured - κ) < 0.05
            # the wrong pairing sign fails by far more than the residual
            rwrong = verify_ddp(v2, v1, -κ, ħ; theta = 0.0)
            @test rwrong.relative_residual > 10 * r.relative_residual
            # ddp_transform maps the minus-side value to the plus-side one - and the
            # opposite direction is visibly wrong (off by (1+V')^{2κ})
            Vp_plus = voros_value(v2, ħ; theta = 0.0, side = :plus)
            Vp_minus = voros_value(v2, ħ; theta = 0.0, side = :minus)
            right = abs(ddp_transform(Vp_minus, r.v_prime, κ) / Vp_plus - 1)
            wrong = abs(ddp_transform(Vp_plus, r.v_prime, κ) / Vp_minus - 1)
            @test right < 1e-3
            @test wrong > 10 * right
        end
    end

    @testset "pairing properties" begin
        κ = intersection_pairing(prob, c2, c1)
        @test intersection_pairing(prob, c1, c2) == -κ            # antisymmetry
        @test intersection_pairing(prob, reverse(c2), c1) == -κ   # orientation flip
        # nested homologous copy: same class, no crossings
        c1in = reverse(encircling_contour(tps[1], tps[2]; margin = 0.2))
        @test intersection_pairing(prob, c1in, c2) == -κ          # homologous to c1
        @test intersection_pairing(prob, c1in, c1) == 0
        # disjoint cycles of the double well pair to zero
        dw = SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0])
        dtps = sort(turning_points(dw); by = t -> real(location(t)))
        d1 = encircling_contour(dtps[1], dtps[2]; margin = 0.15)
        d2 = encircling_contour(dtps[3], dtps[4]; margin = 0.15)
        @test intersection_pairing(dw, d1, d2) == 0
        # a contour vertex on a turning point is refused
        bad = [0.5 - 0.5im, 1.0 + 0.0im, 1.5 + 0.5im, 1.0 + 1.0im]
        @test_throws ContourError intersection_pairing(prob, bad, c1)
    end

    @testset "jump = y-mutation (direction and B-sign pin)" begin
        κ = intersection_pairing(prob, c2, c1)
        P = [0 -κ; κ 0]                                  # P[i,j] = ⟨γ_i, γ_j⟩
        seed = ddp_seed(P)
        @test ClusterAlgebras.cartan_type(ClusterAlgebras.Quiver(P)) == (:A, 2)
        @test ClusterAlgebras.c_vector(seed, 1) == [1, 0]          # green before μ₁
        ħ = hs[1]
        res = verify_ddp_mutation([v1, v2], seed, 1, ħ; theta = 0.0)
        @test res.max_residual < 3e-2
        # wrong B-sign (equivalently wrong mutation direction) fails by ≥ 10×
        reswrong = verify_ddp_mutation([v1, v2], ddp_seed(-P), 1, ħ; theta = 0.0)
        @test reswrong.max_residual > 10 * res.max_residual
        # wrong vertex: γ₂'s wall is not at θ = 0 - the decay guard refuses it
        @test_throws Resurgence.InvalidArgument verify_ddp_mutation(
            [v1, v2], seed, 2, ħ; theta = 0.0)
        # tropical pin: the jump's leading (tropical) part is c-vector mutation
        ms = ClusterAlgebras.mutate(seed, 1)
        Bmat = seed.quiver.B
        @test ClusterAlgebras.c_vector(ms, 1) == [-1, 0]
        @test ClusterAlgebras.c_vector(ms, 2) == [max(Bmat[1, 2], 0), 1]
    end

    @testset "pentagon / DT closure on the WKB chamber data" begin
        κ = intersection_pairing(prob, c2, c1)
        seed = ddp_seed([0 -κ; κ 0])
        # confirmed saddles = the two-state chamber, θ-ordered γ₁ (θ_c=0), γ₂ (θ_c=π/2)
        sads = sort(saddles(prob); by = s -> ExactWKB.theta(s))
        @test [s.pair for s in sads] == [(1, 2), (2, 3)]
        # charge additivity where it is exact (open-path orientations fix the sign):
        # Z₁₃ = Z₁₂ − Z₂₃ for saddle_candidates' straight-path orientations
        cands = saddle_candidates(prob)
        Z = Dict(s.pair => central_charge(s) for s in cands)
        @test Z[(1, 3)] ≈ Z[(1, 2)] - Z[(2, 3)] rtol = 1e-6
        # chamber ↔ MGS identification (ledger item 6): the MGS charge order is the
        # θ-DECREASING order of the confirmed saddle phases - here γ₂ (π/2), γ₁ (0)
        mgs = ClusterAlgebras.maximal_green_sequences(seed)
        @test any(seq -> ClusterAlgebras.ordered_c_vectors(seed, seq) ==
                         [[0, 1], [1, 0]], mgs)
        @test !any(seq -> ClusterAlgebras.ordered_c_vectors(seed, seq) ==
                          [[1, 0], [0, 1]], mgs)          # θ-increasing does NOT match
        # Zamolodchikov periodicity: h + 2 = 5 (the pentagon)
        @test ClusterAlgebras.y_system(seed).period == 5
        # DT transformation: Ω(γ) ≡ 1, tropical y_images = y ↦ y⁻¹ (C_final = −P_σ)
        dt = ClusterAlgebras.dt_transformation(seed)
        word = ClusterAlgebras.quantum_dilog_word(seed, dt.sequence)
        @test all(==(1), values(ClusterAlgebras.omega(word)))
        Cf = ClusterAlgebras.cmatrix(dt.seed)
        @test all(Cf[dt.sigma[j], j] == -1 for j in 1:2)
    end

    @testset "precision from the caller (BigFloat)" begin
        # 128 bits and order 10 keep the pin while the BigFloat quadrature (whose
        # rtol scales with eps(BigFloat)) stays affordable
        setprecision(128) do
            probb = SchrodingerProblem([big(0.0), big(-1.0), big(0.0), big(1.0)])
            tpsb = simple_turning_points(probb)
            wb = wkb_expansion(probb; order = 10)
            c1b = reverse(encircling_contour(tpsb[1], tpsb[2]; margin = big(0.4)))
            c2b = encircling_contour(tpsb[2], tpsb[3]; margin = big(0.4))
            v1b = voros_symbol(wb, c1b)
            v2b = voros_symbol(wb, c2b)
            κ = intersection_pairing(probb, c2b, c1b)
            ħ = abs(classical_period(v1b)) / 5
            r = verify_ddp(v2b, v1b, κ, ħ; theta = 0.0, rtol = 1e-20)
            @test r.jump_measured isa Complex{BigFloat}
            @test r.relative_residual < 3e-2             # Padé-order-limited
            @test abs(r.kappa_measured - κ) < 0.05
        end
    end

    @testset "error paths" begin
        ħ = hs[1]
        # growing wall symbol (un-reversed γ₁ cycle) is refused, never auto-flipped
        v1grow = voros_symbol(w, reverse(c1))
        @test_throws Resurgence.InvalidArgument verify_ddp(v2, v1grow, 1, ħ; theta = 0.0)
        # pairing 0: predicted jump 0, absolute residual, no NaN
        r0 = verify_ddp(v2, v1, 0, ħ; theta = 0.0)
        @test iszero(r0.jump_predicted)
        @test isfinite(r0.relative_residual)
        # invalid side / invalid pairing matrices / bad mutation calls
        @test_throws Resurgence.InvalidArgument voros_value(v1, ħ; side = :sideways)
        @test_throws Resurgence.InvalidArgument ddp_seed([0 1; 1 0])
        @test_throws Resurgence.InvalidArgument ddp_seed([0 1 0; -1 0 0])
        seed = ddp_seed([0 -1; 1 0])
        @test_throws Resurgence.InvalidArgument verify_ddp_mutation([v1], seed, 1, ħ)
        @test_throws Resurgence.InvalidArgument verify_ddp_mutation([v1, v2], seed, 3, ħ)
        @test_throws Resurgence.InvalidArgument verify_ddp_mutation(
            [v1, v2], ClusterAlgebras.mutate(seed, 1), 1, ħ)
    end
end
