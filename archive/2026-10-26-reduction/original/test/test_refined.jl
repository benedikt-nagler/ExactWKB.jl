# The refined (motivic) bridge. The oracles pin the convention ledger of src/refined.jl:
# the skew form is the intersection pairing, the factor order is the sweep order, the
# pentagon is a statement about moduli (not about θ), and the refined jump collapses
# onto the classical DDP formula at v → 1.

@testset "refined bridge" begin
    strong = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])     # z³ − z: 2 BPS states
    weak   = SchrodingerProblem([-1.2, 0.0, 0.0, 1.0])     # z³ − 1.2: 3 BPS states
    sp_strong = bps_spectrum(strong; theta = reference_theta(strong))
    sp_weak   = bps_spectrum(weak; theta = 0.524)

    # ── ledger items 1 and 2: the torus and the ordering are the seed's ──────────────
    @testset "skew form and factor order agree with the bridge seed" begin
        cb = sp_strong.basis
        @test refined_skew(cb) == signed_pairing(cb)
        @test refined_skew(cb) == -bridge_seed(cb).quiver.B[1:2, 1:2]

        w = refined_dilog_word(sp_strong)
        wseed = ClusterAlgebras.quantum_dilog_word(bridge_seed(cb), sp_strong.sequence)
        @test w.charges == wseed.charges          # same charges, same order
        @test w.skew == wseed.skew                # same quantum torus
        @test all(==(1), w.weights)               # simply laced: no q^{d_k} weighting
        @test refined_dt(sp_strong; truncation_degree = 4) ==
              ClusterAlgebras.ks_dilog_product(wseed; truncation_degree = 4)
    end

    # ── the headline oracle: the quantum pentagon, from two potentials ───────────────
    # Across the wall of marginal stability the bound state γ₁+γ₂ appears; the refined
    # DT invariant does not change. This is the A₂ pentagon identity reached from
    # physics rather than from two green sequences of one seed.
    @testset "quantum pentagon across the wall of marginal stability" begin
        @test n_states(sp_strong) == 2
        @test n_states(sp_weak) == 3
        @test charges(sp_weak)[2] == [1, 1]                # the bound state, in the middle

        @test_throws ChamberError refined_frame_map(sp_strong, sp_weak)  # different potentials
        r = verify_refined_wall_crossing(sp_strong, sp_weak; truncation_degree = 5)
        @test r.equal
        @test r.map == [1 0; 0 1]
        @test r.n_states == (2, 3)
        @test !isempty(r.product)

        # the identity is not vacuous: the quantum torus is doing work, i.e. the
        # same word over a commutative torus (Λ = 0) gives a different answer
        commuting = ClusterAlgebras.QuantumDilogWord(charges(sp_weak),
                                                     zeros(Int, 2, 2),
                                                     copy(sp_weak.sequence))
        @test ClusterAlgebras.ks_dilog_product(commuting; truncation_degree = 5) !=
              r.product

        # ledger item 1: with the opposite skew form the two sides disagree
        flipped = ClusterAlgebras.QuantumDilogWord(charges(sp_weak),
                                                   -refined_skew(sp_weak.basis),
                                                   copy(sp_weak.sequence))
        @test ClusterAlgebras.ks_dilog_product(flipped; truncation_degree = 5) !=
              refined_dt(sp_strong; truncation_degree = 5)
    end

    # ── ledger item 3: rotating θ is NOT an invariance ───────────────────────────────
    # Both θ-chambers of z³ − z carry the same two states, in opposite sweep order and
    # with ε flipped on the second basis cycle. Rotating the half-plane past a BPS ray
    # replaces γ by −γ, so the ordered product changes - a wall of the second kind.
    @testset "θ-rotation changes the product (and must)" begin
        sp_other = bps_spectrum(strong; theta = 0.3)
        @test n_states(sp_other) == n_states(sp_strong) == 2
        @test signs(sp_other.basis) != signs(sp_strong.basis)
        M = refined_frame_map(sp_strong, sp_other)
        @test M == [1 0; 0 -1]                    # the second cycle reverses
        rotated = ClusterAlgebras.QuantumDilogWord([M * c for c in charges(sp_other)],
                                                   refined_skew(sp_strong.basis),
                                                   copy(sp_other.sequence))
        @test ClusterAlgebras.ks_dilog_product(rotated; truncation_degree = 4) !=
              refined_dt(sp_strong; truncation_degree = 4)

        # and the trap the ledger warns about: the permutation that matches the skew
        # forms is a basis SWAP, under which the products agree for a reason that has
        # nothing to do with the physical cycles
        swap = refined_quiver_iso(sp_strong, sp_other)
        @test swap == [0 1; 1 0]
        swapped = ClusterAlgebras.QuantumDilogWord([swap * c for c in charges(sp_other)],
                                                    refined_skew(sp_strong.basis),
                                                    copy(sp_other.sequence))
        @test ClusterAlgebras.ks_dilog_product(swapped; truncation_degree = 4) ==
              refined_dt(sp_strong; truncation_degree = 4)
    end

    # ── the refined DDP jump and its classical limit ─────────────────────────────────
    @testset "refined DDP jump collapses onto the classical formula" begin
        cb = sp_strong.basis
        P = signed_pairing(cb)
        for k in 1:n_charges(cb)
            v = verify_refined_ddp(cb, k; truncation_degree = 4)
            @test v.matches
            @test v.exponents == P[:, k]
        end
        # a cycle is continuous across its own wall: ⟨γ_k, γ_k⟩ = 0
        self = ClusterAlgebras.ks_classical_limit(refined_ddp_jump(cb, 1, 1;
                                                                   truncation_degree = 3))
        @test self == Dict([1, 0] => 1 // 1)
        @test_throws Resurgence.InvalidArgument refined_ddp_jump(cb, 1, 3)
    end

    # ── higher rank: the physical sweep is the quiver's refined DT invariant ─────────
    # A₃ and A₄ chambers have no second potential to cross to here, so the invariance
    # tested is Reineke/Keller's: every maximal green sequence of the chamber's own
    # seed gives the same product as the physically swept one.
    @testset "A₃ / A₄ chambers: every green sequence gives the same product" begin
        for (prob, rank) in ((SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0]), 3),
                             (SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0]), 4))
            sp = bps_spectrum(prob)
            @test n_charges(sp.basis) == rank
            seed = bridge_seed(sp.basis)
            mine = refined_dt(sp; truncation_degree = 3)
            others = ClusterAlgebras.maximal_green_sequences(seed)
            @test sp.sequence in others
            for seq in others
                w = ClusterAlgebras.quantum_dilog_word(seed, seq)
                @test ClusterAlgebras.ks_dilog_product(w; truncation_degree = 3) == mine
            end
        end
    end

    @testset "errors" begin
        sp3 = bps_spectrum(SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0]))
        @test_throws ChamberError refined_quiver_iso(sp_strong, sp3)  # ranks differ

        # a state with Ω ≠ 1 is never absorbed silently
        vect = ExactWKB.BPSSpectrum([ExactWKB.BPSState([1, 0], 1.0 + 0.0im, -2)],
                                    [1], sp_strong.basis, 1.0)
        @test_throws Resurgence.InvalidArgument refined_dilog_word(vect)
    end
end
