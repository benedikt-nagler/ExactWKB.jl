# BPS-spectrum oracles: the cluster-algebraic enumeration must reproduce physics.
# State count and (mass, phase) multiset match the confirmed saddles; the cubic
# chamber order is the ledger-item-6 pin [[0,1],[1,0]]; Ω ≡ 1 in finite type;
# Zamolodchikov periodicity h + 2 and DT closure hold on the derived seeds; and the
# spectrum is chamber-independent (θ only reorders the sweep, never the physics).

import ClusterAlgebras

# circular distance of two wall phases (period π)
_phase_dist(a, b) = min(mod(a - b, π), mod(b - a, π))

# every state matches a confirmed saddle in (mass, phase), bijectively
function _matches_saddles(sp, sads; mrtol = 1e-5, ptol = 1e-5)
    length(sp.states) == length(sads) || return false
    used = falses(length(sads))
    for s in sp.states
        i = findfirst(k -> !used[k] &&
                          abs(mass(sads[k]) - mass(s)) ≤ mrtol * mass(sads[k]) &&
                          _phase_dist(ExactWKB.theta(sads[k]), phase(s)) ≤ ptol,
                      eachindex(sads))
        i === nothing && return false
        used[i] = true
    end
    true
end

@testset "bps" begin
    cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])

    @testset "harmonic: a single state" begin
        harm = SchrodingerProblem([-1.0, 0.0, 1.0])
        sp = bps_spectrum(harm)
        @test n_states(sp) == 1
        @test charges(sp) == [[1]]
        @test ExactWKB.omega(sp.states[1]) == 1
        @test phase(sp.states[1]) ≈ π / 2 atol = 1e-8
        @test mass(sp.states[1]) ≈ mass(only(saddles(harm))) rtol = 1e-6
        @test sp.sequence == [1]
    end

    @testset "Airy: the empty spectrum" begin
        sp = bps_spectrum(SchrodingerProblem([0.0, 1.0]))
        @test n_states(sp) == 0
        @test isempty(sp.sequence)
        @test n_charges(sp.basis) == 0
    end

    # Ledger item 6 is pinned in the REFERENCE chamber - the gauge in which it was
    # originally derived (ε ≡ +1, and the chamber-relative order degenerates to the
    # absolute θ-decreasing one as θ₀ → π⁻). The old values survive there unchanged.
    @testset "cubic chamber pin (ledger item 6, reference chamber)" begin
        sp = bps_spectrum(cubic; theta = reference_theta(cubic))
        @test signs(sp.basis) == [1, 1]
        @test n_states(sp) == 2
        @test charges(sp) == [[0, 1], [1, 0]]
        @test phase(sp.states[1]) ≈ π / 2 atol = 1e-8
        @test phase(sp.states[2]) ≈ 0 atol = 1e-8
        @test all(ExactWKB.omega(s) == 1 for s in sp.states)
        @test _matches_saddles(sp, saddles(cubic))
        # the greedy sequence is one of the honest maximal green sequences
        seed = bridge_seed(sp.basis)
        @test sp.sequence in ClusterAlgebras.maximal_green_sequences(seed)
        @test ClusterAlgebras.ordered_c_vectors(seed, sp.sequence) == charges(sp)
    end

    # Away from the reference chamber the SWEEP ORDER differs - that is the whole
    # content of the chamber-relative rule - while the physics does not. At θ₀ = 0.3
    # the wall at θ_c = 0 is the nearest one below, so it is crossed first; in the
    # reference chamber it is crossed last.
    @testset "cubic: order is chamber-relative, physics is not" begin
        for θ in (nothing, 0.3)
            sp = θ === nothing ? bps_spectrum(cubic) : bps_spectrum(cubic; theta = θ)
            @test n_states(sp) == 2
            @test charges(sp) == [[1, 0], [0, 1]]        # the reversed sweep order
            @test phase(sp.states[1]) ≈ 0 atol = 1e-8
            @test phase(sp.states[2]) ≈ π / 2 atol = 1e-8
            @test _matches_saddles(sp, saddles(cubic))   # same physics
            seed = bridge_seed(sp.basis)
            @test sp.sequence in ClusterAlgebras.maximal_green_sequences(seed)
            @test ClusterAlgebras.ordered_c_vectors(seed, sp.sequence) == charges(sp)
        end
    end

    @testset "quartic double well: commuting tie + physics ($(θ))" for θ in (nothing,)
        prob = SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0])
        sp = bps_spectrum(prob)
        sads = saddles(prob)
        @test n_states(sp) == length(sads)
        @test _matches_saddles(sp, sads)
        @test all(all(≥(0), c) && any(>(0), c) for c in charges(sp))
        @test all(ExactWKB.omega(s) == 1 for s in sp.states)
        seed = bridge_seed(sp.basis)
        @test ClusterAlgebras.y_system(seed).period == 6          # A₃: h + 2
        dt = ClusterAlgebras.dt_transformation(seed)
        Cf = ClusterAlgebras.cmatrix(dt.seed)
        @test all(Cf[dt.sigma[j], j] == -1 for j in 1:3)
    end

    @testset "quintic (A₄ stress test): physics in every chamber" begin
        prob = SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0])
        sp = bps_spectrum(prob)
        sads = saddles(prob)
        @test n_states(sp) == length(sads)
        @test _matches_saddles(sp, sads)
        @test all(all(≥(0), c) && any(>(0), c) for c in charges(sp))
        @test all(ExactWKB.omega(s) == 1 for s in sp.states)
        seed = bridge_seed(sp.basis)
        # the chamber seed is a non-bipartite A₄ orientation; Zamolodchikov
        # periodicity h + 2 = 7 is checked on the bipartite representative of the
        # type recognized from OUR seed
        @test ClusterAlgebras.cartan_type(
            ClusterAlgebras.Quiver(-signed_pairing(sp.basis))) == (:A, 4)
        Bbip = [0 1 0 0; -1 0 -1 0; 0 1 0 1; 0 0 -1 0]   # alternating A₄ orientation
        @test ClusterAlgebras.cartan_type(ClusterAlgebras.Quiver(Bbip)) == (:A, 4)
        seed_bip = ClusterAlgebras.extend(ClusterAlgebras.Seed(ClusterAlgebras.Quiver(Bbip)))
        @test ClusterAlgebras.y_system(seed_bip).period == 7
        dt = ClusterAlgebras.dt_transformation(seed)
        Cf = ClusterAlgebras.cmatrix(dt.seed)
        @test all(Cf[dt.sigma[j], j] == -1 for j in 1:4)
    end

    # The flagship of the signed-flip layer: EVERY chamber gives the same physics.
    # Under M4 this could only be asserted for two cubic chambers - the quintic's four
    # cyclic chambers were refused outright and its fan chamber swept 6 states.
    @testset "chamber independence (every chamber, every fixture)" begin
        quartic4 = SchrodingerProblem([-1.0, 0.0, 0.0, 0.0, 1.0])
        quintic5 = SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0])
        for prob in (cubic, quartic4, quintic5)
            sads = saddles(prob)
            ws = [w[1] for w in chamber_walls(prob)]
            mids = [mod((ws[i] + (i == length(ws) ? ws[1] + π : ws[i + 1])) / 2, π)
                    for i in eachindex(ws)]
            @test length(mids) ≥ 2
            for θ in mids
                sp = bps_spectrum(prob; theta = θ)
                @test n_states(sp) == length(sads)
                @test _matches_saddles(sp, sads)
            end
            # ... and the (mass, phase) multiset is literally the same in each
            keys_of(θ) = sort([(round(mass(s), digits = 8), round(phase(s), digits = 8))
                               for s in bps_spectrum(prob; theta = θ).states])
            ref = keys_of(first(mids))
            @test all(keys_of(θ) == ref for θ in mids)
        end
    end

    @testset "loop closure: the automatic basis feeds the DDP layer" begin
        # the flagship consistency: bps_spectrum's chamber data, handed to the
        # hand-pinned DDP machinery, verifies wall-crossing = y-mutation with no
        # hand-built input (margin 0.4 as in test_ddp for clean high-order periods).
        #
        # Run in the REFERENCE chamber, where ε ≡ +1: vertex j of the signed seed
        # carries ε_j·γ_j, so `voros_symbol(w, cb.contours[j])` models y_j only where
        # ε_j == +1 - and `_require_decaying` would (correctly) refuse a reversed wall
        # symbol. The assertion below makes that precondition explicit rather than
        # relying on it holding by luck.
        sp = bps_spectrum(cubic; theta = reference_theta(cubic), margin = 0.4)
        cb = sp.basis
        @test signs(cb) == ones(Int, n_charges(cb))
        w = wkb_expansion(cubic; order = 12)
        vs = [voros_symbol(w, c) for c in cb.contours]
        seed = bridge_seed(cb)
        k = findfirst(Z -> abs(imag(Z)) < 1e-8, central_charges(cb))  # the θ_c = 0 wall
        @test k !== nothing
        @test signs(cb)[k] == 1
        ħ = abs(central_charges(cb)[k]) / 5
        res = verify_ddp_mutation(vs, seed, k, ħ; theta = 0.0)
        @test res.max_residual < 3e-2
    end

    # The two chambers the M4 layer could not handle, as named regressions.
    @testset "the chambers M4 could not do" begin
        quintic = SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0])
        sads = saddles(quintic)
        @test length(sads) == 7

        # (1) a CYCLIC chamber - charge_basis refused it outright under M4
        t = ideal_triangulation(stokes_graph(quintic; theta = 0.542))
        @test any(tri -> all(e -> t.is_diagonal[e], tri), t.triangles)
        sp = bps_spectrum(quintic; theta = 0.542)
        @test n_states(sp) == 7
        @test _matches_saddles(sp, sads)
        @test signs(sp.basis) == [1, 1, -1, -1] || signs(sp.basis) == -[1, 1, -1, -1]

        # (2) the real-axis fan chamber - a TREE chamber whose all-decay sweep was
        # cluster-consistent but unphysical (6 states against 7 saddles)
        sp = bps_spectrum(quintic; theta = 0.052)
        @test n_states(sp) == 7
        @test _matches_saddles(sp, sads)
    end

    @testset "errors" begin
        @test_throws NonGenericGraph bps_spectrum(cubic; theta = 0.0)   # on a wall
    end
end
