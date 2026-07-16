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

    @testset "cubic chamber pin (ledger item 6)" begin
        for θ in (nothing, 0.3)                  # default chamber and the DDP chamber
            sp = θ === nothing ? bps_spectrum(cubic) : bps_spectrum(cubic; theta = θ)
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

    @testset "quintic (A₄ stress test): physics in the tree chamber" begin
        prob = SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0])
        sp = bps_spectrum(prob)                  # default finds a tree chamber
        sads = saddles(prob)
        @test n_states(sp) == length(sads)
        @test _matches_saddles(sp, sads)
        @test all(all(≥(0), c) && any(>(0), c) for c in charges(sp))
        @test all(ExactWKB.omega(s) == 1 for s in sp.states)
        seed = bridge_seed(sp.basis)
        # the chamber seed is a non-bipartite A₄ orientation; Zamolodchikov
        # periodicity h + 2 = 7 is checked on the bipartite representative of the
        # type recognized from OUR seed
        @test ClusterAlgebras.cartan_type(ClusterAlgebras.Quiver(-sp.basis.pairing)) ==
              (:A, 4)
        Bbip = [0 1 0 0; -1 0 -1 0; 0 1 0 1; 0 0 -1 0]   # alternating A₄ orientation
        @test ClusterAlgebras.cartan_type(ClusterAlgebras.Quiver(Bbip)) == (:A, 4)
        seed_bip = ClusterAlgebras.extend(ClusterAlgebras.Seed(ClusterAlgebras.Quiver(Bbip)))
        @test ClusterAlgebras.y_system(seed_bip).period == 7
        dt = ClusterAlgebras.dt_transformation(seed)
        Cf = ClusterAlgebras.cmatrix(dt.seed)
        @test all(Cf[dt.sigma[j], j] == -1 for j in 1:4)
    end

    @testset "chamber independence (cubic, two chambers)" begin
        spa = bps_spectrum(cubic; theta = 0.3)
        spb = bps_spectrum(cubic; theta = 2.0)
        @test n_states(spa) == n_states(spb)
        ka = sort([(mass(s), phase(s)) for s in spa.states])
        kb = sort([(mass(s), phase(s)) for s in spb.states])
        @test all(isapprox(a[1], b[1]; rtol = 1e-6) && _phase_dist(a[2], b[2]) < 1e-6
                  for (a, b) in zip(ka, kb))
    end

    @testset "loop closure: the automatic basis feeds the DDP layer" begin
        # the flagship consistency: bps_spectrum's chamber data, handed to the
        # hand-pinned DDP machinery, verifies wall-crossing = y-mutation with no
        # hand-built input (margin 0.4 as in test_ddp for clean high-order periods)
        sp = bps_spectrum(cubic; theta = 0.3, margin = 0.4)
        cb = sp.basis
        w = wkb_expansion(cubic; order = 12)
        vs = [voros_symbol(w, c) for c in cb.contours]
        seed = bridge_seed(cb)
        k = findfirst(Z -> abs(imag(Z)) < 1e-8, central_charges(cb))  # the θ_c = 0 wall
        @test k !== nothing
        ħ = abs(central_charges(cb)[k]) / 5
        res = verify_ddp_mutation(vs, seed, k, ħ; theta = 0.0)
        @test res.max_residual < 3e-2
    end

    @testset "errors" begin
        quintic = SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0])
        @test_throws ChamberError bps_spectrum(quintic; theta = 0.542)  # cyclic chamber
        @test_throws NonGenericGraph bps_spectrum(cubic; theta = 0.0)   # on a wall
        # M4 finding part 2: the quintic fan chamber is a tree chamber whose all-decay
        # sweep is cluster-consistent but unphysical (6 states vs 7 saddles) — the
        # saddle gate catches it
        @test_throws ChamberError bps_spectrum(quintic; theta = 0.052)
    end
end
