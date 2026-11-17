# Tests for the SU(2) BPS spectrum + Kronecker wall-crossing (`src/sw_bps.jl`).
# Oracles: both-chamber spectrum content + Ω, the affine-Kronecker cluster capability check
# against ClusterAlgebras (Ã₁ recognition + the dyon-tower c-vectors), the self-validating KS
# wall-crossing identity, and the tower's phase accumulation onto the W-boson.

const CA = ExactWKB.ClusterAlgebras

@testset "sw_bps" begin
    @testset "spectrum content & Ω" begin
        sw = SeibergWittenSU2()
        u = 25.0
        strong = su2_bps_states(sw, u; chamber = :strong)
        @test length(strong) == 2
        @test Set(ExactWKB.charge.(strong)) == Set([[1, 0], [-1, 2]])   # monopole, dyon
        @test all(ExactWKB.omega.(strong) .== 1)

        weak = su2_bps_states(sw, u; chamber = :weak, tower = 3)
        # exactly one W-boson, with Ω = −2, charge (0,2)
        wbosons = filter(s -> ExactWKB.omega(s) == -2, weak)
        @test length(wbosons) == 1
        @test ExactWKB.charge(wbosons[1]) == [0, 2]
        # everything else is a hypermultiplet (Ω = 1)
        @test all(s -> ExactWKB.omega(s) == 1, filter(s -> ExactWKB.omega(s) != -2, weak))
        # δ = γ₁ + γ₂
        @test [0, 2] == [1, 0] + [-1, 2]
        # masses positive, central-charge additivity Z_δ = Z_{γ₁} + Z_{γ₂}
        @test all(s -> ExactWKB.mass(s) > 0, weak)
        @test central_charge(sw, u, (0, 2)) ≈ central_charge(sw, u, (1, 0)) +
                                              central_charge(sw, u, (-1, 2))
    end

    @testset "affine-Kronecker cluster capability check (ClusterAlgebras)" begin
        q = su2_bps_quiver()
        @test CA.is_affine_type(q)            # Ã₁
        @test !CA.is_finite_type(q)           # cluster-infinite
        @test CA.is_mutation_finite(q)        # B-class {±B}: mutation-finite but not finite-type
        # mutating the principal Kronecker seed alternately at 1,2,1,2,… reproduces the
        # preprojective tower roots (n+1, n)
        s = CA.extend(CA.Seed(q))
        cur = s
        roots = Vector{Int}[]
        for k in [1, 2, 1, 2, 1]
            push!(roots, CA.c_vector(cur, k))
            cur = CA.mutate(cur, k)
        end
        @test roots == [[1, 0], [2, 1], [3, 2], [4, 3], [5, 4]]
    end

    @testset "KS wall-crossing identity (m = 2, pentagon analog)" begin
        for d in (2, 3, 4, 5)
            res = verify_su2_wall_crossing(; degree = d)
            @test res.closed
            @test res.max_residual == 0
        end
        @test_throws Resurgence.InvalidArgument verify_su2_wall_crossing(; degree = 0)
    end

    @testset "auto-chamber (wall seam)" begin
        sw = SeibergWittenSU2()
        # inside the wall the default resolves to the 2-state strong spectrum …
        @test sw_chamber(sw, 0.2im) == :strong
        strong = su2_bps_states(sw, 0.2im)
        @test length(strong) == 2
        @test Set(ExactWKB.charge.(strong)) == Set([[1, 0], [-1, 2]])
        # … with genuinely strong-coupling central charges: nonzero, distinct phases
        Zs = [central_charge(sw, 0.2im, (1, 0)), central_charge(sw, 0.2im, (-1, 2))]
        @test all(abs.(Zs) .> 1e-2)
        @test abs(angle(Zs[1]) - angle(Zs[2])) > 1e-2
        # outside the wall the default resolves to the weak tower
        weak = su2_bps_states(sw, 10.0; tower = 4)
        @test length(weak) == 2 * 5 + 1
        # an explicit chamber is still honored anywhere
        @test length(su2_bps_states(sw, 10.0; chamber = :strong)) == 2
    end

    @testset "the NS-deformed spectrum" begin
        sw = SeibergWittenSU2()
        u = 1.9 + 0.8im

        # ħ = 0 at order 1 is the classical spectrum exactly, state for state
        for ch in (:strong, :weak)
            s0 = su2_bps_states(sw, u; chamber = ch)
            s1 = su2_bps_states(sw, u; chamber = ch, ħ = 0.0, order = 1)
            @test length(s0) == length(s1)
            @test all(charge(a) == charge(b) for (a, b) in zip(s0, s1))
            @test all(ExactWKB.omega(a) == ExactWKB.omega(b) for (a, b) in zip(s0, s1))
            @test all(central_charge(a) == central_charge(b) for (a, b) in zip(s0, s1))
        end

        # the deformation moves Z but not the BPS quiver: charges and Ω are untouched
        # (ledger item 3), and the mass shift scales as ħ²
        s0 = su2_bps_states(sw, u; chamber = :weak, tower = 4)
        shift(h) = maximum(abs(central_charge(a) - central_charge(b)) for (a, b) in
                           zip(s0, su2_bps_states(sw, u; chamber = :weak, tower = 4,
                                                  ħ = h, order = 1)))
        for h in (0.1, 0.2)
            sh = su2_bps_states(sw, u; chamber = :weak, tower = 4, ħ = h, order = 1)
            @test all(charge(a) == charge(b) for (a, b) in zip(s0, sh))
            @test all(ExactWKB.omega(a) == ExactWKB.omega(b) for (a, b) in zip(s0, sh))
        end
        @test isapprox(log(shift(0.2) / shift(0.1)) / log(2), 2; atol = 0.05)

        # the deformed Z is n_m a_D + n_e a from the *same* deformed pair - which is what
        # makes the twistor layer's basis/state cross-check pass by construction
        qs = quantum_sw_periods(sw, u; order = 1)
        h = 0.3
        aD = Resurgence.evaluate(qs.a_D, h)
        a = Resurgence.evaluate(qs.a, h)
        for s in su2_bps_states(sw, u; chamber = :weak, tower = 2, ħ = h, order = 1)
            nm, ne = charge(s)
            @test central_charge(s) ≈ nm * aD + ne * a
        end

        # :auto resolves against the deformed wall, which at ħ = 0 is the classical one
        auto0 = su2_bps_states(sw, u)
        auto1 = su2_bps_states(sw, u; ħ = 0.0, order = 1)
        @test length(auto0) == length(auto1)
        @test all(charge(a) == charge(b) && central_charge(a) == central_charge(b)
                  for (a, b) in zip(auto0, auto1))
    end

    @testset "phase accumulation onto the W-boson (geometry tie)" begin
        sw = SeibergWittenSU2()
        u = 25.0
        aphase = angle(central_charge(sw, u, (0, 2)))   # arg Z_W = arg(a) (= 0, a real)
        # preprojective dyons (1, 2n): phase → arg(a) from above (converges like 1/n)
        pp = [angle(central_charge(sw, u, (1, 2n))) for n in 1:20]
        @test abs(pp[end] - aphase) < abs(pp[1] - aphase)
        @test abs(pp[end] - aphase) < 5e-2
        # preinjective dyons (−1, 2n+2): phase → arg(a) from below
        pi_ = [angle(central_charge(sw, u, (-1, 2n + 2))) for n in 1:20]
        @test abs(pi_[end] - aphase) < abs(pi_[1] - aphase)
        @test abs(pi_[end] - aphase) < 5e-2
    end
    # ── the refined (motivic) identity ───────────────────────────────────────────
    # The affine case the finite-type refined layer cannot reach. Self-validating in
    # the same sense as the classical one: closure alone pins the W-boson's refined
    # index AND its quadratic-refinement sign - nothing is taken from the literature.
    @testset "refined SU(2) wall-crossing" begin
        for d in 2:6
            r = verify_su2_refined_wall_crossing(; degree = d)
            @test r.closed
            @test r.delta_only                       # the residue lives on δ alone
            @test r.residue == r.vector_factor       # and it IS the vector multiplet
        end
        @test_throws Resurgence.InvalidArgument verify_su2_refined_wall_crossing(; degree = 0)

        # ledger R2: every other candidate for the W-boson factor fails, and already
        # at the first order in ŷ^δ - the sign is not decorative
        D = 4
        rays(fs) = ExactWKB._su2_refined_product(
            [((1, 1), fs)], D)
        measured = verify_su2_refined_wall_crossing(; degree = D).residue
        for wrong in ([(1, 1, -1), (-1, 1, -1)],      # no quadratic-refinement sign
                      [(0, -1, -2)],                  # unshifted, squared
                      [(0, 1, -1)],                   # unshifted, Ω = −1
                      [(2, -1, -1), (-2, -1, -1)])    # shifted by v^{±2}
            @test rays(wrong) != measured
        end
        # the two sub-factors commute (δ is isotropic), so their order is no convention
        @test rays([(1, -1, -1), (-1, -1, -1)]) == rays([(-1, -1, -1), (1, -1, -1)])

        # ledger R3: at y → 1 the refined index is the classical Ω(δ) = −2
        weak = su2_refined_rays(:weak, 3)
        δ = only(r for r in weak if r[1] == (1, 1))
        @test sum(f[3] for f in δ[2]) == -2
        @test all(sum(f[3] for f in r[2]) == 1 for r in weak if r[1] != (1, 1))
        # and the ray content matches the classical list state for state
        @test [r[1] for r in weak] == [r[1] for r in ExactWKB._su2_rays(:weak, 3)] ||
              Set(r[1] for r in weak) == Set(r[1] for r in ExactWKB._su2_rays(:weak, 3))
        @test [r[1] for r in su2_refined_rays(:strong, 3)] == [(0, 1), (1, 0)]
        @test_throws Resurgence.InvalidArgument su2_refined_rays(:nonsense, 2)
        @test_throws Resurgence.InvalidArgument su2_refined_rays(:weak, -1)
    end

end
