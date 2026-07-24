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
end
