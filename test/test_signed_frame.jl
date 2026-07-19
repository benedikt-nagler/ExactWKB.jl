# Signed-frame oracles. The central theorem: ε solves the keystone P = −εBε in EVERY
# chamber - including the cyclic ones (a triangle of the triangulation with all three
# edges internal), where the cocycle condition around the quiver 3-cycle is a genuine
# constraint and the uniform decay frame of M4 could not meet it. The reference (top)
# chamber is pinned to ε ≡ +1, and the combinatorial flip walk from it reproduces the
# traced triangulation in every chamber. ε is chamber-local: it is NOT transported,
# and NOT continuous across a wall even on a diagonal whose cycle persists.

import ClusterAlgebras

# chamber midpoints of a problem, on the wall circle [0, π)
function _chamber_midpoints(prob)
    ws = [w[1] for w in chamber_walls(prob)]
    [mod((ws[i] + (i == length(ws) ? ws[1] + π : ws[i + 1])) / 2, π)
     for i in eachindex(ws)]
end

_is_cyclic(t) = any(tri -> all(e -> t.is_diagonal[e], tri), t.triangles)

const CUBIC   = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
const QUARTIC = SchrodingerProblem([-1.0, 0.0, 0.0, 0.0, 1.0])
const QUINTIC = SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0])

@testset "signed frame" begin

    @testset "chamber_walls reduces onto the wall circle mod π" begin
        # the quintic's saddles at θ_c = 0 (pair (1,2)) and θ_c = π (pair (3,4)) are
        # the SAME wall - phases live mod π - and merge into one degenerate wall
        ws = chamber_walls(QUINTIC)
        @test all(0 ≤ w[1] < π for w in ws)
        @test issorted([w[1] for w in ws])
        @test ws[1][1] == 0
        @test Set(ws[1][2]) == Set([(1, 2), (3, 4)])
        @test length(ws) == 6
        # every confirmed saddle is accounted for exactly once
        @test sum(length(w[2]) for w in ws) == length(saddles(QUINTIC))
        # the quartic's two degenerate walls
        @test Set(w[2] for w in chamber_walls(QUARTIC) if length(w[2]) == 2) ==
              Set([[(1, 2), (3, 4)], [(1, 3), (2, 4)]])
    end

    @testset "reference chamber is the top chamber" begin
        for prob in (CUBIC, QUARTIC, QUINTIC)
            θref = reference_theta(prob)
            θmax = maximum(w[1] for w in chamber_walls(prob))
            @test θmax < θref < π
            # it is wall-free
            @test !any(θmax < w[1] < π for w in chamber_walls(prob))
        end
        @test reference_theta(QUINTIC) ≈ 3.0878 atol = 1e-3
    end

    # ── the central theorem ────────────────────────────────────────────────────────
    # ε solves P = −εBε in EVERY chamber. On a tree chamber the quiver graph has no
    # cycle and solvability is automatic; on a CYCLIC chamber (a triangle of the
    # triangulation with all three edges internal, generic for degree ≥ 5) the single
    # 3-cycle imposes a genuine cocycle condition - the one the uniform decay frame of
    # M4 could not meet. Four of the quintic's six chambers are cyclic.
    @testset "keystone solves in every chamber" begin
        n_cyclic = 0
        for (prob, name) in ((CUBIC, "cubic"), (QUARTIC, "quartic"), (QUINTIC, "quintic"))
            for θ in _chamber_midpoints(prob)
                t = ideal_triangulation(stokes_graph(prob; theta = θ))
                _is_cyclic(t) && (n_cyclic += 1)
                cb = charge_basis(prob, t)
                B = triangulation_quiver(t).B
                @test signed_pairing(cb) == -B
                @test all(in((-1, 1)), signs(cb))
                @test length(signs(cb)) == n_diagonals(t)
                # the signed seed is the finite type the triangulation demands
                @test ClusterAlgebras.is_finite_type(
                    ClusterAlgebras.Quiver(-signed_pairing(cb)))
            end
        end
        @test n_cyclic == 4              # the quintic's four cyclic chambers
    end

    @testset "reference-chamber pin: ε ≡ +1" begin
        for prob in (CUBIC, QUARTIC, QUINTIC)
            t = ideal_triangulation(stokes_graph(prob; theta = reference_theta(prob)))
            cb = charge_basis(prob, t)
            @test signs(cb) == ones(Int, n_diagonals(t))
            # equivalently: the raw decay pairing already IS the triangulation quiver
            @test cb.pairing == -triangulation_quiver(t).B
        end
    end

    # ── what connects the chambers ────────────────────────────────────────────────
    # The combinatorial flip walk from the reference chamber lands on exactly the
    # triangulation the tracer produces at the target - so every chamber's frame is
    # reached from a verified one by mutations, with no re-tracing in between.
    @testset "the flip walk reproduces the traced triangulation" begin
        for prob in (CUBIC, QUARTIC, QUINTIC)
            for θ in _chamber_midpoints(prob)
                r = verify_signed_frame(prob, θ)
                @test r.triangulation_matches
                @test r.walked == signed_frame(prob, θ)
            end
        end
        # the walk really does cross walls (it is not trivially the reference chamber)
        @test maximum(verify_signed_frame(QUINTIC, θ).n_walls
                      for θ in _chamber_midpoints(QUINTIC)) == 5
    end

    # ε is chamber-LOCAL: it is not transported, and it is not continuous across a
    # wall even on a diagonal whose turning-point pair (hence whose cycle) persists.
    # This is why `verify_signed_frame` asserts only the combinatorial claim.
    @testset "ε is chamber-local, not continuous across a wall" begin
        tref = ideal_triangulation(stokes_graph(CUBIC; theta = reference_theta(CUBIC)))
        tlow = ideal_triangulation(stokes_graph(CUBIC; theta = 0.7854))
        εref = signs(charge_basis(CUBIC, tref))
        εlow = signs(charge_basis(CUBIC, tlow))
        @test tref.diagonal_tp_pair == tlow.diagonal_tp_pair    # same cycles...
        @test εref != εlow                                      # ...different ε
    end

    @testset "gauge is inert" begin
        # negating ε globally leaves the seed, the masses and the wall phases fixed
        t = ideal_triangulation(stokes_graph(QUINTIC; theta = 0.5416))
        cb = charge_basis(QUINTIC, t)
        flipped = ChargeBasis{Float64}(cb.problem, cb.triangulation, cb.contours,
                                       cb.central_charges, cb.pairing, -cb.signs)
        @test signed_pairing(flipped) == signed_pairing(cb)
        @test abs.(physical_charges(flipped)) == abs.(physical_charges(cb))
        # wall phases agree as points on the circle of circumference π
        circdist(a, b) = (d = abs(mod(angle(a), π) - mod(angle(b), π));
                          min(d, π - d))
        @test all(circdist(a, b) < 1e-9
                  for (a, b) in zip(physical_charges(flipped), physical_charges(cb)))
    end

    @testset "rank 1 and the isolated-vertex convention" begin
        # harmonic: a single diagonal, no quiver edge - the keystone is vacuous there
        # and ε comes from the convention alone
        harm = SchrodingerProblem([-1.0, 0.0, 1.0])
        t = ideal_triangulation(stokes_graph(harm; theta = 1.0))
        @test n_diagonals(t) == 1
        @test signs(charge_basis(harm, t)) == [1]
    end

    @testset "error paths" begin
        @test_throws ContourError signed_frame([0 2; -2 0], ideal_triangulation(
            stokes_graph(CUBIC; theta = reference_theta(CUBIC))))
        # a sign cocycle that cannot close around an oriented 3-cycle: take P = −B and
        # flip one edge, so the edge signs multiply to −1 around the triangle
        B = [0 1 -1; -1 0 1; 1 -1 0]
        @test ExactWKB._solve_signs(-B, B) == [1, 1, 1]      # consistent
        P = -B
        P[1, 2], P[2, 1] = -P[1, 2], -P[2, 1]
        @test ExactWKB._solve_signs(P, B) === nothing        # cocycle fails
        # a magnitude mismatch is refused too
        @test ExactWKB._solve_signs(2 .* -B, B) === nothing
    end
end
