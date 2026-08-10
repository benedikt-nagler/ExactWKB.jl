# Charge-lattice oracles. The keystone consistency: the numerically measured
# intersection pairing of the decay-oriented cycles reproduces the combinatorial
# triangulation quiver up to cycle reorientation (|P| == |B|), and −P is a
# finite-type exchange matrix of the same Cartan type - with strict equality
# P == −B realized in the quintic's top chamber. The cubic values regress against
# the hand-built contours of test_ddp.jl, which pinned the conventions.

import ClusterAlgebras

@testset "charge lattice" begin
    cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
    tc = ideal_triangulation(stokes_graph(cubic; theta = 0.3))
    cbc = charge_basis(cubic, tc)

    @testset "cubic regression vs the hand-built DDP layer" begin
        @test n_charges(cbc) == 2
        tps = simple_turning_points(cubic)
        # the hand-built basis of test_ddp.jl: γ₁ = reversed (−1,0) cycle (Z₁ < 0),
        # γ₂ = the (0,1) cycle as built (Z₂ ∈ i𝐑₋)
        c1 = reverse(encircling_contour(tps[1], tps[2]; margin = 0.4))
        c2 = encircling_contour(tps[2], tps[3]; margin = 0.4)
        Z1 = period_integral(cubic, c1)
        Z2 = period_integral(cubic, c2)
        @test cbc.central_charges[1] ≈ Z1 rtol = 1e-8
        @test cbc.central_charges[2] ≈ Z2 rtol = 1e-8
        @test real(cbc.central_charges[1]) < 0
        @test imag(cbc.central_charges[2]) < 0
        κ = intersection_pairing(cubic, c2, c1)
        @test cbc.pairing == [0 -κ; κ 0]
        # ledger item 5 on the automatic basis: the seed's B is -P
        seed = bridge_seed(cbc)
        @test seed.quiver.B[1:2, 1:2] == -signed_pairing(cbc)
        @test ClusterAlgebras.cmatrix(seed) == [1 0; 0 1]
        @test ClusterAlgebras.cartan_type(ClusterAlgebras.Quiver(-cbc.pairing)) == (:A, 2)
    end

    @testset "keystone: |P| == |B|, -P finite type of the right Cartan type" begin
        for (coeffs, θ, ct) in [
            ([0.0, -1.0, 0.0, 1.0], 0.3, (:A, 2)),
            ([0.75, 0.0, -2.0, 0.0, 1.0], 0.25, (:A, 3)),
            ([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0], 3.09, (:A, 4)),   # quintic TOP chamber
        ]
            prob = SchrodingerProblem(coeffs)
            t = ideal_triangulation(stokes_graph(prob; theta = θ))
            cb = charge_basis(prob, t)
            B = triangulation_quiver(t).B
            @test abs.(cb.pairing) == abs.(B)
            @test cb.pairing == -transpose(cb.pairing)
            @test ClusterAlgebras.cartan_type(ClusterAlgebras.Quiver(-cb.pairing)) == ct
            # every basis cycle decays along its own wall (ledger item 4): each Voros
            # symbol built on it is a legal wall symbol. Numerically-real Z (wall at
            # θ_c = 0) is resolved by the real part - the same tolerance rule as the
            # production orientation.
            for Z in central_charges(cb)
                tol = sqrt(eps(Float64)) * (1 + abs(Z))
                @test imag(Z) < tol
                if abs(imag(Z)) ≤ tol
                    @test real(Z) < 0
                else
                    θc = mod(angle(Z), π)
                    @test real(Z * cis(-θc)) < 0
                end
            end
            # enclosure: winding 1 around the diagonal's pair, 0 elsewhere
            tps = simple_turning_points(prob)
            for (e, (i, j)) in enumerate(t.diagonal_tp_pair)
                w = [ExactWKB._winding(cb.contours[e], location(tp)) for tp in tps]
                @test abs(w[i]) == 1 && w[j] == w[i]
                @test all(iszero, w[k] for k in eachindex(tps) if k != i && k != j)
            end
        end
    end

    @testset "strict P == -B in the quintic top chamber (regression pin)" begin
        prob = SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0])
        t = ideal_triangulation(stokes_graph(prob; theta = 3.09))
        cb = charge_basis(prob, t)
        @test cb.pairing == -triangulation_quiver(t).B
    end

    @testset "cyclic chambers: the signed frame repairs the M4 finding" begin
        prob = SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0])
        t = ideal_triangulation(stokes_graph(prob; theta = 0.542))
        @test any(tri -> all(e -> t.is_diagonal[e], tri), t.triangles)
        # the M4 finding, kept as a record: the RAW decay frame is cluster-infinite
        # there, which is exactly why the uniform decay rule could not be a seed frame
        cb = charge_basis(prob, t)
        @test !ClusterAlgebras.is_finite_type(ClusterAlgebras.Quiver(-cb.pairing))
        # ... and the signed frame repairs it: the cocycle closes and the signed
        # pairing is the triangulation quiver on the nose
        @test signs(cb) == [1, 1, -1, -1] || signs(cb) == -[1, 1, -1, -1]
        @test signed_pairing(cb) == -triangulation_quiver(t).B
        @test ClusterAlgebras.is_finite_type(
            ClusterAlgebras.Quiver(-signed_pairing(cb)))
        @test ClusterAlgebras.cartan_type(
            ClusterAlgebras.Quiver(-signed_pairing(cb))) == (:A, 4)
    end

    @testset "precision from the caller (BigFloat)" begin
        setprecision(128) do
            probb = SchrodingerProblem([big(0.0), big(-1.0), big(0.0), big(1.0)])
            # topology in Float64 (a discrete datum), periods in BigFloat
            tb = ideal_triangulation(stokes_graph(probb; theta = 0.3))
            cbb = charge_basis(probb, tb; rtol = 1e-20)
            @test cbb.central_charges isa Vector{Complex{BigFloat}}
            @test cbb.pairing == cbc.pairing
            @test cbb.central_charges[1] ≈ cbc.central_charges[1] rtol = 1e-10
        end
    end

    @testset "the graph route agrees with the ellipse route on a disk" begin
        # `charge_basis(prob, g)` builds each cycle from its strip region instead of
        # an ellipse around a turning-point pair. On a disk the pair already names the
        # cycle, so the two routes must land on the same lattice: the same central
        # charges, and the same pairing and orientation signs *exactly*. The Z
        # tolerance is the price of the strip route - its core path runs out to a
        # marked point and back, so the two legs partly cancel.
        # θ = the midpoint of each fixture's widest wall-free gap (see
        # test_triangulation.jl; hard-coded here so this file stands alone)
        for (name, prob, θ) in (("cubic", cubic, 0.785),
                                ("quartic",
                                 SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0]), 0.893),
                                ("quintic",
                                 SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0]),
                                 0.542))
            g = stokes_graph(prob; theta = θ)
            t = ideal_triangulation(g)
            ce = charge_basis(prob, t)
            cg = charge_basis(prob, g)
            @test cg.triangulation.diagonal_tp_pair == t.diagonal_tp_pair
            @test n_charges(cg) == n_charges(ce)
            @test cg.pairing == ce.pairing
            @test signs(cg) == signs(ce)
            @test physical_charges(cg) ≈ physical_charges(ce) rtol = 1e-9
            # the keystone holds on both routes
            @test signed_pairing(cg) == -triangulation_quiver(t).B
        end
    end

    @testset "core paths and racetracks" begin
        g = stokes_graph(cubic; theta = 0.3)
        paths = diagonal_core_paths(g)
        @test length(paths) == n_diagonals(tc)
        tps = simple_turning_points(cubic)
        for (e, path) in enumerate(paths)
            i, j = tc.diagonal_tp_pair[e]
            # a core path runs from one turning point of its strip to the other
            @test path[1] ≈ location(tps[i])
            @test path[end] ≈ location(tps[j])
            @test length(path) ≥ 3
            # its racetrack encircles exactly those two
            c = charge_contour(cubic, path)
            for (k, tp) in enumerate(tps)
                w = ExactWKB._winding(c, location(tp))
                @test (k in (i, j)) ? abs(w) == 1 : w == 0
            end
            # and carries the same central charge as the ellipse around the pair
            @test abs(period_integral(cubic, c)) ≈
                  abs(period_integral(cubic, charge_contour(cubic, i, j))) rtol = 1e-9
        end
        @test_throws Resurgence.InvalidArgument charge_contour(cubic, [1.0 + 0im, 2.0 + 0im])
        # a half-width so large the third turning point is swallowed
        @test_throws ContourError charge_contour(cubic, paths[1]; margin = 5.0)
    end

    @testset "error paths" begin
        @test_throws Resurgence.InvalidArgument charge_contour(cubic, 1, 1)
        @test_throws Resurgence.InvalidArgument charge_contour(cubic, 0, 2)
        # a margin so large the third turning point cannot be excluded
        @test_throws ContourError charge_contour(cubic, 1, 2; margin = 5.0)
        # a doctored triangulation whose combinatorial B has empty support
        tbad = IdealTriangulation(tc.n_marked, tc.edge_endpoints, tc.is_diagonal,
                                  [(1, 3, 4), (2, 5, 6), (3, 6, 7)], tc.triangle_tp,
                                  tc.diagonal_tp_pair, tc.marked_angles, tc.theta)
        @test_throws ContourError charge_basis(cubic, tbad)
    end
end
