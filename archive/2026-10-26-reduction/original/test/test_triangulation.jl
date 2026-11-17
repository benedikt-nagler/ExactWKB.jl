# Triangulation oracles: a generic degree-d Stokes graph triangulates a (d+2)-gon
# (Iwaki–Nakanishi), so the quiver type must be A_{d−1} - Airy → rank 0, harmonic → A₁,
# cubic → A₂, quartic double well → A₃, quintic → A₄ - with all combinatorial
# invariants (marked points d+2, triangles d, diagonals d−1, regions 2d+1) exact, the
# triangulation constant on chambers, and a wall crossing acting as one flip = one
# quiver mutation.

import ClusterAlgebras

# midpoint of the largest wall-free gap of the candidate phases (period π)
function _generic_theta(prob)
    ths = sort(unique([ExactWKB.theta(s) for s in saddle_candidates(prob)]))
    isempty(ths) && return 0.5
    gaps = [i == length(ths) ? ths[1] + π - ths[end] : ths[i + 1] - ths[i]
            for i in eachindex(ths)]
    i = argmax(gaps)
    mod(ths[i] + gaps[i] / 2, float(π))
end

@testset "triangulation" begin
    fixtures = [
        ("Airy", SchrodingerProblem([0.0, 1.0]), 1, nothing),
        ("harmonic", SchrodingerProblem([-1.0, 0.0, 1.0]), 2, (:A, 1)),
        ("cubic", SchrodingerProblem([0.0, -1.0, 0.0, 1.0]), 3, (:A, 2)),
        ("quartic double well", SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0]), 4,
         (:A, 3)),
        # z⁵ − z − 1/2: an asymmetric quintic (z⁵ − z puts turning points exactly on
        # the saddle paths of other pairs)
        ("quintic", SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0]), 5, (:A, 4)),
    ]

    @testset "combinatorial invariants + cluster type ($name)" for (name, prob, n, ct) in fixtures
        θ = _generic_theta(prob)
        g = stokes_graph(prob; theta = θ)
        @test length(infinite_lines(g)) == 3n
        @test length(ray_exit_angles(g)) == 3n
        t = ideal_triangulation(g)
        @test n_marked_points(t) == n + 2
        @test length(t.triangles) == n
        @test t.triangle_tp == collect(1:n)
        @test n_diagonals(t) == n - 1
        @test diagonals(t) == collect(1:(n - 1))
        @test length(t.edge_endpoints) == 2n + 1
        @test t.is_diagonal == vcat(trues(n - 1), falses(n + 2))
        @test length(t.marked_angles) == n + 2
        @test issorted(t.marked_angles)
        @test all(0 ≤ a < 2π for a in t.marked_angles)
        # each triangle has three distinct edges; every diagonal borders exactly two
        # triangles, every boundary edge exactly one
        @test all(allunique, t.triangles)
        usage = zeros(Int, 2n + 1)
        for tri in t.triangles, e in tri
            usage[e] += 1
        end
        @test usage == vcat(fill(2, n - 1), fill(1, n + 2))
        # the boundary edges close up into the (n+2)-gon
        bdry = sort(t.edge_endpoints[n:end])
        @test bdry == sort(vcat([(i, i + 1) for i in 1:(n + 1)], [(1, n + 2)]))
        # diagonals join non-adjacent marked points
        for e in 1:(n - 1)
            (u, v) = t.edge_endpoints[e]
            @test u != v
            @test !((v - u == 1) || (u == 1 && v == n + 2))
        end
        # diagonal ↔ turning-point-pair bookkeeping is canonical
        @test issorted(t.diagonal_tp_pair)
        @test allunique(t.diagonal_tp_pair)
        q = triangulation_quiver(t)
        @test q.n_mutable == n - 1
        @test q.B == -transpose(q.B)
        if ct !== nothing
            @test ClusterAlgebras.cartan_type(q) == ct
        end
        @test q.labels == ["γ($(p[1]),$(p[2]))" for p in t.diagonal_tp_pair]
    end

    @testset "chamber invariance (cubic)" begin
        prob = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
        ta = ideal_triangulation(stokes_graph(prob; theta = 0.2))
        tb = ideal_triangulation(stokes_graph(prob; theta = 0.3))
        @test ta.diagonal_tp_pair == tb.diagonal_tp_pair
        @test ta.edge_endpoints == tb.edge_endpoints
        @test triangulation_quiver(ta).B == triangulation_quiver(tb).B
    end

    # Crossing the wall of the saddle (i, j) flips the diagonal dual to the collapsing
    # strip, which reconnects the SAME turning-point pair on the other side (the two
    # adjacent triangles are the same two turning points before and after) - so the
    # crossed diagonal keeps its label, and the flip is the quiver mutation there.
    #
    # The labelling does NOT in general match vertex by vertex across the wall: two of
    # the four quad sides swap which turning point borders them, so their labels change
    # and the canonical order (diagonals sorted by `diagonal_tp_pair`) permutes. The
    # rank-2 cubic case below cannot see this - at rank 2, μ₁ and μ₂ give the same B.
    # `canonical_reorder(flip(t, k)) == t′` in test_signed_frame.jl is the oracle that
    # pins the permutation, and a degenerate wall (several saddles at one θ_c) flips
    # several commuting diagonals at once.
    @testset "flip = mutation across the θ_c = 0 wall (cubic)" begin
        prob = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
        tp = ideal_triangulation(stokes_graph(prob; theta = 0.15))
        tm = ideal_triangulation(stokes_graph(prob; theta = -0.15))
        @test tp.diagonal_tp_pair == tm.diagonal_tp_pair
        k = findfirst(==((1, 2)), tp.diagonal_tp_pair)   # the crossed wall's saddle
        @test k !== nothing
        Bp, Bm = triangulation_quiver(tp).B, triangulation_quiver(tm).B
        @test Bp != Bm
        @test ClusterAlgebras.mutate(triangulation_quiver(tp), k).B == Bm
    end

    @testset "flip = mutation at rank 3 (quartic double well, wall (2,3))" begin
        prob = SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0])
        # the inner pair (2,3) has real central charge: its wall sits at θ_c = 0
        @test any(s -> ExactWKB.pair(s) == (2, 3) && ExactWKB.theta(s) < 1e-8,
                  saddle_candidates(prob))
        tp = ideal_triangulation(stokes_graph(prob; theta = 0.1))
        tm = ideal_triangulation(stokes_graph(prob; theta = -0.1))
        @test tp.diagonal_tp_pair == tm.diagonal_tp_pair
        k = findfirst(==((2, 3)), tp.diagonal_tp_pair)
        @test k !== nothing
        Bp, Bm = triangulation_quiver(tp).B, triangulation_quiver(tm).B
        @test Bp != Bm
        @test ClusterAlgebras.mutate(triangulation_quiver(tp), k).B == Bm
        # mutating at any other vertex does NOT reach the other chamber
        for j in 1:3
            j == k && continue
            @test ClusterAlgebras.mutate(triangulation_quiver(tp), j).B != Bm
        end
    end

    @testset "combinatorial flip" begin
        prob = SchrodingerProblem([-0.5, -1.0, 0.0, 0.0, 0.0, 1.0])   # quintic, A₄
        t = ideal_triangulation(stokes_graph(prob; theta = 3.09))

        # flipping back the way you came is the identity, and the flip is index-aligned
        # with the quiver mutation. Re-crossing in the SAME sense is not the identity:
        # it exchanges the two quad triangles' turning points (see `flip`'s docstring).
        for k in diagonals(t)
            f = flip(t, k)
            back = flip(f, k; direction = -1)
            @test triangulation_quiver(f).B ==
                  ClusterAlgebras.mutate(triangulation_quiver(t), k).B
            @test back.edge_endpoints == t.edge_endpoints
            @test back.triangles == t.triangles
            @test back.diagonal_tp_pair == t.diagonal_tp_pair
            @test flip(f, k).triangles != t.triangles
            # only the flipped diagonal moves; it keeps its own turning-point pair
            @test f.edge_endpoints[setdiff(1:end, k)] ==
                  t.edge_endpoints[setdiff(1:end, k)]
            @test f.diagonal_tp_pair[k] == t.diagonal_tp_pair[k]
            @test f.is_diagonal == t.is_diagonal
            @test f.theta == t.theta          # provenance is carried, not recomputed
            # every edge is still used by exactly two triangles (diagonal) or one
            for e in eachindex(f.edge_endpoints)
                @test count(tri -> e in tri, f.triangles) == (f.is_diagonal[e] ? 2 : 1)
            end
        end

        # canonical_reorder is the identity on a freshly traced triangulation
        tc, perm = canonical_reorder(t)
        @test perm == 1:length(t.edge_endpoints)
        @test tc.edge_endpoints == t.edge_endpoints
        @test tc.triangles == t.triangles

        @test_throws Resurgence.InvalidArgument flip(t, 0)
        @test_throws Resurgence.InvalidArgument flip(t, n_diagonals(t) + 1)  # boundary
        @test_throws Resurgence.InvalidArgument flip(t, 1; direction = 0)
    end

    @testset "surface invariants: a polynomial problem gives a disk" begin
        # The general construction reduces to the disk when Q has no poles, and the
        # divisor identity n = M + 2b − 4 and the region count 2n + 2 − b reduce to
        # d = (d+2) + 2 − 4 and 2d + 1. Corners are aligned with edges: edge slot i
        # runs from corner i to corner i+1 - the datum `flip` needs on any surface.
        for (name, prob, d, _) in fixtures
            t = ideal_triangulation(stokes_graph(prob; theta = _generic_theta(prob)))
            @test ExactWKB.n_boundaries(t) == 1
            @test all(==(1), t.marked_boundary)
            @test n_marked_points(t) == d + 2
            @test length(t.triangles) == length(t.triangle_corners)
            for s in eachindex(t.triangles)
                for i in 1:3
                    e = t.triangles[s][i]
                    c1 = t.triangle_corners[s][i]
                    c2 = t.triangle_corners[s][mod1(i + 1, 3)]
                    @test t.edge_endpoints[e] == (min(c1, c2), max(c1, c2))
                end
            end
        end
    end

    @testset "non-generic graphs are refused" begin
        prob = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])
        # θ exactly on the (1,2) wall: the finite line is present
        gwall = stokes_graph(prob; theta = 0.0, allow_incomplete = true)
        @test_throws NonGenericGraph ideal_triangulation(gwall)
        err = try
            ideal_triangulation(gwall)
        catch e
            e
        end
        @test err isa NonGenericGraph
        @test err isa ExactWKBError
        @test occursin("wall", sprint(showerror, err))
        # incomplete rays are refused too
        ginc = stokes_graph(prob; theta = 0.3, max_mass = 0.5, allow_incomplete = true)
        @test any(l -> ExactWKB.endpoint(l) === :incomplete, ExactWKB.lines(ginc))
        @test_throws NonGenericGraph ideal_triangulation(ginc)
    end
end
