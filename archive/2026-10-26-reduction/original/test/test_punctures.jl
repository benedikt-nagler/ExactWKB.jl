using ClusterAlgebras

# Regular punctures: a DOUBLE pole of Q. Unlike an irregular pole it is not a boundary
# circle but a single interior vertex of the dual triangulation, so a once-punctured
# m-gon is cluster type D_m - the first non-A, non-affine types the bridge reaches.
#
# The fixtures are `Q = N(z)/z²` with `N(0) = 1` (so the residue √c = 1 and the
# ring-domain wall sits at θ = π/2) and `deg N = m`, which puts m marked points at
# infinity and m simple turning points: the once-punctured m-gon.
_punctured(m_coeffs) = RationalProblem(m_coeffs, [0.0], [2])
const PUNC3 = _punctured([1.0, -0.3, 0.0, 1.0])              # once-punctured triangle
const PUNC4 = _punctured([1.0, 0.0, 0.0, 0.0, 1.0])          # once-punctured square
const PUNC5 = _punctured([1.0, 0.2, 0.0, 0.0, 0.0, 1.0])     # once-punctured pentagon

# Moduli whose chamber structure reaches valence 1, i.e. a SELF-FOLDED triangle. The
# same once-punctured m-gons, so the cluster type is still D_m and the self-folded
# chamber has to reproduce it - which it only does through the FST π-substitution.
const SELF3 = _punctured([0.2, 1.0, 0.0, 1.0])               # valence 1 for θ ∈ (1.19, 1.95)
const SELF4 = _punctured([0.2, 1.0, 0.0, 0.0, 1.0])
const SELF5 = _punctured([0.2, 1.0, 0.0, 0.0, 0.0, 1.0])

@testset "regular punctures (double poles)" begin
    @testset "a double pole is admitted, a simple pole is not" begin
        @test PUNC3 isa RationalProblem
        @test pole_orders(PUNC3) == [2]
        @test puncture_indices(PUNC3) == [1]
        @test n_punctures(PUNC3) == 1
        @test n_punctures(mathieu_problem(1.0, 3.0)) == 0    # two irregular poles
        @test n_punctures(SchrodingerProblem([0, 1])) == 0
        @test puncture_indices(SchrodingerProblem([0, 1])) == Int[]

        # a simple pole is an orbifold point: still refused
        @test_throws InvalidPotential RationalProblem([1.0, 1.0], [0.0], [1])
        # infinity must stay irregular: Q = 1/z² is a puncture at ∞ too
        @test_throws InvalidPotential RationalProblem([1.0], [0.0], [2])

        # mixed orders in one problem
        mixed = RationalProblem([1.0, 0.0, 0.0, 0.0, 0.0, 1.0], [0.0, 1.0], [2, 3])
        @test n_punctures(mixed) == 1
        @test puncture_indices(mixed) == [1]
        @test ExactWKB.degree(mixed) == 0
    end

    @testset "a puncture carries no marked points" begin
        # |e + 2| = 0 at a double pole: the empty vector, NOT an error - the Stokes
        # lines spiral in and approach no direction at all.
        @test asymptotic_directions(PUNC3, 0.4; pole = 1) == Float64[]
        @test length(asymptotic_directions(PUNC3, 0.4; pole = 0)) == 3
        # the divisor identity Σ(zero orders) = M + 2b − 4 with M counting only
        # boundary marked points: 3 = 3 + 2·2 − 4
        @test length(turning_points(PUNC3)) == 3
        @test length(turning_points(PUNC5)) == 5
    end

    @testset "the ring-domain wall of a puncture" begin
        w = only(ring_domain_walls(PUNC3))
        @test w isa RingDomainWall
        @test pole_index(w) == 1
        # Z_p = ∮√Q = 2πi√c with c = N(0) = 1
        @test central_charge(w) ≈ 2π * im
        @test mass(w) ≈ 2π
        @test ExactWKB.theta(w) ≈ π / 2
        @test isempty(ring_domain_walls(mathieu_problem(1.0, 3.0)))
        @test isempty(ring_domain_walls(SchrodingerProblem([0, 1])))

        # the residue sets the wall: scaling N(0) by i rotates it by π/4
        w2 = only(ring_domain_walls(_punctured([im, -0.3 + 0im, 0.0im, 1.0 + 0im])))
        @test ExactWKB.theta(w2) ≈ mod(π / 2 + π / 4, π)

        # ON the wall the trajectories close into a ring domain instead of spiralling,
        # and the winding diverges: the tracer refuses rather than integrating it
        @test_throws ExactWKB.TracingFailed stokes_graph(PUNC3; theta = π / 2)
        # ... and says so with a raisable budget
        @test_throws ExactWKB.TracingFailed stokes_graph(PUNC3; theta = 1.55,
                                                         pole_radius = 0.001)
        @test stokes_graph(PUNC3; theta = 1.55, pole_radius = 0.001,
                           max_winding = 1e6) isa StokesGraph
    end

    @testset "rays spiral into the puncture in finite mass" begin
        g = stokes_graph(PUNC3; theta = 0.4)
        into = filter(l -> endpoint(l) === :pole, ExactWKB.lines(g))
        @test length(into) == 3                       # one per turning point
        @test all(l -> ExactWKB.target(l) == 1, into)
        @test length(ExactWKB.infinite_lines(g)) == 6
        @test length(boundary_lines(g)) == 9          # every ray reaches a singularity

        # the mass to the pole circle is |√c|·log(R/r)/|cos(α−θ)| - LOGARITHMIC, which
        # is the branch the ζ^{1−m/2} budget of an irregular pole cannot express. The
        # bound is generous (the ray does not start at the escape radius), so this is
        # an order-of-magnitude oracle on the branch, not on the constant.
        r = 0.95 / 20                                  # the default pole_radius here
        @test all(l -> mass(l) < 1.05 * log(5 * 1.9 / r) / abs(cos(0.4)), into)
        @test all(l -> mass(l) > 0.5 * log(1 / r) / abs(cos(0.4)), into)
        # each ray ends ON the circle it was stopped at
        @test all(l -> 0.9r < abs(ExactWKB.points(l)[end]) < 1.1r, into)
    end

    @testset "the once-punctured m-gon is cluster type D_m" begin
        # The oracle of the whole layer. D_3 = A_3 so it is only a check of the count
        # there; D_4 and D_5 are types no polynomial potential can reach.
        for (prob, m, θ) in ((PUNC3, 3, 0.4), (PUNC4, 4, 0.3), (PUNC5, 5, 0.3))
            t = ideal_triangulation(stokes_graph(prob; theta = θ))
            @test n_punctures(t) == 1
            @test n_boundaries(t) == 2                 # infinity and the puncture
            @test n_marked_points(t) == m + 1          # m on the circle, plus the
            @test count(t.marked_is_puncture) == 1     # puncture itself
            @test t.marked_boundary == vcat(ones(Int, m), [2])
            @test t.marked_is_puncture == vcat(falses(m), [true])
            # diagonals = n − 2 + nb = m, which is FST's 6g + 3b + 3p + |M| − 6
            @test n_diagonals(t) == m
            @test count(!, t.is_diagonal) == m         # boundary edges = marked points
            @test length(t.triangles) == m
            q = triangulation_quiver(t)
            @test size(q.B, 1) == m
            @test string(ClusterAlgebras.mutation_type(q)) ==
                  (m == 3 ? "A3" : "D$m")              # D_3 = A_3
        end
    end

    @testset "the puncture is a vertex, and its valence is the chamber" begin
        # θ = 0.4: every turning point sends one ray in, so the triangulation is the
        # m radii from the puncture and the puncture has valence m.
        t = ideal_triangulation(stokes_graph(PUNC3; theta = 0.4))
        p = findfirst(t.marked_is_puncture)
        @test p == 4
        @test puncture_valence(t, p) == 3
        @test sort(t.edge_endpoints[t.is_diagonal]) == [(1, 4), (2, 4), (3, 4)]
        @test sort(t.edge_endpoints[.!t.is_diagonal]) == [(1, 2), (1, 3), (2, 3)]

        # θ = 1.5 is a different chamber (across the saddle wall at θ ≈ 1.431): only
        # two rays reach the puncture, so it has valence 2 and one strip region reverts
        # to two boundary arcs. Both region types in one graph - this is the case the
        # "ends" classification exists for.
        t2 = ideal_triangulation(stokes_graph(PUNC3; theta = 1.5))
        p2 = findfirst(t2.marked_is_puncture)
        @test puncture_valence(t2, p2) == 2
        @test n_diagonals(t2) == 3                     # the count does not move
        @test string(ClusterAlgebras.mutation_type(triangulation_quiver(t2))) == "A3"
        # one diagonal now joins two BOUNDARY marked points, going round the puncture:
        # it is a different arc from the boundary edge with the same endpoints
        @test count(e -> !(p2 in e), t2.edge_endpoints[t2.is_diagonal]) == 1
        @test length(unique(t2.edge_endpoints)) < length(t2.edge_endpoints)
    end

    @testset "flip is mutation on a punctured surface" begin
        for (prob, θ) in ((PUNC3, 0.4), (PUNC4, 0.3), (PUNC5, 0.3))
            t = ideal_triangulation(stokes_graph(prob; theta = θ))
            q = triangulation_quiver(t)
            for k in 1:n_diagonals(t)
                @test triangulation_quiver(flip(t, k)).B ==
                      ClusterAlgebras.mutate(q, k).B
                @test flip(flip(t, k), k; direction = -1).edge_endpoints ==
                      t.edge_endpoints
            end
        end
    end

    @testset "what the layer refuses, it refuses loudly" begin
        t = ideal_triangulation(stokes_graph(PUNC3; theta = 0.4))
        # The charge lattice is NOT built on a punctured surface: a cycle would have to
        # be told apart from the residue cycle around the puncture (horizon.md, Tier B).
        @test_throws ExactWKB.ContourError charge_basis(PUNC3,
                                                        stokes_graph(PUNC3; theta = 0.4))
        @test_throws ExactWKB.ContourError ExactWKB._require_unpunctured(t)
        # ... and it still runs on the unpunctured surfaces it was built for
        @test ExactWKB._require_unpunctured(
            ideal_triangulation(stokes_graph(mathieu_problem(1.0, 3.0);
                                             theta = 0.3))) === nothing

    end
end

# -- self-folded triangles: a puncture of valence 1 -----------------------------------
#
# A single Stokes ray spiralling into the puncture makes the face walk send it straight
# back out along the previous ray of its OWN turning point, so the region it bounds is
# a strip with both corners at that one turning point: FST's radius. Its dual triangle
# carries the radius twice and the enclosing loop once, and the exchange matrix needs
# the π-substitution π(radius) = loop or the once-punctured m-gon stops being D_m.
#
# (problem, θ in a valence-1 chamber, m)
const SELF_CASES = ((SELF3, 1.4, 3), (SELF4, 1.06, 4), (SELF5, 1.38, 5))

@testset "self-folded triangles (puncture valence 1)" begin
    @testset "the dual triangle carries the same edge twice" begin
        for (prob, θ, m) in SELF_CASES
            t = ideal_triangulation(stokes_graph(prob; theta = θ))
            p = findfirst(t.marked_is_puncture)
            @test puncture_valence(t, p) == 1
            s = only(selffolded_triangles(t))
            @test is_selffolded(t, s)
            @test count(x -> is_selffolded(t, x), eachindex(t.triangles)) == 1
            r, ℓ = selffold_arcs(t, s)
            # the radius appears twice in its triangle, the loop once
            @test count(==(r), t.triangles[s]) == 2
            @test count(==(ℓ), t.triangles[s]) == 1
            # the radius joins a boundary marked point to the puncture; the loop runs
            # from that marked point back to itself
            @test p in t.edge_endpoints[r]
            M = only(x for x in t.edge_endpoints[r] if x != p)
            @test t.edge_endpoints[ℓ] == (M, M)
            # the triangle's corners are (M, P, M) up to rotation
            @test sort(collect(t.triangle_corners[s])) == sort([M, M, p])
            # the radius borders ONE triangle, on both sides, which is what the
            # turning-point pair records
            @test t.diagonal_tp_pair[r] == (t.triangle_tp[s], t.triangle_tp[s])
            # nothing else in the layer moved: the counts of M9 all survive
            @test n_diagonals(t) == m
            @test n_marked_points(t) == m + 1
            @test length(t.triangles) == m
        end
    end

    @testset "the FST π-substitution is what keeps the type D_m" begin
        for (prob, θ, m) in SELF_CASES
            t = ideal_triangulation(stokes_graph(prob; theta = θ))
            r, ℓ = selffold_arcs(t, only(selffolded_triangles(t)))
            B = triangulation_quiver(t).B
            @test B == -transpose(B)
            # THE rule: no arrow between the radius and its loop, and the radius'
            # row is a copy of the loop's everywhere else ([FST08] Def. 4.1)
            @test B[r, ℓ] == 0
            @test all(B[r, j] == B[ℓ, j] for j in 1:m if j != r && j != ℓ)
            # and therefore the chamber still has the type of the surface. D_3 = A_3,
            # but D_4 and D_5 are the real oracles - dropping π gives A_2 × A_1 etc.
            @test string(ClusterAlgebras.mutation_type(triangulation_quiver(t))) ==
                  (m == 3 ? "A3" : "D$m")
        end
    end

    @testset "flip is mutation, and the radius has no ideal flip" begin
        for (prob, θ, m) in SELF_CASES
            t = ideal_triangulation(stokes_graph(prob; theta = θ))
            q = triangulation_quiver(t)
            r, ℓ = selffold_arcs(t, only(selffolded_triangles(t)))
            # the radius borders its one triangle on both sides, so there is no quad
            # to flip it across; the error names the loop, whose flip undoes the fold
            err = try
                flip(t, r)
                nothing
            catch e
                e
            end
            @test err isa ExactWKB.NonGenericGraph
            @test occursin("radius", sprint(showerror, err))
            @test occursin("flip the enclosing loop $ℓ", sprint(showerror, err))
            # every other diagonal flips, and flipping the loop removes the self-fold
            for k in 1:n_diagonals(t)
                k == r && continue
                tf = flip(t, k)
                @test triangulation_quiver(tf).B == ClusterAlgebras.mutate(q, k).B
                @test flip(tf, k; direction = -1).edge_endpoints == t.edge_endpoints
                @test flip(tf, k; direction = -1).triangles == t.triangles
                @test flip(tf, k; direction = -1).triangle_corners == t.triangle_corners
            end
            @test isempty(selffolded_triangles(flip(t, ℓ)))
        end
    end

    @testset "a flip creates the self-fold, and it is the wall crossing" begin
        # SELF3 has the valence-2 chamber below θ ≈ 1.19 and the valence-1 chamber
        # above it. Crossing the wall by hand reproduces the traced graph exactly -
        # not just the quiver, but the triangles, corners and turning-point pairs.
        t1 = ideal_triangulation(stokes_graph(SELF3; theta = 1.15))
        t2 = ideal_triangulation(stokes_graph(SELF3; theta = 1.25))
        @test isempty(selffolded_triangles(t1))
        @test length(selffolded_triangles(t2)) == 1
        k = only(k for k in 1:n_diagonals(t1)
                 if first(canonical_reorder(flip(t1, k; direction = 1))).triangles ==
                    t2.triangles)
        tf, perm = canonical_reorder(flip(t1, k; direction = 1))
        @test perm == collect(1:length(t1.edge_endpoints))
        @test tf.edge_endpoints == t2.edge_endpoints
        @test tf.triangle_corners == t2.triangle_corners
        @test tf.diagonal_tp_pair == t2.diagonal_tp_pair
        # the mirror crossing at the far end of the valence-1 chamber, θ ≈ 1.95
        t3 = ideal_triangulation(stokes_graph(SELF3; theta = 1.90))
        t4 = ideal_triangulation(stokes_graph(SELF3; theta = 2.05))
        @test length(selffolded_triangles(t3)) == 1
        @test isempty(selffolded_triangles(t4))
        ℓ = selffold_arcs(t3, only(selffolded_triangles(t3)))[2]
        tg, _ = canonical_reorder(flip(t3, ℓ; direction = 1))
        @test tg.edge_endpoints == t4.edge_endpoints
        @test tg.triangles == t4.triangles
        @test tg.triangle_corners == t4.triangle_corners
    end

    @testset "the tagged triangulation names the two arcs apart" begin
        for (prob, θ, m) in SELF_CASES
            t = ideal_triangulation(stokes_graph(prob; theta = θ))
            p = findfirst(t.marked_is_puncture)
            r, ℓ = selffold_arcs(t, only(selffolded_triangles(t)))
            ta = tagged_arcs(t)
            @test length(ta) == m
            # τ(t) replaces the loop by the radius NOTCHED at the puncture, so the two
            # now share endpoints and are told apart only by the tag - which is why the
            # exchange matrix has no arrow between them
            @test ta[r] == TaggedArc(t.edge_endpoints[r], (:plain, :plain))
            @test ta[ℓ].endpoints == t.edge_endpoints[r]
            notched = only(i for i in 1:2 if ta[ℓ].tags[i] == :notched)
            @test ta[ℓ].endpoints[notched] == p
            @test length(unique(ta)) == m
            # every other arc is plain at both ends
            @test all(a.tags == (:plain, :plain) for (k, a) in enumerate(ta) if k != ℓ)
        end
        # with no self-fold the tagged form is the diagonals themselves, all plain
        t = ideal_triangulation(stokes_graph(PUNC3; theta = 0.4))
        @test tagged_arcs(t) ==
              [TaggedArc(e, (:plain, :plain)) for e in t.edge_endpoints[t.is_diagonal]]
    end

    @testset "the unpunctured layer is untouched" begin
        # no disk or annulus triangulation has a self-folded triangle, and the tagged
        # form there carries nothing the edge endpoints did not already
        for g in (stokes_graph(SchrodingerProblem([-1.0, 0.0, 0.0, 1.0]); theta = 0.3),
                  stokes_graph(mathieu_problem(1.0, 3.0); theta = 0.3))
            t = ideal_triangulation(g)
            @test isempty(selffolded_triangles(t))
            @test all(s -> !is_selffolded(t, s), eachindex(t.triangles))
            @test all(a -> a.tags == (:plain, :plain), tagged_arcs(t))
        end
    end
end
