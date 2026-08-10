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

        # A self-folded triangle (puncture valence 1) needs a tagged triangulation.
        # `flip` will not manufacture one silently: flipping a diagonal of the
        # valence-2 chamber that would drop the valence to 1 is refused.
        t2 = ideal_triangulation(stokes_graph(PUNC3; theta = 1.5))
        p2 = findfirst(t2.marked_is_puncture)
        radii = [k for k in 1:n_diagonals(t2) if p2 in t2.edge_endpoints[k]]
        @test length(radii) == 2
        @test any(radii) do k
            try
                flip(t2, k)
                false
            catch err
                err isa ExactWKB.NonGenericGraph &&
                    occursin("self-folded", sprint(showerror, err))
            end
        end
    end
end
