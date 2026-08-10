# Stokes graph → ideal triangulation → quiver: the combinatorial half of the
# cluster bridge.
#
# A generic Stokes graph with all-simple turning points is dual to an ideal
# triangulation of a marked surface (Iwaki–Nakanishi, Gaiotto–Moore–Neitzke): each
# turning point carries a triangle, each Stokes region an edge of the triangulation -
# a *diagonal* when the region is a strip (two turning points on its boundary), a
# *boundary edge* when it is a half-plane (one turning point). The quiver of the
# triangulation is the exchange matrix the DDP layer's `ddp_seed` expects (`B = −P`,
# ledger item 5 of src/ddp.jl).
#
# DEGENERATE turning points are the same walk, not a second one. An order-`m` point
# emits `m + 2` rays and is dual to an `(m+2)`-gon, so the dual of a degenerate graph
# is a *polygon decomposition* - a positive-dimensional face of the associahedron,
# refined by several triangulations. `_build_decomposition` below is written in terms
# of the valences `m + 2` throughout, never the constant 3, and
# [`ideal_triangulation`](@ref) is its all-triangles case; the general object and its
# refinements live in `src/polygon_decomposition.jl`.
#
# WHICH surface. Every singularity of the quadratic differential `Q dz²` on ℙ¹ - each
# pole of `Q` and the point at infinity - contributes one boundary circle carrying
# `|e + 2|` marked points, where `Q ~ c·ζ^e` locally (`e = −m` at a pole of order `m`,
# `e = degree(prob)` at infinity). A degree-`d` polynomial has one singularity, at
# infinity, with `d + 2` marked points: the disk of the original bridge layer. Two
# irregular poles give the *annulus*, i.e. `Ã(p,q)`; more give more. So poles - not
# higher-order turning points - are the route out of type A (see the M4 scope
# boundary in PLANNING/roadmap.md, corrected by M7a).
#
# PUNCTURES are the case `e = −2`, and they are NOT "a boundary circle with zero
# marked points". A double pole has no marked points and no boundary at all: the
# trajectories spiral into it (`src/stokes_graph.jl`), the rays that reach it all end
# at the same place, and in the planar map that place is a **vertex** of degree = the
# number of rays spiralling in - the puncture's valence in the dual triangulation.
# So a puncture contributes no boundary arcs, no hole face, and exactly one vertex,
# which is where the type-D surfaces come from (a once-punctured `m`-gon is `D_m`).
# Everything else in this file is unchanged, because every count below was already
# written over "singularities" rather than "boundary circles":
#
#   * interior faces = L − n + 2 − nb still, since a puncture removes one hole face
#     and one boundary vertex but adds one interior vertex;
#   * diagonals = n − 2 + nb still, which is FST's `6g + 3b + 3p + m − 6` on ℙ¹;
#   * the divisor identity Σ(zero orders) = M + 2b − 4 still, since a double pole
#     contributes `0 + 2` exactly as an order-m pole contributes `(m−2) + 2`.
#
# What DOES change is the region classification: a Stokes region can now end at a
# puncture instead of on a boundary arc, so a half-plane is 1 turning point and 1
# *end*, a strip 2 turning points and 2 *ends*, where an end is a boundary arc OR a
# visit to a puncture. The once-punctured triangle shows both: at θ = 0.4 its three
# strips each have one arc and one puncture end, while in the chamber θ ≈ 1.5 the
# puncture has valence 2 and one strip reverts to two arcs.
#
# The construction is purely combinatorial - no geometry is reconstructed. The traced
# graph is turned into a planar map (rotation system) whose faces are enumerated by
# the standard dual-walk φ(dart) = cw-next(reverse(dart)); the only numerical input
# beyond the traced polylines is the cyclic order of the ray endpoints on each circle,
# which is exact because distinct Stokes lines of one foliation never cross. Every
# combinatorial invariant is asserted - with `n` turning points, `b` boundary circles
# and `M` marked points, `n = M + 2b − 4`, the interior faces number `2n + 2 − b` and
# the diagonals `(3n − M)/2` - and any failure throws `NonGenericGraph`.
#
# ── Convention ledger ──────────────────────────────────────────────────────────────
#  • BOUNDARY ORIENTATION. A circle is traversed so that the *surface* is on the left.
#    The escape circle at infinity encloses the surface, so that is counterclockwise
#    (increasing `angle(z)`); a pole's circle is enclosed *by* the surface, so that is
#    clockwise (decreasing `angle(z − p)`). This one reversal is the whole difference
#    between the disk case and the general case: the local clockwise dart order at an
#    exit point, (arc-forward, arc-backward, inward ray), is the same either way.
#    Pinned by: the face walk finds exactly `b` hole faces of the right sizes, and the
#    Euler counts above; getting it backwards merges faces and fails immediately.
#  • BOUNDARY COMPONENTS ARE INDEXED `1 = ∞`, then the finite poles in the order of
#    `poles(prob)`. Marked points are numbered by component and, within a component,
#    in that component's own surface-counterclockwise order. A puncture is a component
#    that owns exactly one marked point - itself - and `marked_is_puncture` says which.
#  • A PUNCTURE'S ROTATION IS COUNTERCLOCKWISE, like a turning point's and unlike a
#    pole's boundary circle. That is not an inconsistency: a pole's circle is a *face*
#    boundary traversed with the surface on the left (hence clockwise in the plane),
#    while a puncture is a *vertex* whose darts are ordered by the plane's own
#    counterclockwise sense. The rays are ordered by their first crossing of the pole
#    circle, which is well defined because `r` decreases monotonically along a spiral
#    and Stokes lines never cross. Pinned by: the face walk closes and finds exactly
#    `n − 2 + nb` diagonals on the once-punctured triangle; the opposite sense turns
#    every strip inside out and fails on the first count.

# -- additive StokesGraph helpers ----------------------------------------------------

"""
    infinite_lines(g::StokesGraph) -> Vector{StokesLine}

The Stokes lines of `g` escaping to infinity (endpoint `:infinity`).
"""
infinite_lines(g::StokesGraph) = filter(l -> l.endpoint === :infinity, g.lines)

"""
    boundary_lines(g::StokesGraph) -> Vector{StokesLine}

The Stokes lines of `g` that reach a boundary circle of the surface - escaping to
infinity or running into a pole of `Q`. For a polynomial problem this is
[`infinite_lines`](@ref); every ray of a generic graph is one of these.
"""
boundary_lines(g::StokesGraph) =
    filter(l -> l.endpoint === :infinity || l.endpoint === :pole, g.lines)

"""
    ray_exit_angles(g::StokesGraph) -> Vector

`angle(z)` of the point where each infinite ray of `g` crosses the escape circle, in
the order of [`infinite_lines`](@ref). Distinct Stokes lines never cross, so sorting
these angles gives the exact cyclic order of the rays on the boundary of the disk -
the one numerical datum [`ideal_triangulation`](@ref) reads beyond the topology.
"""
ray_exit_angles(g::StokesGraph) = [angle(l.points[end]) for l in infinite_lines(g)]

# -- the triangulation type ----------------------------------------------------------

"""
    IdealTriangulation

The ideal triangulation of the marked surface dual to a generic Stokes graph
(build with [`ideal_triangulation`](@ref)). The surface is a sphere with
`n_boundaries` boundary circles, carrying `n_marked` marked points in total;
`marked_boundary[i]` is the singularity marked point `i` belongs to,
`marked_is_puncture[i]` says whether it *is* that singularity (a double pole of `Q`,
an interior vertex with no boundary) rather than a point on its boundary circle, and
`triangle_corners[s]`
lists triangle `s`'s corners counterclockwise, aligned with `triangles[s]` so that
edge `triangles[s][i]` runs from corner `i` to corner `i+1`. The corners are not
redundant: once two edges of a triangle share *both* marked points - which is exactly
what happens on the annulus, where both arcs of `Ã(1,1)` join the same two points -
the endpoints alone no longer say which corner two edges meet at. Marked points are
numbered by circle and, within a circle, in that circle's surface-counterclockwise
order (at the asymptotic directions `marked_angles`, display data). One boundary
circle is the disk of a polynomial problem, where the labels `1:n_marked` are the
familiar counterclockwise polygon labels.

Edges are indexed canonically - diagonals first, sorted by `diagonal_tp_pair` (the
pair of turning points carrying the two adjacent triangles) with the turning points'
ray sectors breaking ties, then boundary edges sorted by their endpoints - and
`edge_endpoints[e]` gives edge `e`'s marked-point pair. `triangles[s]` lists the three
edges of the triangle dual to turning point `triangle_tp[s]`, counterclockwise. The
quiver on the diagonals is [`triangulation_quiver`](@ref). `theta` records the phase
of the Stokes graph this triangulation came from - provenance, so the chamber is
recoverable downstream (`Float64` because graph topology is Float64 by design);
[`flip`](@ref) carries it through unchanged.
"""
struct IdealTriangulation
    n_marked::Int
    edge_endpoints::Vector{Tuple{Int,Int}}
    is_diagonal::Vector{Bool}
    triangles::Vector{NTuple{3,Int}}
    triangle_tp::Vector{Int}
    diagonal_tp_pair::Vector{Tuple{Int,Int}}
    marked_angles::Vector{Float64}
    theta::Float64
    marked_boundary::Vector{Int}
    triangle_corners::Vector{NTuple{3,Int}}
    marked_is_puncture::Vector{Bool}
end

# Interior constructor: a surface with no punctures unless said otherwise, so the
# annulus/disk call sites keep their arity.
IdealTriangulation(n_marked, ee, isd, tris, ttp, dtp, ma, θ, mb, tc) =
    IdealTriangulation(n_marked, ee, isd, tris, ttp, dtp, ma, θ, mb, tc,
                       falses(n_marked))

# Disk constructor: one boundary circle, and corners read off the edge endpoints -
# unambiguous on a disk, where no two edges of a triangle share both endpoints.
function IdealTriangulation(n_marked, ee, isd, tris, ttp, dtp, ma, θ)
    corners = map(tris) do tri
        ntuple(3) do i
            a, b = ee[tri[i]]
            c, d = ee[tri[mod1(i - 1, 3)]]
            a == c || a == d ? a : b
        end
    end
    IdealTriangulation(n_marked, ee, isd, tris, ttp, dtp, ma, θ,
                       ones(Int, n_marked), corners)
end

"""
    n_marked_points(t::IdealTriangulation) -> Int

Total number of marked points, i.e. of vertices of the triangulation: `d + 2` for a
degree-`d` polynomial potential, whose surface is a disk. **Punctures are counted** -
they are vertices like any other - so this is (boundary marked points) +
[`n_punctures`](@ref).
"""
n_marked_points(t::IdealTriangulation) = t.n_marked

"""
    n_boundaries(t::IdealTriangulation) -> Int

Number of singularities of `Q dz²`, so `1` for a polynomial problem (the disk) and `2`
for the annulus. Each is a boundary circle unless it is a puncture (see
[`n_punctures`](@ref)), which is why the diagonal count `n − 2 + n_boundaries` is one
formula for both.
"""
n_boundaries(t::IdealTriangulation) = maximum(t.marked_boundary; init = 1)

"""
    n_punctures(t::IdealTriangulation) -> Int

Number of punctures of the surface - marked points that are interior vertices (double
poles of `Q`) rather than points on a boundary circle. `0` for a disk or an annulus.
"""
n_punctures(t::IdealTriangulation) = count(t.marked_is_puncture)

"""
    puncture_valence(t::IdealTriangulation, i::Integer) -> Int

The number of edges of `t` incident to marked point `i`, counting a loop twice. At a
puncture this is its valence in the dual triangulation - the number of Stokes rays
that spiral into it - and `1` is the self-folded case that
[`ideal_triangulation`](@ref) refuses.
"""
puncture_valence(t::IdealTriangulation, i::Integer) =
    sum(count(==(i), e) for e in t.edge_endpoints; init = 0)

"""
    n_diagonals(t::IdealTriangulation) -> Int

Number of diagonals (internal edges) - the rank of [`triangulation_quiver`](@ref).
"""
n_diagonals(t::IdealTriangulation) = count(t.is_diagonal)

"""
    diagonals(t::IdealTriangulation) -> Vector{Int}

The edge indices of the diagonals (canonically `1:n_diagonals(t)`).
"""
diagonals(t::IdealTriangulation) = findall(t.is_diagonal)

# -- small union-find (over ray indices, for gluing triangle corners) ----------------

function _uf_find(parent::Vector{Int}, i::Int)
    while parent[i] != i
        parent[i] = parent[parent[i]]
        i = parent[i]
    end
    i
end

function _uf_union!(parent::Vector{Int}, i::Int, j::Int)
    ri, rj = _uf_find(parent, i), _uf_find(parent, j)
    ri != rj && (parent[ri] = rj)
    nothing
end

# circular distance of two angles
_circ_dist(a, b) = abs(_wrap(a - b))

# A short arc on a boundary circle from `za` to `zb` about `centre`, the short way
# round (the two exits are adjacent positions on the circle, so the short way is the
# strip's own boundary arc). Radius is interpolated too, since the escape circle and a
# pole circle are only nominally circles for the traced endpoints.
function _boundary_arc(za::Complex{F}, zb::Complex{F}, centre::Complex{F};
                       npts::Int = 8) where {F}
    a, b = za - centre, zb - centre
    φ0, dφ = angle(a), _wrap(angle(b) - angle(a))
    r0, r1 = abs(a), abs(b)
    [centre + (r0 + (r1 - r0) * s) * cis(φ0 + dφ * s)
     for s in range(zero(F), one(F), length = npts)[2:(end - 1)]]
end

# -- the construction ----------------------------------------------------------------

"""
    ideal_triangulation(g::StokesGraph) -> IdealTriangulation

The ideal triangulation of the marked surface dual to the generic Stokes graph `g` -
a disk for a polynomial `Q`, a sphere with one boundary circle per pole (plus one at
infinity) for a [`RationalProblem`](@ref).

Requires genericity: every Stokes line must reach a boundary circle (a finite line
means `θ` sits on a wall - perturb it away from the [`saddle_candidates`](@ref)
phases; an `:incomplete` line means the trace gave up). Throws
[`NonGenericGraph`](@ref) otherwise, and whenever a combinatorial invariant of a
generic graph fails.
"""
function ideal_triangulation(g::StokesGraph)
    _require_simple(g)
    _as_triangulation(_build_decomposition(g)[1])
end

# The degeneracy gate of `ideal_triangulation`. The walk itself handles any order
# (see the header) - it is the *triangulation* that does not exist, because an
# order-m point is dual to an (m+2)-gon with several triangulations refining it.
function _require_simple(g::StokesGraph)
    is_degenerate(g) || return nothing
    bad = findfirst(!is_simple, g.turning_points)::Int
    throw(NonGenericGraph(
        "turning point $bad has order $(g.turning_points[bad].order): a degenerate " *
        "Stokes graph is dual to a polygon decomposition, not a triangulation - " *
        "perturb the energy off the critical value, or use polygon_decomposition"))
end

# A decomposition all of whose cells are triangles, as an `IdealTriangulation`.
function _as_triangulation(pd::PolygonDecomposition)
    all(c -> length(c) == 3, pd.cells) || throw(NonGenericGraph(
        "the dual cell decomposition has a $(maximum(length, pd.cells))-gon cell: " *
        "it is not a triangulation"))
    IdealTriangulation(pd.n_marked, pd.edge_endpoints, pd.is_diagonal,
                       [NTuple{3,Int}(c) for c in pd.cells], pd.cell_tp,
                       pd.diagonal_tp_pair, pd.marked_angles, pd.theta,
                       pd.marked_boundary, [NTuple{3,Int}(c) for c in pd.cell_corners],
                       pd.marked_is_puncture)
end

"""
    diagonal_core_paths(g::StokesGraph) -> Vector{Vector{Complex}}

One open polyline per diagonal of `ideal_triangulation(g)`, in the same canonical
order: a path from the first turning point of the diagonal's strip region to the
second, **running inside that strip**. This is the geometric datum a charge cycle
needs and a turning-point pair does not carry - see [`charge_basis`](@ref).

A strip region is a quadrilateral `z_s → M₁ → z_{s'} → M₂` with the two turning points
at *opposite* corners, so its two sides each join `z_s` to `z_{s'}` through one marked
point. The path returned is the cheaper of those two sides (cheaper in accumulated
Stokes [`mass`](@ref), which is what the period integral along it costs): the traced
ray out of `z_s`, the strip's own boundary arc across the circle, and the traced ray
back into `z_{s'}`.

Why not a straight segment: on a punctured surface two diagonals can join the *same*
pair of turning points and differ only by a loop around a hole - on the annulus that
is exactly the difference between the monopole and the dyon - so the pair is not a
name for the cycle and only the strip decides.
"""
diagonal_core_paths(g::StokesGraph) = _build_decomposition(g)[2]

# The shared construction. Returns the dual polygon decomposition and, aligned with its
# diagonals, the strip core paths of `diagonal_core_paths`. Split this way because the
# paths are geometry and must NOT live on the decomposition: `flip` returns a valid
# triangulation and could not maintain them (nor the ray sectors, for the same reason).
function _build_decomposition(g::StokesGraph{F}) where {F}
    n = length(g.turning_points)
    lines = g.lines
    nb = length(g.poles) + 1                    # singularities: ∞, then the poles
    # which singularities are punctures (double poles): vertices, not boundary circles
    is_punc = falses(nb)
    for j in eachindex(g.poles)
        is_punc[1 + j] = g.pole_orders[j] == 2
    end
    n_punc = count(is_punc)

    # 0. valences. An order-m turning point emits m + 2 rays and carries an (m+2)-gon
    # cell; every count below is written with `val`, never the constant 3, which is
    # what makes this walk the general one. `L` is the total number of rays (= darts
    # at turning points = cell sides = boundary exit points).
    val = [order(t) + 2 for t in g.turning_points]
    L = sum(val)

    # 1. genericity gate: all L rays must reach a boundary circle.
    for l in lines
        if l.endpoint === :turning_point
            throw(NonGenericGraph(
                "the graph contains a finite Stokes line between turning points " *
                "$(l.source) and $(l.target): θ = $(Float64(g.theta)) sits on a wall - " *
                "perturb θ away from the saddle_candidates phases"))
        elseif l.endpoint !== :infinity && l.endpoint !== :pole
            throw(NonGenericGraph(
                "Stokes line (direction $(l.direction)) from turning point " *
                "$(l.source) is :$(l.endpoint); retrace without allow_incomplete " *
                "or with a larger max_mass"))
        end
    end
    length(lines) == L || throw(NonGenericGraph(
        "expected $L rays (Σ (order + 2) over the turning points), " *
        "found $(length(lines))"))

    # rays of each turning point, indexed by direction k = 0 … order+1 - this IS the
    # counterclockwise rotation at the turning point (the seed directions
    # (2θ − arg c_m + 2πk)/(m+2) increase counterclockwise; exact, no numerics).
    rays = [zeros(Int, val[s]) for s in 1:n]
    for (ℓ, l) in enumerate(lines)
        rays[l.source][l.direction + 1] = ℓ
    end
    all(r -> all(>(0), r), rays) || throw(NonGenericGraph(
        "some turning point is missing one of its rays"))

    # 2. which circle each ray lands on, and its angle *about that circle's centre*.
    circle = zeros(Int, L)
    ang = zeros(F, L)
    for (ℓ, l) in enumerate(lines)
        if l.endpoint === :infinity
            circle[ℓ] = 1
            ang[ℓ] = angle(l.points[end])
        else
            j = l.target::Int
            circle[ℓ] = 1 + j
            ang[ℓ] = angle(l.points[end] - g.poles[j])
        end
    end

    # Cyclic order of the arrival points at each singularity. On a boundary circle it
    # is surface-counterclockwise: increasing angle at infinity (the surface is
    # inside), decreasing angle at an irregular pole (the surface is outside). At a
    # puncture the rays arrive at a *vertex*, whose rotation is the plane's own
    # counterclockwise sense, i.e. increasing angle - see the ledger at the top of
    # this file. Only boundary circles get exit *positions* and arcs; a puncture's
    # rays are recorded in `pring` and never leave the vertex.
    bpos = Int[]                                # bpos[p] = line at boundary position p
    cranges = Vector{UnitRange{Int}}(undef, nb) # only for boundary circles
    pring = [Int[] for _ in 1:nb]               # only for punctures
    tol = sqrt(eps(F))
    for c in 1:nb
        mem = findall(==(c), circle)
        isempty(mem) && throw(NonGenericGraph(
            "singularity $c (" * (c == 1 ? "infinity" : "pole $(c - 1)") *
            ") receives no Stokes rays - the graph is not generic"))
        ccw = c == 1 || is_punc[c]
        sort!(mem; by = ℓ -> (ccw ? ang[ℓ] : -ang[ℓ]))
        for i in eachindex(mem)
            length(mem) == 1 && break
            a, bb = ang[mem[i]], ang[mem[mod1(i + 1, length(mem))]]
            gap = ccw ? mod(bb - a, 2 * F(π)) : mod(a - bb, 2 * F(π))
            gap ≥ tol || throw(NonGenericGraph(
                "two rays reach singularity $c at nearly the same angle " *
                "(gap ≈ $(Float64(gap))): θ = $(Float64(g.theta)) is numerically on " *
                "a wall - perturb θ"))
        end
        if is_punc[c]
            cranges[c] = 1:0
            pring[c] = mem
        else
            cranges[c] = (length(bpos) + 1):(length(bpos) + length(mem))
            append!(bpos, mem)
        end
    end
    Lb = length(bpos)                           # exit points on boundary circles
    prevpos = zeros(Int, L)
    for c in 1:nb
        r = cranges[c]
        for (i, p) in enumerate(r)
            prevpos[p] = r[mod1(i - 1, length(r))]
        end
    end

    # 3. face walk of the planar map. Darts (half-edges):
    #      ℓ ∈ 1:L         ray outward  (turning point → arrival)
    #      L + ℓ           ray inward   (arrival → turning point)
    #      2L + p          arc forward  (exit position p → next on its circle)
    #      3L + p          arc backward (the reverse)
    # only `p ≤ Lb` exist, since a puncture has no arcs: there the ray-inward darts
    # are based at the puncture VERTEX and rotate among themselves, which is the one
    # structural change punctures make. Faces are the orbits of φ(d) = cwnext[rev(d)]
    # (rotations counterclockwise, faces on the left, interior faces traversed
    # counterclockwise).
    rev(d) = (d ≤ L || (2L < d ≤ 3L)) ? d + L : d - L
    cwnext = zeros(Int, 4L)
    for s in 1:n, k in 1:val[s]
        cwnext[rays[s][k]] = rays[s][mod1(k - 1, val[s])]
    end
    for p in 1:Lb
        pm1 = prevpos[p]
        # clockwise order at exit p: [arc-forward p, arc-backward p−1, ray inward]
        cwnext[2L + p] = 3L + pm1
        cwnext[L + bpos[p]] = 2L + p
        cwnext[3L + pm1] = L + bpos[p]
    end
    for c in 1:nb                               # the puncture vertices
        is_punc[c] || continue
        r = pring[c]
        for i in eachindex(r)
            cwnext[L + r[i]] = L + r[mod1(i - 1, length(r))]
        end
    end

    live = vcat(1:(2L), (2L + 1):(2L + Lb), (3L + 1):(3L + Lb))
    visited = trues(4L)
    visited[live] .= false
    faces = Vector{Vector{Int}}()               # interior faces (lists of darts)
    face_of = zeros(Int, 4L)
    holes = 0
    for d0 in live
        visited[d0] && continue
        orbit = Int[]
        d = d0
        while !visited[d]
            visited[d] = true
            push!(orbit, d)
            d = cwnext[rev(d)]
        end
        d == d0 || throw(NonGenericGraph("face walk failed to close (tracing bug)"))
        if any(x -> x > 3L, orbit)              # arc-backward darts ⇒ a hole face
            (!all(x -> x > 3L, orbit) ||
             !any(c -> !is_punc[c] && length(orbit) == length(cranges[c]), 1:nb)) &&
                throw(NonGenericGraph("a boundary face is malformed - the graph is " *
                                      "not a generic planar Stokes graph"))
            holes += 1
        else
            push!(faces, orbit)
            for d in orbit
                face_of[d] = length(faces)
            end
        end
    end
    holes == nb - n_punc || throw(NonGenericGraph(
        "expected $(nb - n_punc) boundary circles, the face walk found $holes"))
    # Euler on the planar map (V = n + Lb + n_punc vertices, E = L + Lb edges): the
    # interior faces number L − n + 2 − nb, punctures and all - each one trades its
    # hole face and its boundary vertices for a single interior vertex. The all-simple
    # unpunctured case L = 3n is the familiar 2n + 2 − b.
    length(faces) == L - n + 2 - nb || throw(NonGenericGraph(
        "expected $(L - n + 2 - nb) Stokes regions, found $(length(faces))"))

    # 4. classify faces: one turning point on the boundary ⇒ half-plane region ⇒
    # boundary edge; two ⇒ strip region ⇒ diagonal. The other index is the number of
    # ENDS - the places where the region reaches a singularity - and an end is a
    # boundary arc OR a visit to a puncture vertex (an inward ray dart based there).
    # On an unpunctured surface every end is an arc and this is the original rule.
    face_tps = Vector{Vector{Int}}(undef, length(faces))
    for (f, orbit) in enumerate(faces)
        tps = Int[]
        ends = 0
        for d in orbit
            if d ≤ L
                lines[d].source in tps || push!(tps, lines[d].source)
            elseif d ≤ 2L
                ℓ = d - L
                lines[ℓ].source in tps || push!(tps, lines[ℓ].source)
                is_punc[circle[ℓ]] && (ends += 1)
            else
                ends += 1
            end
        end
        sort!(tps)
        face_tps[f] = tps
        ok = (length(tps) == 1 && ends == 1) || (length(tps) == 2 && ends == 2)
        ok || throw(NonGenericGraph(
            "Stokes region with $(length(tps)) turning points and $ends ends " *
            "(boundary arcs and puncture visits) - not a half-plane or strip region " *
            "of a generic graph"))
    end
    isdiag = [length(t) == 2 for t in face_tps]
    # boundary marked points, read off the problem: |e + 2| per singularity, which is
    # 0 at a puncture (see the header). The vertices of the triangulation are these
    # plus one per puncture.
    nbmarked = sum(length(asymptotic_directions(g.problem, g.theta; pole = c - 1))
                   for c in 1:nb)
    nmarked = nbmarked + n_punc
    # The divisor identity for `Q dz²` on ℙ¹: Σ (zero orders) = M + 2b − 4. The
    # all-simple case is the familiar n = M + 2b − 4.
    total_order = L - 2n
    total_order == nbmarked + 2nb - 4 || throw(NonGenericGraph(
        "turning-point orders summing to $total_order but $nbmarked boundary marked " *
        "points on $nb singularities: the divisor identity Σ orders = M + 2b − 4 fails"))
    # h + 2d = L (each interior face contributes 1 or 2 ends, and the ends are the L
    # rays: every ray reaches a singularity exactly once) together with
    # h + d = L − n + 2 − b gives d = n − 2 + b, independent of the valences and of
    # how many singularities are punctures: a degeneration merges cells, it does not
    # remove diagonals from the count. This is FST's 6g + 3b + 3p + m − 6 on ℙ¹.
    ndiag = n - 2 + nb
    count(isdiag) == ndiag || throw(NonGenericGraph(
        "expected $ndiag strip regions (diagonals), found $(count(isdiag))"))

    # 5. cells. The corner of turning point s between its counterclockwise-consecutive
    # rays (r_k, r_{k+1}) lies in the face of the outward dart r_k, so cell s has
    # edge-faces (F(r₁), …, F(r_{m+2})) counterclockwise, and the corner between
    # consecutive edge-faces (F(r_k), F(r_{k+1})) sits at the shared ray r_{k+1} -
    # corners are labelled by rays (L corners, one per ray).
    cell_faces = Vector{Vector{Int}}(undef, n)
    for s in 1:n
        fs = [face_of[r] for r in rays[s]]
        allunique(fs) || throw(NonGenericGraph(
            "two sectors of turning point $s lie in the same Stokes region"))
        cell_faces[s] = fs
    end

    # 6. glue cells along diagonals: identify corners with reversed orientation.
    # If diagonal f occupies sector (r_a, r_{a+1}) of s and (r_b, r_{b+1}) of s′, the
    # gluing is corner(r_a) ~ corner(r′_{b+1}) and corner(r_{a+1}) ~ corner(r′_b).
    sector_index(s, f) = begin
        idx = findall(==(f), cell_faces[s])
        length(idx) == 1 || throw(NonGenericGraph(
            "the Stokes region dual to an edge is adjacent to turning point $s more " *
            "than once: its dual triangle has the same edge twice, i.e. a " *
            "SELF-FOLDED triangle at a puncture of valence 1. That needs a tagged " *
            "triangulation ([FST08]), which this layer does not have - perturb θ " *
            "into a neighbouring chamber, where the puncture has valence ≥ 2"))
        idx[1]
    end
    parent = collect(1:L)                       # union-find over rays (= corners)
    # Every ray spiralling into a puncture has its corner AT the puncture, so those
    # corners are one class before any gluing. The gluing along diagonals would find
    # most of this on its own - around a puncture consecutive rays are separated by
    # strips - but not all of it, and the puncture is a fact about the problem rather
    # than an output of the walk.
    for c in 1:nb
        is_punc[c] || continue
        for i in 2:length(pring[c])
            _uf_union!(parent, pring[c][1], pring[c][i])
        end
    end
    for f in eachindex(faces)
        isdiag[f] || continue
        s, s′ = face_tps[f]
        a, b = sector_index(s, f), sector_index(s′, f)
        _uf_union!(parent, rays[s][a], rays[s′][mod1(b + 1, val[s′])])
        _uf_union!(parent, rays[s][mod1(a + 1, val[s])], rays[s′][b])
    end
    roots = unique(_uf_find(parent, ℓ) for ℓ in 1:L)
    length(roots) == nmarked || throw(NonGenericGraph(
        "expected $nmarked marked points, found $(length(roots))"))

    # 7. marked points ↔ the asymptotic directions of their boundary circle. Each
    # class of rays must sit on one singularity and match one of its directions. A
    # puncture has no directions to match: its single class IS the marked point, and
    # the angle recorded for it is the mean arrival angle - display data only, since
    # a spiral has no asymptotic direction (that is the whole point).
    θ = g.theta
    phis = [asymptotic_directions(g.problem, θ; pole = c - 1) for c in 1:nb]
    class_circle = Dict{Int,Int}()
    class_phi = Dict{Int,F}()
    for r in roots
        members = [ℓ for ℓ in 1:L if _uf_find(parent, ℓ) == r]
        c = circle[members[1]]
        all(ℓ -> circle[ℓ] == c, members) || throw(NonGenericGraph(
            "a marked point's rays reach two different singularities - retrace " *
            "with a larger escape_radius or a smaller pole_radius"))
        mean = angle(sum(cis(ang[ℓ]) for ℓ in members))
        class_circle[r] = c
        if is_punc[c]
            class_phi[r] = mean
            continue
        end
        k = argmin([_circ_dist(mean, φ) for φ in phis[c]])
        for ℓ in members
            _circ_dist(ang[ℓ], phis[c][k]) < F(π) / length(phis[c]) ||
                throw(NonGenericGraph(
                    "a ray's exit angle is not close to its marked point's " *
                    "asymptotic direction - retrace with a larger escape_radius " *
                    "or a smaller pole_radius"))
        end
        class_phi[r] = phis[c][k]
    end
    for c in 1:nb
        is_punc[c] || continue
        count(r -> class_circle[r] == c, roots) == 1 || throw(NonGenericGraph(
            "the rays reaching the puncture at $(g.poles[c - 1]) fall into more than " *
            "one marked point - the gluing and the trace disagree"))
    end
    length(unique((class_circle[r], class_phi[r]) for r in roots)) == nmarked ||
        throw(NonGenericGraph(
            "two marked points matched the same asymptotic direction - retrace " *
            "with a larger escape_radius"))
    # number by singularity, then in its own surface-counterclockwise order
    roots_sorted = sort(roots;
        by = r -> (class_circle[r],
                   class_circle[r] == 1 ? class_phi[r] : -class_phi[r]))
    marked_label = Dict(r => i for (i, r) in enumerate(roots_sorted))
    marked_angles = Float64[Float64(class_phi[r]) for r in roots_sorted]
    marked_boundary = Int[class_circle[r] for r in roots_sorted]
    marked_is_puncture = Bool[is_punc[class_circle[r]] for r in roots_sorted]

    # 8. edge endpoints: the two corners flanking the sector a face occupies.
    endpoint_pair = Vector{Tuple{Int,Int}}(undef, length(faces))
    for f in eachindex(faces)
        s = face_tps[f][1]
        a = sector_index(s, f)
        u = marked_label[_uf_find(parent, rays[s][a])]
        v = marked_label[_uf_find(parent, rays[s][mod1(a + 1, val[s])])]
        endpoint_pair[f] = (min(u, v), max(u, v))
    end

    # 9. canonical edge order: diagonals sorted by turning-point pair, then boundary
    # edges sorted by endpoints. The pair alone is NOT a key once the surface has a
    # hole - on the annulus both arcs of Ã(1,1) join the same two turning points -
    # so the two ray sectors the strip occupies break the tie. Sectors are indexed by
    # the exact ray directions at a turning point, so this stays canonical.
    diag_key(f) = (face_tps[f][1], face_tps[f][2],
                   sector_index(face_tps[f][1], f), sector_index(face_tps[f][2], f))
    diag_faces = findall(isdiag)
    allunique(diag_key(f) for f in diag_faces) || throw(NonGenericGraph(
        "two strip regions share the same turning-point pair and sectors"))
    sort!(diag_faces; by = diag_key)
    bdry_faces = sort!(findall(!, isdiag); by = f -> endpoint_pair[f])
    new_id = zeros(Int, length(faces))
    for (i, f) in enumerate(diag_faces)
        new_id[f] = i
    end
    for (i, f) in enumerate(bdry_faces)
        new_id[f] = ndiag + i
    end
    edge_endpoints = [endpoint_pair[f] for f in vcat(diag_faces, bdry_faces)]
    is_diagonal = vcat(trues(ndiag), falses(length(bdry_faces)))
    diagonal_tp_pair = [(face_tps[f][1], face_tps[f][2]) for f in diag_faces]
    # Edges and corners of a cell are aligned: edge slot i is the face of ray i, whose
    # two flanking corners are the classes of rays i and i+1 (step 8), so corner slot
    # i is the class of ray i and edge i runs from corner i to corner i+1.
    cells = Vector{Vector{Int}}(undef, n)
    corners = Vector{Vector{Int}}(undef, n)
    for s in 1:n
        v = val[s]
        e = [new_id[f] for f in cell_faces[s]]
        c = [marked_label[_uf_find(parent, r)] for r in rays[s]]
        shift = argmin(e) - 1                   # rotate to start at the smallest edge
        cells[s] = [e[mod1(i + shift, v)] for i in 1:v]
        corners[s] = [c[mod1(i + shift, v)] for i in 1:v]
    end

    tri = PolygonDecomposition(nmarked, edge_endpoints, is_diagonal, cells,
                               collect(1:n), diagonal_tp_pair, marked_angles,
                               Float64(g.theta), marked_boundary, corners,
                               marked_is_puncture)

    # 10. strip core paths (see `diagonal_core_paths`). A strip is a quadrilateral
    # z_s → M₁ → z_{s'} → M₂ with the turning points at OPPOSITE corners, so each of
    # its two sides already joins them through one marked point. Step 6's gluing names
    # the sides: corner(r_a of s) ~ corner(r_{b+1} of s′) is one, corner(r_{a+1} of s)
    # ~ corner(r_b of s′) the other. Take the cheaper in accumulated Stokes mass - the
    # two legs contribute ±e^{iθ}t to the period and partly cancel, so the excursion is
    # what the integral costs.
    tpz = [tp.z for tp in g.turning_points]
    core_paths = Vector{Vector{Complex{F}}}(undef, length(diag_faces))
    for (i, f) in enumerate(diag_faces)
        s, s′ = face_tps[f]
        a, b = sector_index(s, f), sector_index(s′, f)
        sides = ((rays[s][a], rays[s′][mod1(b + 1, val[s′])]),
                 (rays[s][mod1(a + 1, val[s])], rays[s′][b]))
        cost(sd) = lines[sd[1]].mass + lines[sd[2]].mass
        p, q = cost(sides[1]) ≤ cost(sides[2]) ? sides[1] : sides[2]
        circle[p] == circle[q] || throw(NonGenericGraph(
            "the two rays bounding a side of the strip region of diagonal $i reach " *
            "different singularities - the gluing and the trace disagree"))
        centre = circle[p] == 1 ? zero(Complex{F}) : g.poles[circle[p] - 1]
        path = Complex{F}[tpz[s]]
        append!(path, lines[p].points)
        append!(path, _boundary_arc(lines[p].points[end], lines[q].points[end], centre))
        append!(path, reverse(lines[q].points))
        push!(path, tpz[s′])
        core_paths[i] = path
    end

    tri, core_paths
end

# -- the quiver ----------------------------------------------------------------------

"""
    triangulation_quiver(t::IdealTriangulation) -> ClusterAlgebras.Quiver

The adjacency quiver of the triangulation: one vertex per diagonal (in the canonical
order of `t`, labelled `"γ(i,j)"` by the turning-point pair), and for each triangle
one arrow between each counterclockwise-consecutive pair of its diagonals. Internal
triangles always come out as oriented 3-cycles, so this quiver is finite type
`A_{d−1}` for every chamber.

It equals the physical seed quiver of [`charge_basis`](@ref) **exactly**:
`B == −signed_pairing(cb)`, the tightened keystone. (The original bridge layer could only assert
`abs.(B) == abs.(P)` against the raw decay pairing; the signed frame of
`src/signed_frame.jl` supplies the missing orientation.)
"""
function triangulation_quiver(t::IdealTriangulation)
    m = n_diagonals(t)
    B = zeros(Int, m, m)
    for tri in t.triangles
        for c in 1:3
            i, j = tri[c], tri[mod1(c + 1, 3)]
            (t.is_diagonal[i] && t.is_diagonal[j]) || continue
            B[i, j] += 1
            B[j, i] -= 1
        end
    end
    labels = ["γ($(p[1]),$(p[2]))" for p in t.diagonal_tp_pair]
    ClusterAlgebras.Quiver(B, m, ones(Int, m), labels)
end

# -- flips ---------------------------------------------------------------------------
#
# Crossing a wall in θ flips the diagonals whose `diagonal_tp_pair` is a saddle pair of
# that wall (one diagonal generically; several commuting ones at a degenerate wall) and
# mutates the quiver at those vertices. `flip` performs the move combinatorially, so a
# chamber walk costs no Stokes traces.
#
# Edge indices are deliberately PRESERVED - the flipped diagonal keeps index `k`, so
# `triangulation_quiver(flip(t, k)).B == mutate(triangulation_quiver(t), k).B` is an
# index-aligned identity. The canonical order of `ideal_triangulation` (diagonals
# sorted by `diagonal_tp_pair`) is generally a different indexing, because the flip
# permutes the labels of two of the four quad sides; `canonical_reorder` applies that
# permutation when a flipped triangulation must be compared with a freshly traced one.

"""
    flip(t::IdealTriangulation, k::Integer; direction = 1) -> IdealTriangulation

The ideal triangulation obtained by flipping diagonal `k`: the two triangles adjacent
to `k` form a quad, and `k` is replaced by the quad's other diagonal. Edge indices are
preserved (the new diagonal keeps index `k`), so the quiver mutates at vertex `k` -
`triangulation_quiver(flip(t, k)).B == ClusterAlgebras.mutate(triangulation_quiver(t), k).B`.

`direction` is the sense in which the wall is crossed: `+1` for increasing `θ`, `-1`
for decreasing. It is *not* redundant - the Stokes graph rotates with `θ`, so the two
senses reconnect the collapsing strip oppositely and exchange which of the two new
triangles is dual to which turning point. `flip(flip(t, k), k; direction = -1) == t`.

`diagonal_tp_pair` is recomputed: the flipped diagonal keeps its own turning-point
pair (the strip re-forms between the same two turning points), while two of the quad's
four sides swap which turning point borders them. `theta` is carried through unchanged -
it records where the triangulation was traced, not where the flip put it. Use
[`canonical_reorder`](@ref) to compare the result with a freshly traced triangulation.
"""
function flip(t::IdealTriangulation, k::Integer; direction::Integer = 1)
    (1 ≤ k ≤ length(t.edge_endpoints) && t.is_diagonal[k]) || throw(
        Resurgence.InvalidArgument(
            "flip: edge $k is not a diagonal of the triangulation (diagonals are " *
            "$(diagonals(t)))"))
    direction in (-1, 1) || throw(
        Resurgence.InvalidArgument("flip: direction must be ±1, got $direction"))

    adj = findall(tri -> k in tri, t.triangles)
    length(adj) == 2 || throw(NonGenericGraph(
        "diagonal $k borders $(length(adj)) triangles, expected 2"))
    s1, s2 = adj

    # Orient the quad (u, p, v, q) counterclockwise from the triangles' OWN
    # counterclockwise edge/corner lists, not from the marked-point labels: triangle
    # (u, p, v) lists edges (u–p, p–v, v–u) against corners (u, p, v), so with `k` in
    # edge slot `i` the corners fall out as v = c[i], u = c[i+1], p = c[i+2]. This
    # works on any surface; the original rule `u < p < v` read the labels as a convex
    # polygon's counterclockwise order and is meaningless once the surface has a hole,
    # and reading corners off `edge_endpoints` is ambiguous the moment two edges of a
    # triangle share both endpoints (the annulus).
    i1 = findfirst(==(k), collect(t.triangles[s1]))::Int
    i2 = findfirst(==(k), collect(t.triangles[s2]))::Int
    e_up = t.triangles[s1][mod1(i1 + 1, 3)]     # u - p
    e_pv = t.triangles[s1][mod1(i1 + 2, 3)]     # p - v
    e_vq = t.triangles[s2][mod1(i2 + 1, 3)]     # v - q
    e_qu = t.triangles[s2][mod1(i2 + 2, 3)]     # q - u
    v, u, p = t.triangle_corners[s1][i1], t.triangle_corners[s1][mod1(i1 + 1, 3)],
              t.triangle_corners[s1][mod1(i1 + 2, 3)]
    q = t.triangle_corners[s2][mod1(i2 + 2, 3)]
    other = (t.triangle_corners[s2][i2], t.triangle_corners[s2][mod1(i2 + 1, 3)])
    other == (u, v) || throw(NonGenericGraph(
        "the triangles adjacent to diagonal $k are not consistently oriented: it " *
        "runs $u → $v in triangle $s1 but $(other[1]) → $(other[2]) in triangle $s2"))

    endpoints = copy(t.edge_endpoints)
    endpoints[k] = (min(p, q), max(p, q))

    # the two new triangles, counterclockwise: (u, p, q) and (p, v, q)
    rotate(e) = argmin(collect(e)) - 1
    rot(e, shift) = ntuple(i -> e[mod1(i + shift, 3)], 3)
    triangles = copy(t.triangles)
    corners = copy(t.triangle_corners)
    # `triangle_tp` is unchanged - each slot keeps its turning point - so the handedness
    # of the reconnection is expressed by WHICH slot each new triangle lands in. Crossing
    # the wall towards increasing θ sends the u-side triangle to the slot whose apex was
    # p; crossing the other way sends it to the slot whose apex was q. Pinned in both
    # senses by `canonical_reorder(flip(t, k; direction)) == t′` against a re-traced
    # neighbouring chamber, at every wall of the cubic/quartic/quintic fixtures.
    a_slot, b_slot = direction == 1 ? (s1, s2) : (s2, s1)
    sa = rotate((e_up, k, e_qu))
    triangles[a_slot] = rot((e_up, k, e_qu), sa)
    corners[a_slot] = rot((u, p, q), sa)
    sb = rotate((e_pv, e_vq, k))
    triangles[b_slot] = rot((e_pv, e_vq, k), sb)
    corners[b_slot] = rot((p, v, q), sb)
    # A flip that lowers a puncture's valence to 1 produces a self-folded triangle,
    # which this layer cannot represent (see `sector_index`), so refuse it here rather
    # than hand back a triangle with a repeated edge.
    for (s, tri) in enumerate(triangles)
        allunique(tri) || throw(NonGenericGraph(
            "flipping diagonal $k makes triangle $s self-folded (its edges are " *
            "$(tri)): the puncture it wraps would have valence 1, which needs a " *
            "tagged triangulation ([FST08])"))
    end
    tp_pair = _diagonal_tp_pairs(endpoints, t.is_diagonal, triangles, t.triangle_tp)

    IdealTriangulation(t.n_marked, endpoints, t.is_diagonal, triangles, t.triangle_tp,
                       tp_pair, t.marked_angles, t.theta, t.marked_boundary, corners,
                       t.marked_is_puncture)
end

# a diagonal's turning-point pair = the turning points of the two triangles it borders
function _diagonal_tp_pairs(endpoints, is_diagonal, triangles, triangle_tp)
    pairs = Tuple{Int,Int}[]
    for e in eachindex(endpoints)
        is_diagonal[e] || continue
        tps = sort!([triangle_tp[s] for s in eachindex(triangles) if e in triangles[s]])
        length(tps) == 2 || throw(NonGenericGraph(
            "diagonal $e borders $(length(tps)) triangles, expected 2"))
        push!(pairs, (tps[1], tps[2]))
    end
    pairs
end

"""
    canonical_reorder(t::IdealTriangulation) -> (IdealTriangulation, Vector{Int})

Re-index the edges of `t` into the canonical order of [`ideal_triangulation`](@ref) -
diagonals first sorted by `diagonal_tp_pair`, then boundary edges sorted by endpoints -
and return the reordered triangulation together with the permutation `perm`, where
`perm[new] = old`. A freshly traced triangulation is already canonical; the output of
[`flip`](@ref) generally is not.

The traced order breaks a tie between two diagonals with the same turning-point pair
by their ray *sectors*, which are not recoverable from an `IdealTriangulation`; here
the tie is broken by the edge endpoints, and beyond that the (stable) sort keeps the
input order. On a disk, where turning-point pairs are already distinct, the two orders
always agree.
"""
function canonical_reorder(t::IdealTriangulation)
    diag = findall(t.is_diagonal)
    pair_of = Dict(e => p for (e, p) in zip(diag, t.diagonal_tp_pair))
    diag = sort(diag; by = e -> (pair_of[e], t.edge_endpoints[e]))
    bdry = sort(findall(!, t.is_diagonal); by = e -> t.edge_endpoints[e])
    perm = vcat(diag, bdry)
    new_id = zeros(Int, length(perm))
    for (i, e) in enumerate(perm)
        new_id[e] = i
    end
    triangles = similar(t.triangles)
    corners = similar(t.triangle_corners)
    for s in eachindex(t.triangles)
        e = map(x -> new_id[x], collect(t.triangles[s]))
        shift = argmin(e) - 1
        triangles[s] = ntuple(i -> e[mod1(i + shift, 3)], 3)
        corners[s] = ntuple(i -> t.triangle_corners[s][mod1(i + shift, 3)], 3)
    end
    # `is_diagonal` is permuted like every other edge-indexed field. That is a no-op
    # when the input already lists diagonals first (a traced triangulation, or the
    # index-preserving output of `flip`) and load-bearing when it does not - the
    # refinements of a polygon decomposition append their new diagonals at the end.
    tri_out = IdealTriangulation(t.n_marked, t.edge_endpoints[perm], t.is_diagonal[perm],
                                 triangles, t.triangle_tp,
                                 [pair_of[e] for e in diag], t.marked_angles, t.theta,
                                 t.marked_boundary, corners, t.marked_is_puncture)
    tri_out, perm
end
