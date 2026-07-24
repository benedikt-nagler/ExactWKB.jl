# Stokes graph → ideal triangulation → quiver: the combinatorial half of the
# cluster bridge.
#
# A generic Stokes graph of a degree-d polynomial Q with all-simple turning points,
# drawn on the disk (the z-plane compactified at the escape circle), is dual to an
# ideal triangulation of a (d+2)-gon (Iwaki–Nakanishi): each turning point carries a
# triangle, each Stokes region an edge of the triangulation - a *diagonal* when the
# region is a strip (two turning points on its boundary), a *boundary edge* when it is
# a half-plane (one turning point). The quiver of the triangulation is the exchange
# matrix the DDP layer's `ddp_seed` expects (`B = −P`, ledger item 5 of src/ddp.jl).
#
# The construction is purely combinatorial - no geometry is reconstructed. The traced
# graph is turned into a planar map (rotation system) whose faces are enumerated by
# the standard dual-walk φ(dart) = cw-next(reverse(dart)); the only numerical input
# beyond the traced polylines is the cyclic order of the ray endpoints on the escape
# circle, which is exact because distinct Stokes lines of one foliation never cross.
# Every combinatorial invariant of a generic graph (region count 2d+1, marked-point
# count d+2, one-or-two turning points per region, …) is asserted, and any failure
# throws `NonGenericGraph`.

# -- additive StokesGraph helpers ----------------------------------------------------

"""
    infinite_lines(g::StokesGraph) -> Vector{StokesLine}

The Stokes lines of `g` escaping to infinity (endpoint `:infinity`).
"""
infinite_lines(g::StokesGraph) = filter(l -> l.endpoint === :infinity, g.lines)

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

The ideal triangulation of the `n_marked`-gon dual to a generic Stokes graph
(build with [`ideal_triangulation`](@ref)). Marked points are labelled `1:n_marked`
counterclockwise (at the asymptotic directions `marked_angles`, display data). Edges
are indexed canonically - diagonals first, sorted by `diagonal_tp_pair` (the pair of
turning points carrying the two adjacent triangles), then boundary edges sorted by
their endpoints - and `edge_endpoints[e]` gives edge `e`'s marked-point pair.
`triangles[s]` lists the three edges of the triangle dual to turning point
`triangle_tp[s]`, counterclockwise. The quiver on the diagonals is
[`triangulation_quiver`](@ref). `theta` records the phase of the Stokes graph this
triangulation came from - provenance, so the chamber is recoverable downstream
(`Float64` because graph topology is Float64 by design); [`flip`](@ref) carries it
through unchanged.
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
end

"""
    n_marked_points(t::IdealTriangulation) -> Int

Number of marked points of the polygon (`= d + 2` for a degree-`d` potential).
"""
n_marked_points(t::IdealTriangulation) = t.n_marked

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

# -- the construction ----------------------------------------------------------------

"""
    ideal_triangulation(g::StokesGraph) -> IdealTriangulation

The ideal triangulation of the marked disk dual to the generic Stokes graph `g`.
Requires genericity: every Stokes line escapes to infinity (a finite line means `θ`
sits on a wall - perturb it away from the [`saddle_candidates`](@ref) phases; an
`:incomplete` line means the trace gave up). Throws [`NonGenericGraph`](@ref)
otherwise, and whenever a combinatorial invariant of a generic graph fails.
"""
function ideal_triangulation(g::StokesGraph{F}) where {F}
    n = length(g.turning_points)
    lines = g.lines

    # 1. genericity gate: all 3n rays must escape.
    for l in lines
        if l.endpoint === :turning_point
            throw(NonGenericGraph(
                "the graph contains a finite Stokes line between turning points " *
                "$(l.source) and $(l.target): θ = $(Float64(g.theta)) sits on a wall - " *
                "perturb θ away from the saddle_candidates phases"))
        elseif l.endpoint !== :infinity
            throw(NonGenericGraph(
                "Stokes line (direction $(l.direction)) from turning point " *
                "$(l.source) is :$(l.endpoint); retrace without allow_incomplete " *
                "or with a larger max_mass"))
        end
    end
    length(lines) == 3n || throw(NonGenericGraph(
        "expected 3·$n rays, found $(length(lines))"))

    # rays of each turning point, indexed by direction k = 0, 1, 2 - this IS the
    # counterclockwise rotation at the turning point (the seed directions
    # (2θ − arg Q′ + 2πk)/3 increase counterclockwise; exact, no numerics).
    rays = [zeros(Int, 3) for _ in 1:n]
    for (ℓ, l) in enumerate(lines)
        rays[l.source][l.direction + 1] = ℓ
    end
    all(r -> all(>(0), r), rays) || throw(NonGenericGraph(
        "some turning point is missing one of its three rays"))

    # 2. cyclic order of the ray endpoints on the escape circle.
    ang = [angle(l.points[end]) for l in lines]
    order = sortperm(ang)                       # order[p] = line at boundary position p
    posof = invperm(order)
    tol = sqrt(eps(F))
    for p in 1:3n
        gap = p == 3n ? (ang[order[1]] + 2 * F(π) - ang[order[3n]]) :
                        (ang[order[p + 1]] - ang[order[p]])
        gap ≥ tol || throw(NonGenericGraph(
            "two rays exit the escape circle at nearly the same angle " *
            "(gap ≈ $(Float64(gap))): θ = $(Float64(g.theta)) is numerically on a " *
            "wall - perturb θ"))
    end

    # 3. face walk of the planar map. Darts (half-edges), 12n total:
    #      ℓ ∈ 1:3n        ray outward  (turning point → exit)
    #      3n + ℓ          ray inward   (exit → turning point)
    #      6n + p          arc forward  (exit position p → p+1, counterclockwise)
    #      9n + p          arc backward (exit position p+1 → p)
    # Faces are the orbits of φ(d) = cwnext[rev(d)] (rotations counterclockwise,
    # faces on the left, interior faces traversed counterclockwise).
    rev(d) = (d ≤ 3n || (6n < d ≤ 9n)) ? d + 3n : d - 3n
    cwnext = zeros(Int, 12n)
    for s in 1:n, k in 1:3
        cwnext[rays[s][k]] = rays[s][mod1(k - 1, 3)]
    end
    for p in 1:3n
        pm1 = mod1(p - 1, 3n)
        # counterclockwise order at exit p: [arc-forward p, ray inward, arc-backward p−1]
        cwnext[6n + p] = 9n + pm1
        cwnext[3n + order[p]] = 6n + p
        cwnext[9n + pm1] = 3n + order[p]
    end

    visited = falses(12n)
    faces = Vector{Vector{Int}}()               # interior faces (lists of darts)
    face_of = zeros(Int, 12n)
    outer_seen = false
    for d0 in 1:12n
        visited[d0] && continue
        orbit = Int[]
        d = d0
        while !visited[d]
            visited[d] = true
            push!(orbit, d)
            d = cwnext[rev(d)]
        end
        d == d0 || throw(NonGenericGraph("face walk failed to close (tracing bug)"))
        if any(x -> x > 9n, orbit)              # arc-backward darts ⇒ the outer face
            (outer_seen || length(orbit) != 3n || !all(x -> x > 9n, orbit)) &&
                throw(NonGenericGraph("outer face is malformed - the graph is not a " *
                                      "generic planar Stokes graph"))
            outer_seen = true
        else
            push!(faces, orbit)
            for d in orbit
                face_of[d] = length(faces)
            end
        end
    end
    length(faces) == 2n + 1 || throw(NonGenericGraph(
        "expected $(2n + 1) Stokes regions, found $(length(faces))"))

    # 4. classify faces: one turning point on the boundary ⇒ half-plane region ⇒
    # boundary edge; two ⇒ strip region ⇒ diagonal.
    face_tps = Vector{Vector{Int}}(undef, length(faces))
    for (f, orbit) in enumerate(faces)
        tps = Int[]
        arcs = 0
        for d in orbit
            if d ≤ 3n
                lines[d].source in tps || push!(tps, lines[d].source)
            elseif d ≤ 6n
                lines[d - 3n].source in tps || push!(tps, lines[d - 3n].source)
            else
                arcs += 1
            end
        end
        sort!(tps)
        face_tps[f] = tps
        ok = (length(tps) == 1 && arcs == 1) || (length(tps) == 2 && arcs == 2)
        ok || throw(NonGenericGraph(
            "Stokes region with $(length(tps)) turning points and $arcs boundary " *
            "arcs - not a half-plane or strip region of a generic graph"))
    end
    isdiag = [length(t) == 2 for t in face_tps]
    count(isdiag) == n - 1 || throw(NonGenericGraph(
        "expected $(n - 1) strip regions (diagonals), found $(count(isdiag))"))

    # 5. triangles. The corner of turning point s between its counterclockwise-
    # consecutive rays (r_k, r_{k+1}) lies in the face of the outward dart r_k, so
    # triangle s has edge-faces (F(r₁), F(r₂), F(r₃)) counterclockwise, and the
    # corner between consecutive edge-faces (F(r_k), F(r_{k+1})) sits at the shared
    # ray r_{k+1} - corners are labelled by rays (3n corners, one per ray).
    tri_faces = Vector{NTuple{3,Int}}(undef, n)
    for s in 1:n
        fs = (face_of[rays[s][1]], face_of[rays[s][2]], face_of[rays[s][3]])
        allunique(fs) || throw(NonGenericGraph(
            "two sectors of turning point $s lie in the same Stokes region"))
        tri_faces[s] = fs
    end

    # 6. glue triangles along diagonals: identify corners with reversed orientation.
    # If diagonal f occupies sector (r_a, r_{a+1}) of s and (r_b, r_{b+1}) of s′, the
    # gluing is corner(r_a) ~ corner(r′_{b+1}) and corner(r_{a+1}) ~ corner(r′_b).
    sector_index(s, f) = begin
        idx = findall(==(f), collect(tri_faces[s]))
        length(idx) == 1 || throw(NonGenericGraph(
            "diagonal region is adjacent to turning point $s more than once " *
            "(self-folded triangle - impossible on the disk)"))
        idx[1]
    end
    parent = collect(1:3n)                      # union-find over rays (= corners)
    for f in eachindex(faces)
        isdiag[f] || continue
        s, s′ = face_tps[f]
        a, b = sector_index(s, f), sector_index(s′, f)
        _uf_union!(parent, rays[s][a], rays[s′][mod1(b + 1, 3)])
        _uf_union!(parent, rays[s][mod1(a + 1, 3)], rays[s′][b])
    end
    roots = unique(_uf_find(parent, ℓ) for ℓ in 1:3n)
    length(roots) == n + 2 || throw(NonGenericGraph(
        "expected $(n + 2) marked points, found $(length(roots))"))

    # 7. marked points ↔ the d+2 asymptotic directions φ_k = (2θ − arg a_d + 2πk)/(d+2).
    qc = q_coefficients(g.problem)
    d = length(qc) - 1
    d == n || throw(NonGenericGraph(
        "degree $d potential with $n simple turning points - the graph is degenerate"))
    θ = g.theta
    phis = [mod((2θ - angle(Complex{F}(last(qc))) + 2 * F(π) * k) / (d + 2), 2 * F(π))
            for k in 0:(d + 1)]
    class_phi = Dict{Int,F}()
    for r in roots
        members = [ℓ for ℓ in 1:3n if _uf_find(parent, ℓ) == r]
        mean = angle(sum(cis(ang[ℓ]) for ℓ in members))
        k = argmin([_circ_dist(mean, φ) for φ in phis])
        for ℓ in members
            _circ_dist(ang[ℓ], phis[k]) < F(π) / (d + 2) || throw(NonGenericGraph(
                "a ray's exit angle is not close to its marked point's asymptotic " *
                "direction - retrace with a larger escape_radius"))
        end
        class_phi[r] = phis[k]
    end
    length(unique(values(class_phi))) == n + 2 || throw(NonGenericGraph(
        "two marked points matched the same asymptotic direction - retrace with a " *
        "larger escape_radius"))
    roots_sorted = sort(roots; by = r -> class_phi[r])
    marked_label = Dict(r => i for (i, r) in enumerate(roots_sorted))
    marked_angles = Float64[Float64(class_phi[r]) for r in roots_sorted]

    # 8. edge endpoints: the two corners flanking the sector a face occupies.
    endpoint_pair = Vector{Tuple{Int,Int}}(undef, length(faces))
    for f in eachindex(faces)
        s = face_tps[f][1]
        a = sector_index(s, f)
        u = marked_label[_uf_find(parent, rays[s][a])]
        v = marked_label[_uf_find(parent, rays[s][mod1(a + 1, 3)])]
        endpoint_pair[f] = (min(u, v), max(u, v))
    end

    # 9. canonical edge order: diagonals sorted by turning-point pair, then boundary
    # edges sorted by endpoints.
    diag_faces = findall(isdiag)
    pairs = [(face_tps[f][1], face_tps[f][2]) for f in diag_faces]
    allunique(pairs) || throw(NonGenericGraph(
        "two strip regions share the same turning-point pair"))
    sort!(diag_faces; by = f -> (face_tps[f][1], face_tps[f][2]))
    bdry_faces = sort!(findall(!, isdiag); by = f -> endpoint_pair[f])
    new_id = zeros(Int, length(faces))
    for (i, f) in enumerate(diag_faces)
        new_id[f] = i
    end
    for (i, f) in enumerate(bdry_faces)
        new_id[f] = n - 1 + i
    end
    edge_endpoints = [endpoint_pair[f] for f in vcat(diag_faces, bdry_faces)]
    is_diagonal = vcat(trues(n - 1), falses(n + 2))
    diagonal_tp_pair = [(face_tps[f][1], face_tps[f][2]) for f in diag_faces]
    triangles = Vector{NTuple{3,Int}}(undef, n)
    for s in 1:n
        e = map(f -> new_id[f], tri_faces[s])
        shift = argmin(e) - 1                   # rotate to start at the smallest edge
        triangles[s] = ntuple(i -> e[mod1(i + shift, 3)], 3)
    end

    IdealTriangulation(n + 2, edge_endpoints, is_diagonal, triangles, collect(1:n),
                       diagonal_tp_pair, marked_angles, Float64(g.theta))
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

# the marked point shared by two edges (the corner they meet at)
function _shared_corner(t::IdealTriangulation, e1::Int, e2::Int)
    a, b = t.edge_endpoints[e1]
    c, d = t.edge_endpoints[e2]
    a == c || a == d ? a : (b == c || b == d ? b :
        throw(NonGenericGraph("edges $e1 and $e2 of a triangle share no corner")))
end

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
    u, v = t.edge_endpoints[k]

    # the apex of each adjacent triangle: the corner where its two non-`k` edges meet
    others(s) = Tuple(e for e in t.triangles[s] if e != k)
    o1, o2 = others(s1)
    o3, o4 = others(s2)
    p = _shared_corner(t, o1, o2)
    q = _shared_corner(t, o3, o4)

    # orient the quad counterclockwise as (u, p, v, q): marked points are labelled
    # counterclockwise, so the apex lying strictly between u and v comes second.
    if !(u < p < v)
        s1, s2 = s2, s1
        p, q = q, p
        (o1, o2), (o3, o4) = (o3, o4), (o1, o2)
    end

    # name the four quad sides by the corners they join
    has(e, x) = x in t.edge_endpoints[e]
    e_up = has(o1, u) ? o1 : o2          # u - p
    e_pv = e_up == o1 ? o2 : o1          # p - v
    e_vq = has(o3, v) ? o3 : o4          # v - q
    e_qu = e_vq == o3 ? o4 : o3          # q - u

    endpoints = copy(t.edge_endpoints)
    endpoints[k] = (min(p, q), max(p, q))

    # the two new triangles, counterclockwise: (u, p, q) and (p, v, q)
    rotate(e) = (shift = argmin(collect(e)) - 1;
                 ntuple(i -> e[mod1(i + shift, 3)], 3))
    triangles = copy(t.triangles)
    # `triangle_tp` is unchanged - each slot keeps its turning point - so the handedness
    # of the reconnection is expressed by WHICH slot each new triangle lands in. Crossing
    # the wall towards increasing θ sends the u-side triangle to the slot whose apex was
    # p; crossing the other way sends it to the slot whose apex was q. Pinned in both
    # senses by `canonical_reorder(flip(t, k; direction)) == t′` against a re-traced
    # neighbouring chamber, at every wall of the cubic/quartic/quintic fixtures.
    a_slot, b_slot = direction == 1 ? (s1, s2) : (s2, s1)
    triangles[a_slot] = rotate((e_up, k, e_qu))
    triangles[b_slot] = rotate((e_pv, e_vq, k))
    tp_pair = _diagonal_tp_pairs(endpoints, t.is_diagonal, triangles, t.triangle_tp)

    IdealTriangulation(t.n_marked, endpoints, t.is_diagonal, triangles, t.triangle_tp,
                       tp_pair, t.marked_angles, t.theta)
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
"""
function canonical_reorder(t::IdealTriangulation)
    diag = findall(t.is_diagonal)
    pair_of = Dict(e => p for (e, p) in zip(diag, t.diagonal_tp_pair))
    diag = sort(diag; by = e -> pair_of[e])
    bdry = sort(findall(!, t.is_diagonal); by = e -> t.edge_endpoints[e])
    perm = vcat(diag, bdry)
    new_id = zeros(Int, length(perm))
    for (i, e) in enumerate(perm)
        new_id[e] = i
    end
    triangles = map(t.triangles) do tri
        e = map(x -> new_id[x], collect(tri))
        shift = argmin(e) - 1
        ntuple(i -> e[mod1(i + shift, 3)], 3)
    end
    tri_out = IdealTriangulation(t.n_marked, t.edge_endpoints[perm], t.is_diagonal,
                                 triangles, t.triangle_tp,
                                 [pair_of[e] for e in diag], t.marked_angles, t.theta)
    tri_out, perm
end
