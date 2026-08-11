# The dual of a DEGENERATE Stokes graph, and the triangulations that refine it.
#
# An order-`m` turning point emits `m + 2` rays (`src/stokes_graph.jl`), so its dual
# cell is an `(m+2)`-gon rather than a triangle and the dual of the whole graph is an
# *ideal polygon decomposition* of the marked surface. In cluster terms that is a
# positive-dimensional face of the associahedron: the `Cat(m)` triangulations of the
# `(m+2)`-gon are its vertices, and they are exactly the ideal triangulations refining
# it. A degenerate graph is therefore a **wall**, and the simplest case is the one that
# names it - a *double* turning point is a square cell with two refinements, and the two
# differ by a single flip, so **a double turning point is a flip wall** ([IN14] §3,
# [BNR82] for the local geometry).
#
# The face walk that builds this is the one in `src/triangulation.jl`, written in terms
# of the valences `m + 2` rather than the constant 3; there is no second code path.
#
# ── Convention ledger ──────────────────────────────────────────────────────────────
#  • CELL ORIENTATION AND CORNER ALIGNMENT are inherited unchanged from
#    `IdealTriangulation`: `cells[s]` lists cell `s`'s edges counterclockwise starting
#    at its smallest edge index, aligned with `cell_corners[s]` so that edge
#    `cells[s][i]` runs from corner `i` to corner `i+1` (cyclically). Pinned by: the
#    refinements of a decomposition reproduce the freshly traced triangulations on
#    either side of the wall, which reads both fields.
#  • REFINEMENT ORDER. `refinements` returns the product over cells in cell order, and
#    within a cell the recursion of `_polygon_triangulations` (which fixes the base
#    edge `(1, v)` and splits at the apex, ascending). New diagonals are appended after
#    the decomposition's own edges and the result is then put in the canonical order of
#    `ideal_triangulation` by `canonical_reorder`, so the *set* is canonical even
#    though the enumeration order is ours. Pinned by: the two square refinements of the
#    critical-energy double well come out as the E ∓ δ triangulations, in that order.

"""
    PolygonDecomposition

The ideal polygon decomposition of the marked surface dual to a Stokes graph (build
with [`polygon_decomposition`](@ref)). Identical in structure to
[`IdealTriangulation`](@ref) - same marked points, same canonical edge order (diagonals
first), same `theta` provenance - except that `cells[s]`, the cell dual to turning
point `cell_tp[s]`, is an `(m+2)`-gon for a turning point of order `m` instead of a
triangle, listed counterclockwise and aligned with `cell_corners[s]` so that edge
`cells[s][i]` runs from corner `i` to corner `i+1`.

An all-simple graph gives an all-triangle decomposition, which is the
[`IdealTriangulation`](@ref) that `ideal_triangulation` returns. A degenerate one does
not: use [`refinements`](@ref) for the triangulations refining it.
"""
struct PolygonDecomposition
    n_marked::Int
    edge_endpoints::Vector{Tuple{Int,Int}}
    is_diagonal::Vector{Bool}
    cells::Vector{Vector{Int}}
    cell_tp::Vector{Int}
    diagonal_tp_pair::Vector{Tuple{Int,Int}}
    marked_angles::Vector{Float64}
    theta::Float64
    marked_boundary::Vector{Int}
    cell_corners::Vector{Vector{Int}}
    marked_is_puncture::Vector{Bool}
end

PolygonDecomposition(nm, ee, isd, cells, ctp, dtp, ma, θ, mb, cc) =
    PolygonDecomposition(nm, ee, isd, cells, ctp, dtp, ma, θ, mb, cc, falses(nm))

"""
    polygon_decomposition(g::StokesGraph) -> PolygonDecomposition

The ideal polygon decomposition dual to the Stokes graph `g`, degenerate or not. Same
genericity requirements as [`ideal_triangulation`](@ref) (every Stokes line must reach
a boundary circle) - and the same [`NonGenericGraph`](@ref) when a combinatorial
invariant fails - but turning points of any order are accepted.
"""
polygon_decomposition(g::StokesGraph) = _build_decomposition(g)[1]

"""
    cells(d::PolygonDecomposition) -> Vector{Vector{Int}}

The cells, one per turning point: `cells(d)[s]` lists the edges of the cell dual to
turning point `cell_tp(d)[s]`, counterclockwise.
"""
cells(d::PolygonDecomposition) = d.cells

"""
    cell_tp(d::PolygonDecomposition) -> Vector{Int}

The turning point each cell is dual to (index into the Stokes graph's
`turning_points`).
"""
cell_tp(d::PolygonDecomposition) = d.cell_tp

"""
    cell_sizes(d::PolygonDecomposition) -> Vector{Int}

The number of sides of each cell: `m + 2` at a turning point of order `m`, so `3` for
every cell of a triangulation.
"""
cell_sizes(d::PolygonDecomposition) = length.(d.cells)

"""
    n_cells(d::PolygonDecomposition) -> Int

Number of cells (= number of turning points).
"""
n_cells(d::PolygonDecomposition) = length(d.cells)

n_marked_points(d::PolygonDecomposition) = d.n_marked
n_boundaries(d::PolygonDecomposition) = maximum(d.marked_boundary; init = 1)
n_punctures(d::PolygonDecomposition) = count(d.marked_is_puncture)
n_diagonals(d::PolygonDecomposition) = count(d.is_diagonal)
diagonals(d::PolygonDecomposition) = findall(d.is_diagonal)

"""
    is_triangulation(d::PolygonDecomposition) -> Bool

`true` when every cell is a triangle - i.e. the Stokes graph was not degenerate, and
`d` is an [`IdealTriangulation`](@ref) in disguise.
"""
is_triangulation(d::PolygonDecomposition) = all(c -> length(c) == 3, d.cells)

# -- refinements ---------------------------------------------------------------------

# Catalan numbers, for `n_refinements`.
_catalan(n::Integer) = n ≤ 0 ? 1 : binomial(2n, n) ÷ (n + 1)

"""
    n_refinements(d::PolygonDecomposition) -> Int

How many ideal triangulations refine `d`: `∏ Cat(v − 2)` over the cells, `v` the
number of sides. `1` for a triangulation, `2` for a single double turning point (the
flip wall), `5` for a single triple one.
"""
n_refinements(d::PolygonDecomposition) = prod(_catalan(length(c) - 2) for c in d.cells;
                                              init = 1)

# All triangulations of the convex polygon on corner positions `lo … hi` whose base is
# the chord (lo, hi), as lists of ccw position triples. `Cat(hi − lo − 1)` of them.
function _polygon_triangulations(lo::Int, hi::Int)
    hi - lo < 2 && return [NTuple{3,Int}[]]
    out = Vector{Vector{NTuple{3,Int}}}()
    for k in (lo + 1):(hi - 1), A in _polygon_triangulations(lo, k),
        B in _polygon_triangulations(k, hi)

        push!(out, vcat(A, B, [(lo, k, hi)]))
    end
    out
end

_polygon_triangulations(v::Int) = _polygon_triangulations(1, v)

"""
    refinements(d::PolygonDecomposition) -> Vector{IdealTriangulation}

Every ideal triangulation refining `d`: each `(m+2)`-gon cell is triangulated in each
of its `Cat(m)` ways, independently, so there are [`n_refinements`](@ref) of them, all
in the canonical edge order of [`ideal_triangulation`](@ref). A triangulation refines
itself, so a non-degenerate decomposition gives the one-element vector.

This is the combinatorial content of "a degenerate Stokes graph is a wall": for a
single **double** turning point the square cell has two triangulations, they differ by
one [`flip`](@ref), and they are the triangulations traced on either side of the wall -
so [`triangulation_quiver`](@ref) of the two are related by mutation at the new
diagonal.

The new diagonals carry `diagonal_tp_pair = (s, s)`: both triangles beside a diagonal
internal to a cell are dual to the *same* turning point, which is exactly why the pair
is not a name for a cycle (see [`diagonal_core_paths`](@ref)). For the same reason
`triangle_tp` repeats a degenerate turning point once per triangle of its cell, and the
triangles are listed cell by cell - a different slot order from a freshly traced
triangulation, whose edges and quiver nevertheless agree.
"""
function refinements(d::PolygonDecomposition)
    per_cell = [_polygon_triangulations(length(c)) for c in d.cells]
    out = IdealTriangulation[]
    for choice in Iterators.product(per_cell...)
        push!(out, _refine(d, collect(choice)))
    end
    out
end

# Build the triangulation refining `d` that cuts cell `s` along `choice[s]` (position
# triples into `cells[s]` / `cell_corners[s]`).
function _refine(d::PolygonDecomposition, choice::Vector{Vector{NTuple{3,Int}}})
    endpoints = copy(d.edge_endpoints)
    is_diag = copy(d.is_diagonal)
    triangles = NTuple{3,Int}[]
    corners = NTuple{3,Int}[]
    triangle_tp = Int[]
    for (s, tris) in enumerate(choice)
        edges, cs = d.cells[s], d.cell_corners[s]
        v = length(edges)
        # side (p, q) of the cell polygon: an existing edge when q = p + 1 (cyclically),
        # otherwise a new diagonal, minted once and shared by its two triangles.
        new_edge = Dict{Tuple{Int,Int},Int}()
        function side(p, q)
            mod1(p + 1, v) == q && return edges[p]
            key = (min(p, q), max(p, q))
            haskey(new_edge, key) && return new_edge[key]
            push!(endpoints, (min(cs[p], cs[q]), max(cs[p], cs[q])))
            push!(is_diag, true)
            new_edge[key] = length(endpoints)
        end
        for (i, j, k) in tris
            e = (side(i, j), side(j, k), side(k, i))
            c = (cs[i], cs[j], cs[k])
            shift = _canonical_rotation(e, c)   # the canonical rotation of the walk
            push!(triangles, ntuple(x -> e[mod1(x + shift, 3)], 3))
            push!(corners, ntuple(x -> c[mod1(x + shift, 3)], 3))
            push!(triangle_tp, d.cell_tp[s])
        end
    end
    tp_pair = _diagonal_tp_pairs(endpoints, is_diag, triangles, triangle_tp)
    t = IdealTriangulation(d.n_marked, endpoints, is_diag, triangles, triangle_tp,
                           tp_pair, d.marked_angles, d.theta, d.marked_boundary,
                           corners, d.marked_is_puncture)
    first(canonical_reorder(t))
end
