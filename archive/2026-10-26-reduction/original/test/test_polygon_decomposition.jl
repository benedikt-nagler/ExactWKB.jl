# The degenerate half of the dual construction: an order-m turning point is dual to an
# (m+2)-gon, so a degenerate Stokes graph is dual to a polygon decomposition and the
# Cat(m) triangulations of that cell are the refinements. The headline oracle is the
# name of the thing - **a double turning point is a flip wall**: the square cell's two
# refinements are exactly the triangulations traced on either side of the critical
# energy, they differ by one flip, and their quivers by one mutation.

import ClusterAlgebras

@testset "polygon decomposition (degenerate turning points)" begin

    @testset "a non-degenerate graph decomposes into triangles" begin
        # the same walk, so the decomposition must agree field for field with the
        # triangulation it is the all-simple case of
        prob = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])          # cubic, A₂
        g = stokes_graph(prob; theta = 0.3)
        d = polygon_decomposition(g)
        t = ideal_triangulation(g)
        @test is_triangulation(d)
        @test cell_sizes(d) == [3, 3, 3]
        @test n_cells(d) == 3
        @test cell_tp(d) == collect(1:3)
        @test n_diagonals(d) == n_diagonals(t) == 2
        @test diagonals(d) == diagonals(t)
        @test n_marked_points(d) == n_marked_points(t) == 5
        @test n_boundaries(d) == 1
        @test d.edge_endpoints == t.edge_endpoints
        @test cells(d) == [collect(tri) for tri in t.triangles]
        @test d.cell_corners == [collect(c) for c in t.triangle_corners]
        @test n_refinements(d) == 1
        r = refinements(d)
        @test length(r) == 1
        @test r[1].edge_endpoints == t.edge_endpoints
        @test r[1].triangles == t.triangles
        @test r[1].triangle_corners == t.triangle_corners
        @test triangulation_quiver(r[1]).B == triangulation_quiver(t).B
    end

    @testset "Q = z²: one square cell, two refinements" begin
        # A double turning point emits four rays and IS the whole graph: the dual is a
        # single square, the disk's own quadrilateral, with no diagonal at all.
        prob = SchrodingerProblem([0.0, 0.0, 1.0])
        g = stokes_graph(prob; theta = 0.3)
        @test is_degenerate(g)
        @test turning_point_orders(g) == [2]
        @test length(ExactWKB.lines(g)) == 4
        d = polygon_decomposition(g)
        @test !is_triangulation(d)
        @test cell_sizes(d) == [4]
        @test n_marked_points(d) == 4          # degree 2 ⇒ e + 2 = 4 directions
        @test n_diagonals(d) == 0              # d = n − 2 + b = 1 − 2 + 1
        @test length(d.edge_endpoints) == 4
        @test n_refinements(d) == 2            # Cat(2)
        r = refinements(d)
        @test length(r) == 2
        # the two diagonals of the square, and nothing else changed
        @test [t.edge_endpoints[1] for t in r] == [(2, 4), (1, 3)]
        @test all(t -> n_diagonals(t) == 1, r)
        @test all(t -> length(t.triangles) == 2, r)
        @test all(t -> sort(t.edge_endpoints[2:end]) ==
                       [(1, 2), (1, 4), (2, 3), (3, 4)], r)
        # rank-1 quiver: a single diagonal, no arrows
        @test all(t -> triangulation_quiver(t).B == zeros(Int, 1, 1), r)
        # one flip apart - the wall in its barest form
        @test first(canonical_reorder(flip(r[1], 1))).edge_endpoints ==
              r[2].edge_endpoints
    end

    @testset "Q = z³: one pentagon cell, five refinements" begin
        # Cat(3) = 5, and the divisor identity Σ orders = M + 2b − 4 with a single
        # order-3 zero: 3 = 5 + 2 − 4.
        prob = SchrodingerProblem([0.0, 0.0, 0.0, 1.0])
        g = stokes_graph(prob; theta = 0.3)
        @test turning_point_orders(g) == [3]
        @test length(ExactWKB.lines(g)) == 5
        d = polygon_decomposition(g)
        @test cell_sizes(d) == [5]
        @test n_marked_points(d) == 5
        @test n_diagonals(d) == 0
        @test n_refinements(d) == 5
        r = refinements(d)
        @test length(r) == 5
        @test all(t -> n_diagonals(t) == 2, r)
        @test all(t -> length(t.triangles) == 3, r)
        # the five triangulations of the pentagon are distinct, and each is A₂
        @test length(unique(sort(t.edge_endpoints[1:2]) for t in r)) == 5
        @test all(t -> abs.(triangulation_quiver(t).B) == [0 1; 1 0], r)
    end

    @testset "a double turning point is a flip wall" begin
        # The quartic double well at the barrier top E = 1: Q = z²(z² − 2), so the
        # merging pair has become one order-2 point. Cells (3, 4, 3), and the square's
        # two triangulations are the two chambers of the E-wall.
        dw = SchrodingerProblem([1.0, 0.0, -2.0, 0.0, 1.0])
        gc = stokes_graph(with_energy(dw, 1.0); theta = 0.3)
        @test turning_point_orders(gc) == [1, 1, 2]
        d = polygon_decomposition(gc)
        @test cell_sizes(d) == [3, 4, 3]        # cell 2 is the double turning point
        @test n_cells(d) == 3
        @test n_marked_points(d) == 6           # quartic ⇒ 6, and 4 = 6 + 2 − 4
        @test n_diagonals(d) == 2               # n − 2 + b = 3 − 2 + 1
        @test n_refinements(d) == 2
        r = refinements(d)
        # generic energies either side of the wall, at the same phase
        tm = ideal_triangulation(stokes_graph(with_energy(dw, 0.98); theta = 0.3))
        tp = ideal_triangulation(stokes_graph(with_energy(dw, 1.02); theta = 0.3))
        @test all(t -> n_diagonals(t) == 3, (r[1], r[2], tm, tp))   # A₃ either side
        # THE ORACLE: the refinements are the traced triangulations, in order
        @test r[1].edge_endpoints == tm.edge_endpoints
        @test r[2].edge_endpoints == tp.edge_endpoints
        # as SETS: a refinement lists its triangles cell by cell (so `triangle_tp`
        # repeats the double point), which is a different slot order from a traced
        # triangulation's one-triangle-per-turning-point order
        @test sort(collect(r[1].triangles)) == sort(collect(tm.triangles))
        @test sort(collect(r[2].triangles)) == sort(collect(tp.triangles))
        @test r[1].triangle_tp == r[2].triangle_tp == [1, 2, 2, 3]
        @test triangulation_quiver(r[1]).B == triangulation_quiver(tm).B
        @test triangulation_quiver(r[2]).B == triangulation_quiver(tp).B
        # ... and the wall acts as one mutation, at the diagonal internal to the cell
        k = findfirst(p -> p[1] == p[2], r[1].diagonal_tp_pair)::Int
        @test r[1].diagonal_tp_pair[k] == (2, 2)   # both sides dual to the double point
        Bm = triangulation_quiver(tm)
        @test ClusterAlgebras.mutate(Bm, k).B == triangulation_quiver(tp).B
        @test first(canonical_reorder(flip(r[1], k))).edge_endpoints ==
              r[2].edge_endpoints
        # the other two diagonals are untouched by the wall
        @test [r[1].edge_endpoints[e] for e in 1:3 if e != k] ==
              [r[2].edge_endpoints[e] for e in 1:3 if e != k]
    end

    @testset "the triangulation layer still refuses a degenerate graph" begin
        gc = stokes_graph(with_energy(SchrodingerProblem([1.0, 0.0, -2.0, 0.0, 1.0]),
                                      1.0); theta = 0.3)
        @test_throws NonGenericGraph ideal_triangulation(gc)
        err = try
            ideal_triangulation(gc)
        catch e
            e
        end
        @test occursin("polygon decomposition", sprint(showerror, err))
        # the charge lattice is a triangulation consumer, so it refuses too
        @test_throws NonGenericGraph charge_basis(
            with_energy(SchrodingerProblem([1.0, 0.0, -2.0, 0.0, 1.0]), 1.0), gc)
        # ... but the decomposition and its strip core paths are available
        @test length(diagonal_core_paths(gc)) == 2
    end
end
