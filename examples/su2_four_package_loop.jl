# The four-package loop: pure SU(2), computed twice by disjoint routes and compared.
#
# This is the only place the four packages meet in one program, and it is why it is a
# script rather than a test: ExactWKB.jl does not depend on ClusterSurfaces.jl and
# must not start (it is unregistered, so a test dependency would break CI). Each leg
# is asserted inside the suite of the package that owns it; this script runs the whole
# circle at once.
#
#   Resurgence.jl / ExactWKB.jl   Mathieu potential on the w = e^{2ix} plane
#          ↓ stokes_graph            two order-3 poles ⇒ the annulus
#          ↓ ideal_triangulation     one marked point on each boundary circle
#   ClusterAlgebras.jl            triangulation_quiver = Kronecker = Ã(1,1)
#          ↓ mutation_type
#   ClusterSurfaces.jl            surface_realization sends it back to Annulus(1,1)
#          ↓
#   ExactWKB.jl / sw_curve.jl     the SAME theory from the Seiberg–Witten curve:
#                                 the two classical periods must agree with a, a_D
#
# The last step is the one that had never been done. `sw_curve.jl` is firewalled from
# the WKB engine - it reaches a and a_D through Carlson elliptic integrals on the SW
# curve - while this route traces a Schrödinger problem's Stokes geometry and does a
# branch-tracked quadrature. Nothing is shared but the physics.
#
# Run:  julia examples/su2_four_package_loop.jl
#
# It provisions its own temporary environment, because ClusterSurfaces.jl is
# unregistered and must not appear in `test/Project.toml`.

using Pkg
Pkg.activate(temp = true)
let here = @__DIR__
    Pkg.develop([PackageSpec(path = joinpath(here, "..")),
                 PackageSpec(path = joinpath(here, "..", "..", "ClusterSurfaces.jl"))])
end

using ExactWKB
import ClusterAlgebras
import ClusterSurfaces

Λ, u = 1.0, 3.0

# ── 1. the potential ────────────────────────────────────────────────────────────────
prob = mathieu_problem(Λ, u)
show(stdout, MIME"text/plain"(), prob); println("\n")
println("turning points : ", location.(turning_points(prob)))
println("   (closed form: (u ± √(u²−4Λ⁴))/2Λ² = ",
        ((u - sqrt(u^2 - 4Λ^4)) / (2Λ^2), (u + sqrt(u^2 - 4Λ^4)) / (2Λ^2)), ")")

# ── 2. Stokes graph on the annulus ──────────────────────────────────────────────────
g = stokes_graph(prob; theta = 0.3)
println("\n", g)
println("rays to ∞      : ", count(l -> ExactWKB.endpoint(l) === :infinity, ExactWKB.lines(g)))
println("rays to w = 0  : ", count(l -> ExactWKB.endpoint(l) === :pole, ExactWKB.lines(g)))

t = ideal_triangulation(g)
println("\nboundary circles : ", ExactWKB.n_boundaries(t))
println("marked points    : ", n_marked_points(t), " on circles ", t.marked_boundary)
println("arcs             : ", t.edge_endpoints[1:n_diagonals(t)])
println("boundary segments: ", t.edge_endpoints[(n_diagonals(t) + 1):end])

# ── 3. the quiver, and back to the surface ──────────────────────────────────────────
B = triangulation_quiver(t).B
println("\nquiver B         : ", B)
println("mutation type    : ", ClusterAlgebras.mutation_type(ClusterAlgebras.Quiver(B)))
println("= SU(2) BPS quiver: ",
        ClusterAlgebras.canonical_form(ClusterAlgebras.Quiver(B)).B ==
        ClusterAlgebras.canonical_form(su2_bps_quiver()).B)

surf = ClusterSurfaces.surface_realization(ClusterAlgebras.Quiver(B))
println("surface_realization: genus ", ClusterSurfaces.genus(surf),
        ", ", ClusterSurfaces.n_boundary_components(surf), " boundary components, ",
        ClusterSurfaces.marked_points(surf), " marked points")

# ── 4. the periods, against the Seiberg–Witten curve ────────────────────────────────
# electric = the unit circle |w| = 1; magnetic = the cycle around the turning points.
zs = location.(turning_points(prob))
circle = [cis(2π * j / 512) for j in 0:511]
m = 0.4 * min(abs(zs[1]), abs(zs[2] - zs[1]) / 2)
a_wkb = period_integral(prob, circle) / (im * π)
aD_wkb = -period_integral(prob, encircling_contour(zs[1], zs[2]; margin = m, n = 256)) / (im * π)
per = sw_periods(SeibergWittenSU2(Λ = Λ), u)

row(name, x, y) = println(rpad(name, 7), rpad(round(x; sigdigits = 15), 26),
                          rpad(round(y; sigdigits = 15), 26), abs(x - y))
println("\n", rpad("", 7), rpad("WKB route", 34), rpad("SW-curve route", 26), "|Δ|")
row("a", a_wkb, complex(per.a))
row("a_D", aD_wkb, complex(per.a_D))

# ── 4b. the charge lattice, built from the strip regions ────────────────────────────
# Both diagonals join turning points 1 and 2, so an ellipse around the pair cannot name
# either cycle; `charge_basis(prob, g)` builds each from its own strip instead.
cb = charge_basis(prob, g)
println("\ndiagonal turning-point pairs : ", t.diagonal_tp_pair, "  (both the same!)")
println("measured pairing             : ", signed_pairing(cb))
println("keystone −B                  : ", -B)
for (e, Z) in enumerate(physical_charges(cb))
    println("  γ$e : Z = ", round(Z; sigdigits = 12))
end
println("  iπ·a_D (monopole)          = ", round(im * π * per.a_D; sigdigits = 12))
println("  iπ(−a_D + 2a) (dyon)       = ", round(im * π * (-per.a_D + 2per.a); sigdigits = 12))

# ── 5. the saddle is the monopole ───────────────────────────────────────────────────
s = only(saddle_candidates(prob))
println("\nsaddle           : ", s)
println("|Z| vs π|a_D|    : ", ExactWKB.mass(s), " vs ", π * abs(per.a_D))
println("\nstrong-coupling BPS states from the curve side:")
for st in su2_bps_states(SeibergWittenSU2(Λ = Λ), 1.0; chamber = :strong)
    println("  ", st)
end
