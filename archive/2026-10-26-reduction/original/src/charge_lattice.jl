# The numerical half of the cluster bridge: one closed contour per diagonal of the
# ideal triangulation - the charge lattice basis γ₁ … γ_{d−1}. Each γ_j encircles the
# turning-point pair of its diagonal, is oriented by the DDP decay rule (ledger item 4
# of src/ddp.jl: `Re(Z_j e^{-iθ_c}) < 0` along its own wall, i.e. `Z_j` in the closed
# lower half-plane minus the positive real axis), and carries the central charge
# `Z_j = ∮√Q` and the intersection pairing `P[i,j] = ⟨γ_i, γ_j⟩`.
#
# In physics language: Γ = ⊕ ℤγ_j is the electromagnetic charge lattice of the
# associated N = 2 Coulomb branch, ⟨,⟩ its Dirac pairing, and Z : Γ → ℂ the N = 2
# central charge - computed here literally as WKB period integrals.
#
# The decay rule fixes a *representative* of ±γ_e, not the physical charge. The
# orientation signs ε that promote it to one - `γ_e^phys = ε_e·γ_e` - are the signed
# frame of `src/signed_frame.jl`, solved here at construction. So the keystone
# consistency of the bridge is the STRICT identity
#
#     P[i,j] == −ε_i ε_j·B[i,j],     B = triangulation_quiver(t).B,
#
# i.e. `signed_pairing(cb) == −B` exactly, and `B_seed = −signed_pairing` (what
# `ddp_seed` builds, ledger item 5) is the triangulation quiver on the nose.
#
# ── Convention ledger: the punctured-surface route (ours, not forced) ──────────────
# `charge_basis(prob, t)` names a cycle by its diagonal's turning-point PAIR and puts
# an ellipse around it. That stops being a name as soon as two diagonals share a pair -
# on the annulus both arcs of `Ã(1,1)` do, and they differ by the loop around the hole,
# which is the whole difference between the monopole and the dyon. `charge_basis(prob,
# g)` therefore builds each cycle from its own strip region instead.
#  • WHICH SIDE OF THE STRIP. A strip is a quadrilateral `z_s → M₁ → z_{s'} → M₂` with
#    the turning points at opposite corners, so either side is a core path. We take the
#    cheaper in accumulated Stokes `mass`: the two legs contribute ±e^{iθ}t to the
#    period and partly cancel, so the excursion is what the integral costs. Pinned by:
#    the annulus reproducing `sw_curve.jl`'s a and a_D to 1e-11, and the polynomial
#    fixtures agreeing with the ellipse route to 1e-9.
#  • DISTINCT RACETRACK WIDTHS. Two strips routinely share a boundary ray - on the
#    annulus both do - and equal-width racetracks around it are *collinear*, which
#    `intersection_pairing` refuses (rightly: the crossing is not transversal). Widths
#    are therefore scaled by 1/(1 + 0.13·e), strictly decreasing, so no two cycles ever
#    run collinearly. Legitimate because the pairing is a homological invariant of the
#    lifts, so perturbing a representative cannot change it. Pinned by: the graph route
#    reproducing the ellipse route's P and ε *exactly* on every polynomial fixture.
#  • WINDING NUMBERS DO NOT PIN THIS. They are homology in the punctured plane and
#    cannot see which way a core path wraps; both annulus cycles wind (1,1,0) about
#    (z₁,z₂,pole). What pins the construction is the keystone below.
#
# ── this closed the original-bridge finding (2026-07-19) ───────────────────────────────────
# The original bridge layer checked only `|P| == |B|`, and that `abs` was exactly the slack hiding the
# missing layer. With the uniform all-decay frame the seed structure survived only on
# *tree chambers*; on a *cyclic chamber* (a triangle of the triangulation with all
# three edges internal, generic for degree ≥ 5) the decay-oriented pairing orients
# that quiver triangle acyclically, so `−P(decay)` is cluster-INFINITE and cannot be
# a seed - `charge_basis` used to refuse outright. Solving for ε repairs it: the
# cocycle condition around the quiver 3-cycle is a genuine constraint and it closes
# in every chamber measured (see test_signed_frame.jl). `−P(decay)` remains
# cluster-infinite there - that fact is kept as a regression test - but it is no
# longer the seed.

"""
    ChargeBasis{F}

The charge-lattice basis attached to an [`IdealTriangulation`](@ref): per diagonal a
closed `contour` (decay-oriented, ledger item 4), its central charge
`Z_j = ∮√Q ∈ central_charges`, the integer intersection-pairing matrix
`pairing[i,j] = ⟨γ_i, γ_j⟩`, and the orientation `signs` of the signed frame
(`γ_j^phys = signs[j] · γ_j`, see `src/signed_frame.jl`). Build with
[`charge_basis`](@ref); feed to [`bridge_seed`](@ref) / [`voros_symbol`](@ref).

`contours` and `central_charges` are the **decay-frame representatives** and are never
reoriented - ledger item 4 stays literal, and so does the spectral orientation ledger. The
physical objects are [`physical_charges`](@ref) and [`signed_pairing`](@ref).
"""
struct ChargeBasis{F}
    problem::AbstractSchrodingerProblem
    triangulation::IdealTriangulation
    contours::Vector{Vector{Complex{F}}}
    central_charges::Vector{Complex{F}}
    pairing::Matrix{Int}
    signs::Vector{Int}
end

"""
    n_charges(cb::ChargeBasis) -> Int

The rank of the charge lattice: one basis cycle per diagonal of the triangulation.
"""
n_charges(cb::ChargeBasis) = length(cb.contours)
central_charges(cb::ChargeBasis) = cb.central_charges

"""
    signs(cb::ChargeBasis) -> Vector{Int}

The orientation signs `ε` of the signed frame: `γ_j^phys = ε_j · γ_j`. See
[`signed_frame`](@ref).
"""
signs(cb::ChargeBasis) = cb.signs

"""
    physical_charges(cb::ChargeBasis) -> Vector{Complex}

The central charges of the *physical* cycles, `Z_j^phys = ε_j · Z_j`. These carry the
BPS masses and wall phases; `central_charges` gives the decay-frame representatives.
"""
physical_charges(cb::ChargeBasis) = cb.signs .* cb.central_charges

"""
    signed_pairing(cb::ChargeBasis) -> Matrix{Int}

The intersection pairing of the physical cycles, `ε_i ε_j ⟨γ_i, γ_j⟩`. This is the `P`
of ledger item 5, so the seed is `B = −signed_pairing(cb)` ([`bridge_seed`](@ref)),
and it equals `−triangulation_quiver(t).B` exactly (the tightened keystone).
"""
signed_pairing(cb::ChargeBasis) = cb.signs .* cb.pairing .* cb.signs'

# winding number of a closed polygonal contour about z
function _winding(verts::AbstractVector{Complex{F}}, z) where {F}
    s = zero(F)
    for k in eachindex(verts)
        a = verts[k] - z
        b = verts[mod1(k + 1, length(verts))] - z
        s += angle(b / a)
    end
    round(Int, s / (2 * F(π)))
end

# distance from z to the segment [a, b]
function _seg_dist(z, a, b)
    t = clamp(real((z - a) * conj(b - a)) / abs2(b - a), 0, 1)
    abs(z - (a + t * (b - a)))
end

"""
    charge_contour(prob::SchrodingerProblem, i, j; margin = nothing, n = 32) -> Vector

A closed contour encircling exactly the simple turning points `i` and `j` of `prob`
(indices into `simple_turning_points(prob)`): an [`encircling_contour`](@ref) ellipse
whose `margin` defaults to `min(sep/4, 0.4·dist)` with `sep` the pair separation and
`dist` the distance of the nearest third turning point to the segment. The enclosure
is validated by winding numbers (1 around `i` and `j`, 0 around every other turning
point); on failure the margin is halved once, then [`ContourError`](@ref) is thrown.
"""
function charge_contour(prob::SchrodingerProblem, i::Integer, j::Integer;
                        margin = nothing, n::Integer = 32)
    tps = simple_turning_points(prob)
    m = length(tps)
    (1 ≤ i ≤ m && 1 ≤ j ≤ m && i != j) || throw(Resurgence.InvalidArgument(
        "need two distinct turning-point indices in 1:$m, got ($i, $j)"))
    zi, zj = location(tps[i]), location(tps[j])
    F = typeof(real(zi))
    sep = abs(zj - zi)
    others = [location(tps[k]) for k in 1:m if k != i && k != j]
    pad = if margin === nothing
        d = minimum((_seg_dist(z, zi, zj) for z in others); init = F(Inf))
        min(sep / 4, F(2) / 5 * d)
    else
        F(margin)
    end
    pad > 0 || throw(ContourError(
        "cannot enclose turning points ($i, $j) alone: another turning point sits " *
        "on their segment"))
    for _ in 1:2
        c = encircling_contour(zi, zj; margin = pad, n)
        wi, wj = _winding(c, zi), _winding(c, zj)
        wothers = [_winding(c, z) for z in others]
        (abs(wi) == 1 && wj == wi && all(iszero, wothers)) && return c
        pad /= 2
    end
    throw(ContourError(
        "no valid encircling contour for turning points ($i, $j): a third turning " *
        "point stays inside the ellipse - pass an explicit margin"))
end

# -- racetrack contours around an explicit core path --------------------------------
#
# On a punctured surface a turning-point pair is no longer a name for a cycle - two
# diagonals can share one and differ by a loop around a hole (the annulus: monopole vs
# dyon). The cycle is the double-cover lift of a path running through the diagonal's
# own strip region ([`diagonal_core_paths`](@ref)), so the contour is a thin closed
# loop around that path: offset by ±δ, capped by half-circles that pass on the far
# side of each turning point, which is what gives winding 1 about both.

# drop points closer than `tol` to the previous kept one (keeps the offset construction
# well-conditioned and the contour affordable; the endpoints are always kept)
function _decimate(pts::Vector{Complex{F}}, tol) where {F}
    out = Complex{F}[pts[1]]
    for z in @view pts[2:end]
        abs(z - out[end]) ≥ tol && push!(out, z)
    end
    length(out) == 1 || out[end] == pts[end] || push!(out, pts[end])
    out
end

function _racetrack(zs::Complex{F}, pts::Vector{Complex{F}}, ze::Complex{F},
                    δ::F, cap::Int) where {F}
    k = length(pts)
    nrm = Vector{Complex{F}}(undef, k)
    for i in 1:k
        a = i == 1 ? zs : pts[i - 1]
        b = i == k ? ze : pts[i + 1]
        nrm[i] = im * (b - a) / abs(b - a)
    end
    u = (pts[1] - zs) / abs(pts[1] - zs)
    v = (ze - pts[k]) / abs(ze - pts[k])
    # half-circle of radius δ about `c`, from angle `from` sweeping −π (interior points
    # only - the two ends are supplied by the offsets)
    halfcap(c, from) = [c + δ * cis(from - F(π) * s)
                        for s in range(zero(F), one(F), length = cap)[2:(end - 1)]]
    out = Complex{F}[zs + δ * im * u]
    append!(out, (pts[i] + δ * nrm[i] for i in 1:k))
    push!(out, ze + δ * im * v)
    append!(out, halfcap(ze, angle(v) + F(π) / 2))      # over the far side of ze
    push!(out, ze - δ * im * v)
    append!(out, (pts[i] - δ * nrm[i] for i in k:-1:1))
    push!(out, zs - δ * im * u)
    append!(out, halfcap(zs, angle(u) - F(π) / 2))      # behind zs
    out
end

# Split a core path into (all vertices, start turning point, end turning point,
# interior points, the singularities the racetrack must NOT enclose). The two turning
# points the path runs between are its endpoints, identified by proximity because the
# traced graph's float type need not be the problem's.
function _core_path_parts(prob::AbstractSchrodingerProblem, path::AbstractVector)
    length(path) ≥ 3 || throw(Resurgence.InvalidArgument(
        "charge_contour: a core path needs ≥ 3 vertices (both turning points and at " *
        "least one interior point), got $(length(path))"))
    verts = [_as_complex(z) for z in path]
    F = _contour_float(verts)
    zs, ze = verts[1], verts[end]
    pts = Complex{F}[verts[i] for i in 2:(length(verts) - 1)]
    tps = Complex{F}[Complex{F}(location(t)) for t in simple_turning_points(prob)]
    ends = isempty(tps) ? Int[] : [argmin(abs.(tps .- zs)), argmin(abs.(tps .- ze))]
    others = vcat(Complex{F}[tps[k] for k in eachindex(tps) if !(k in ends)],
                  Complex{F}[Complex{F}(p) for p in poles(prob)])
    verts, zs, ze, pts, others
end

# Default racetrack half-width: clear of every other singularity, and thin enough that
# the end caps sit inside the spacing of the path's own first and last steps.
function _path_margin(zs, ze, pts, verts, others)
    F = _contour_float(verts)
    clearance = minimum((abs(z - o) for z in verts, o in others); init = F(Inf))
    min(F(2) / 5 * clearance, abs(pts[1] - zs) / 2, abs(ze - pts[end]) / 2)
end

"""
    charge_contour(prob, path::AbstractVector; margin = nothing, cap = 8) -> Vector

A closed contour around the open core `path` - a thin racetrack whose two ends cap
around `path[1]` and `path[end]`, which must be the two turning points the cycle
encircles. Use it when a turning-point *pair* does not determine the cycle, i.e. on
any surface with a puncture; [`diagonal_core_paths`](@ref) supplies the paths.

`margin` is the half-width `δ` (default: a fraction of the path's clearance from every
other turning point and pole, and of its own end spacing). Validated by the same
winding-number test as the pair method - 1 about each end, 0 about every other turning
point *and* every pole - halving `δ` once before throwing [`ContourError`](@ref).

Winding numbers are homology in the punctured plane and cannot see which way the core
path wraps; the construction is pinned instead by `charge_basis`'s keystone.
"""
function charge_contour(prob::AbstractSchrodingerProblem, path::AbstractVector;
                        margin = nothing, cap::Integer = 8)
    verts, zs, ze, pts, others = _core_path_parts(prob, path)
    F = _contour_float(verts)
    δ = margin === nothing ? _path_margin(zs, ze, pts, verts, others) : F(margin)
    δ > 0 || throw(ContourError(
        "cannot thicken the core path from $zs to $ze: another turning point or pole " *
        "lies on it"))
    for _ in 1:2
        c = _racetrack(zs, _decimate(pts, δ), ze, δ, Int(cap))
        ws, we = _winding(c, zs), _winding(c, ze)
        if abs(ws) == 1 && we == ws && all(o -> iszero(_winding(c, o)), others)
            return c
        end
        δ /= 2
    end
    throw(ContourError(
        "no valid racetrack around the core path from $zs to $ze: it encloses another " *
        "turning point or pole - pass an explicit margin"))
end

# Decay-oriented contour + central charge: the cycle is reversed unless Z lies in the
# closed lower half-plane minus the positive real axis, so that Re(Z e^{-iθ_c}) < 0
# along the cycle's own wall (ledger item 4) - this is what makes every basis symbol a
# legal wall symbol for verify_ddp/verify_ddp_mutation.
function _decay_orient(c, Z)
    # numerically real Z (wall at θ_c = 0) must resolve by the real part, not by the
    # sign of a roundoff-sized imaginary residue
    tol = sqrt(eps(typeof(real(Z)))) * (1 + abs(Z))
    (imag(Z) > tol || (abs(imag(Z)) ≤ tol && real(Z) > 0)) ? (reverse(c), -Z) : (c, Z)
end

function _oriented_contour(prob, i, j; margin, n, rtol)
    c = charge_contour(prob, i, j; margin, n)
    _decay_orient(c, period_integral(prob, c; rtol))
end

function _oriented_path_contour(prob, path; margin, cap, rtol)
    c = charge_contour(prob, path; margin, cap)
    _decay_orient(c, period_integral(prob, c; rtol))
end

# The charge lattice does NOT yet run on a punctured surface, and the reason is
# specific: a cycle here is a homology class of the double cover, built either from an
# ellipse around a turning-point pair or (M8b) from the diagonal's own strip region.
# A puncture is a point the double cover branches around and a residue cycle nobody has
# put in the basis - `∮√Q = 2πi√c ≠ 0` around it - so neither construction is known to
# represent the right class. This is the deferred Tier B item in
# ../../PLANNING/horizon.md; refusing is the honest state, not a missing line of code.
function _require_unpunctured(t::IdealTriangulation)
    p = findfirst(t.marked_is_puncture)
    p === nothing && return nothing
    throw(ContourError(
        "the surface has a puncture (marked point $p): the charge lattice is not " *
        "built on a punctured surface, because a cycle would have to be told apart " *
        "from the residue cycle around the puncture. Use ideal_triangulation and " *
        "triangulation_quiver, which are"))
end

"""
    charge_basis(prob::SchrodingerProblem, t::IdealTriangulation;
                 margin = nothing, n = 32, rtol = nothing, verify = true) -> ChargeBasis

The charge-lattice basis of the triangulation `t`: for each diagonal (in `t`'s
canonical order) a closed contour around its turning-point pair
([`charge_contour`](@ref)), oriented to decay along its own wall (ledger item 4 -
`Re(Z_j e^{-iθ_c}) < 0`), with central charges `Z_j = ∮√Q` and the intersection
pairing `P[i,j] = ⟨γ_i, γ_j⟩` from [`intersection_pairing`](@ref).

The orientation signs of the signed frame ([`signs`](@ref)) are solved from the
pairing at construction; [`physical_charges`](@ref) and [`signed_pairing`](@ref) give
the physical objects. **Every chamber is supported** - cyclic chambers included.

With `verify = true` (default) the keystone consistency of the bridge is enforced in
its strict form: there must exist `ε` with `P == −εBε` against
`B = triangulation_quiver(t).B`, and `signed_pairing(cb) == −B` must then hold
entrywise. A failure throws [`ContourError`](@ref) - either the magnitudes disagree
(the numerical contours do not represent the combinatorial classes) or the sign
cocycle does not close around a cycle of the quiver graph.
"""
function charge_basis(prob::SchrodingerProblem, t::IdealTriangulation;
                      margin = nothing, n::Integer = 32, rtol = nothing,
                      verify::Bool = true)
    _require_unpunctured(t)
    m = n_diagonals(t)
    F = _wkb_float(prob)
    contours = Vector{Vector{Complex{F}}}(undef, m)
    charges = Vector{Complex{F}}(undef, m)
    for e in 1:m
        (i, j) = t.diagonal_tp_pair[e]
        contours[e], charges[e] = _oriented_contour(prob, i, j; margin, n, rtol)
    end
    remesh(a, b) = (_oriented_contour(prob, t.diagonal_tp_pair[a]...;
                                      margin, n = n + 15, rtol)[1],
                    _oriented_contour(prob, t.diagonal_tp_pair[b]...;
                                      margin, n = n + 21, rtol)[1])
    _assemble_charge_basis(prob, t, contours, charges, remesh, verify)
end

# Shared tail of both routes: the measured pairing, the signed frame, the keystone.
function _assemble_charge_basis(prob, t, contours::Vector{Vector{Complex{F}}},
                                charges, remesh, verify) where {F}
    m = length(contours)
    P = zeros(Int, m, m)
    for a in 1:m, b in (a + 1):m
        κ = try
            intersection_pairing(prob, contours[a], contours[b])
        catch err
            # a vertex-crossing / near-parallel mesh artefact: retry once on
            # re-meshed (homologous) copies
            err isa ContourError || rethrow()
            ca, cb = remesh(a, b)
            intersection_pairing(prob, ca, cb)
        end
        P[a, b] = κ
        P[b, a] = -κ
    end
    # The signed frame (src/signed_frame.jl): ε with P = −εBε. Solving this IS the
    # keystone check - it subsumes the original `|P| == |B|` and strengthens it, since
    # a solution exists only if the magnitudes agree edgewise AND the sign cocycle
    # closes around every cycle of the quiver graph.
    B = triangulation_quiver(t).B
    ε = _solve_signs(P, B)
    if ε === nothing
        verify && throw(ContourError(
            "keystone mismatch: no orientation ε solves P = −εBε with P = $P and " *
            "B = $B - either the magnitudes disagree (the numerical contours do not " *
            "represent the combinatorial classes) or the sign cocycle fails to " *
            "close around a cycle of the quiver graph"))
        ε = ones(Int, m)
    end
    if verify
        ε .* P .* ε' == -B || throw(ContourError(
            "keystone mismatch: the signed pairing is $(ε .* P .* ε') but the " *
            "triangulation quiver demands −B = $(-B)"))
    end
    ChargeBasis{F}(prob, t, contours, charges, P, ε)
end

"""
    charge_basis(prob::AbstractSchrodingerProblem, g::StokesGraph;
                 margin = nothing, cap = 8, rtol = nothing, verify = true) -> ChargeBasis

The charge-lattice basis read off the Stokes graph `g` directly - the route that works
on **any** surface, punctures included. It builds the triangulation and the strip core
paths together ([`diagonal_core_paths`](@ref)) and encircles each path with a racetrack
([`charge_contour`](@ref)), instead of putting an ellipse around a turning-point pair.

Use this whenever `Q` has poles: there the pair does not name the cycle (on the annulus
both arcs of `Ã(1,1)` join the same two turning points and differ by the loop around
the hole, which is the whole difference between the monopole and the dyon), and the
two-argument method cannot express the difference.

On a disk the two routes agree - same central charges, and the same pairing and signs
exactly - so `charge_basis(prob, t)` remains the polynomial path and is untouched. The
cost here is that the core path runs out to a marked point and back, so the two legs
contribute `±e^{iθ}t` and partly cancel; the cheaper side of each strip is chosen for
exactly that reason, but a fixture whose strips are long relative to `|Z|` will lose
digits.

**Each cycle gets a distinct racetrack half-width.** Two strips routinely share a
boundary ray - on the annulus both do - and equal-width racetracks around it would be
*collinear*, which `intersection_pairing` rightly refuses. Nested widths make the
representatives transversal, which is legitimate because the pairing is a homological
invariant of the lifts and so is unchanged by the perturbation.
"""
function charge_basis(prob::AbstractSchrodingerProblem, g::StokesGraph;
                      margin = nothing, cap::Integer = 8, rtol = nothing,
                      verify::Bool = true)
    _require_simple(g)
    pd, paths = _build_decomposition(g)
    t = _as_triangulation(pd)
    _require_unpunctured(t)
    m = n_diagonals(t)
    F = _wkb_float(prob)
    # strictly decreasing, so no two cycles run collinearly along a shared ray
    width(e, shrink) = if margin === nothing
        _, zs, ze, pts, others = _core_path_parts(prob, paths[e])
        _path_margin(zs, ze, pts, paths[e], others) * shrink / (1 + F(13) / 100 * e)
    else
        F(margin) * shrink / (1 + F(13) / 100 * e)
    end
    contours = Vector{Vector{Complex{F}}}(undef, m)
    charges = Vector{Complex{F}}(undef, m)
    for e in 1:m
        contours[e], charges[e] =
            _oriented_path_contour(prob, paths[e]; margin = width(e, 1), cap, rtol)
    end
    remesh(a, b) = (_oriented_path_contour(prob, paths[a]; margin = width(a, 3 // 4),
                                           cap = cap + 4, rtol)[1],
                    _oriented_path_contour(prob, paths[b]; margin = width(b, 5 // 8),
                                           cap = cap + 7, rtol)[1])
    _assemble_charge_basis(prob, t, contours, charges, remesh, verify)
end

"""
    bridge_seed(cb::ChargeBasis) -> ClusterAlgebras.Seed

The principal-coefficient cluster seed of the charge basis - [`ddp_seed`](@ref) of its
[`signed_pairing`](@ref) (`B = -P`, ledger item 5, applied to the *physical* cycles).

Vertex `j` of the seed carries the physical cycle `ε_j · cb.contours[j]` with
`ε = signs(cb)`, so its y-variable is `y_j = V_γ^{ε_j}` - the Voros symbol of the
stored contour raised to `ε_j`. Feeding `voros_symbol(w, cb.contours[j])` to the DDP
layer is therefore only correct where `ε_j == +1` (always true in the reference
chamber, [`reference_theta`](@ref)); elsewhere the symbol must be inverted first.
"""
bridge_seed(cb::ChargeBasis) = ddp_seed(signed_pairing(cb))
