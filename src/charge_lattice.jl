# The numerical half of the M4 cluster bridge: one closed contour per diagonal of the
# ideal triangulation — the charge lattice basis γ₁ … γ_{d−1}. Each γ_j encircles the
# turning-point pair of its diagonal, is oriented by the DDP decay rule (ledger item 4
# of src/ddp.jl: `Re(Z_j e^{-iθ_c}) < 0` along its own wall, i.e. `Z_j` in the closed
# lower half-plane minus the positive real axis), and carries the central charge
# `Z_j = ∮√Q` and the intersection pairing `P[i,j] = ⟨γ_i, γ_j⟩`.
#
# The keystone consistency of the bridge — checked at construction (`verify = true`)
# — is that the numerically measured pairing reproduces the combinatorial exchange
# matrix of the triangulation: `|P| == |B|` entrywise, so `B_seed = −P` (what
# `ddp_seed` builds, ledger item 5) is the triangulation quiver up to reorientation
# of basis cycles.
#
# ── M4 finding (2026-07-16, discovered by this oracle) ─────────────────────────────
# The all-decay orientation and the seed structure are compatible exactly on *tree
# chambers* — chambers whose diagonal-adjacency graph is a tree (no triangle of the
# triangulation has all three edges internal). There `−P(decay)` is finite type and
# is the physical seed (cubic-verified, ledger item 5), though it may differ from the
# fixed-counterclockwise-rule `triangulation_quiver` by a harmless sign conjugation.
# On *cyclic chambers* (a fully internal triangle — generic for degree ≥ 5) the
# decay-oriented pairing orients that quiver triangle acyclically, so `−P(decay)` is
# cluster-INFINITE and cannot be a seed: the correct frame needs the Iwaki–Nakanishi
# signed-flip conventions, deferred beyond M4 (decision 2026-07-16). `charge_basis`
# therefore throws `ChamberError` on cyclic chambers; every polynomial has tree
# chambers (e.g. the quintic's top chamber).

"""
    ChargeBasis{F}

The charge-lattice basis attached to an [`IdealTriangulation`](@ref): per diagonal a
closed `contour` (decay-oriented, ledger item 4), its central charge
`Z_j = ∮√Q ∈ central_charges`, and the integer intersection-pairing matrix
`pairing[i,j] = ⟨γ_i, γ_j⟩`. Build with [`charge_basis`](@ref); feed to
[`bridge_seed`](@ref) / [`voros_symbol`](@ref).
"""
struct ChargeBasis{F}
    problem::SchrodingerProblem
    triangulation::IdealTriangulation
    contours::Vector{Vector{Complex{F}}}
    central_charges::Vector{Complex{F}}
    pairing::Matrix{Int}
end

n_charges(cb::ChargeBasis) = length(cb.contours)
central_charges(cb::ChargeBasis) = cb.central_charges

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
        "point stays inside the ellipse — pass an explicit margin"))
end

# Decay-oriented contour + central charge for a turning-point pair: the cycle is
# reversed unless Z lies in the closed lower half-plane minus the positive real axis,
# so that Re(Z e^{-iθ_c}) < 0 along the cycle's own wall (ledger item 4) — this is
# what makes every basis symbol a legal wall symbol for verify_ddp/verify_ddp_mutation.
function _oriented_contour(prob, i, j; margin, n, rtol)
    c = charge_contour(prob, i, j; margin, n)
    Z = period_integral(prob, c; rtol)
    # numerically real Z (wall at θ_c = 0) must resolve by the real part, not by the
    # sign of a roundoff-sized imaginary residue
    tol = sqrt(eps(typeof(real(Z)))) * (1 + abs(Z))
    if imag(Z) > tol || (abs(imag(Z)) ≤ tol && real(Z) > 0)
        c = reverse(c)
        Z = -Z
    end
    c, Z
end

"""
    charge_basis(prob::SchrodingerProblem, t::IdealTriangulation;
                 margin = nothing, n = 32, rtol = nothing, verify = true) -> ChargeBasis

The charge-lattice basis of the triangulation `t`: for each diagonal (in `t`'s
canonical order) a closed contour around its turning-point pair
([`charge_contour`](@ref)), oriented to decay along its own wall (ledger item 4 —
`Re(Z_j e^{-iθ_c}) < 0`), with central charges `Z_j = ∮√Q` and the intersection
pairing `P[i,j] = ⟨γ_i, γ_j⟩` from [`intersection_pairing`](@ref).

With `verify = true` (default) the keystone consistency of the bridge is enforced:
the chamber must be a *tree chamber* (no triangle of `t` has all three edges
internal — otherwise the decay frame is not a seed frame and [`ChamberError`](@ref)
is thrown; pick a different `θ`, e.g. in the top chamber), and the numerical pairing
must reproduce the combinatorial exchange matrix up to cycle reorientation,
`abs.(P) == abs.(triangulation_quiver(t).B)` — a mismatch throws
[`ContourError`](@ref).
"""
function charge_basis(prob::SchrodingerProblem, t::IdealTriangulation;
                      margin = nothing, n::Integer = 32, rtol = nothing,
                      verify::Bool = true)
    if verify && any(tri -> all(e -> t.is_diagonal[e], tri), t.triangles)
        throw(ChamberError(
            "cyclic chamber: a triangle of the triangulation has all three edges " *
            "internal, so the decay-oriented cycle basis is not a seed basis (the " *
            "Iwaki–Nakanishi signed-flip frame is deferred beyond M4) — trace the " *
            "Stokes graph at a θ in a tree chamber (the top chamber always is one)"))
    end
    m = n_diagonals(t)
    F = _wkb_float(prob)
    contours = Vector{Vector{Complex{F}}}(undef, m)
    charges = Vector{Complex{F}}(undef, m)
    for e in 1:m
        (i, j) = t.diagonal_tp_pair[e]
        contours[e], charges[e] = _oriented_contour(prob, i, j; margin, n, rtol)
    end
    P = zeros(Int, m, m)
    for a in 1:m, b in (a + 1):m
        κ = try
            intersection_pairing(prob, contours[a], contours[b])
        catch err
            # a vertex-crossing / near-parallel mesh artefact: retry once on
            # re-meshed (homologous) copies with coprime vertex counts
            err isa ContourError || rethrow()
            ca, _ = _oriented_contour(prob, t.diagonal_tp_pair[a]...;
                                      margin, n = n + 15, rtol)
            cb, _ = _oriented_contour(prob, t.diagonal_tp_pair[b]...;
                                      margin, n = n + 21, rtol)
            intersection_pairing(prob, ca, cb)
        end
        P[a, b] = κ
        P[b, a] = -κ
    end
    if verify
        B = triangulation_quiver(t).B
        abs.(P) == abs.(B) || throw(ContourError(
            "keystone mismatch: the intersection pairing of the decay-oriented " *
            "cycles is P = $P, but the triangulation quiver demands |P| = |B| with " *
            "B = $B — the numerical contours do not represent the combinatorial " *
            "classes"))
    end
    ChargeBasis{F}(prob, t, contours, charges, P)
end

"""
    bridge_seed(cb::ChargeBasis) -> ClusterAlgebras.Seed

The principal-coefficient cluster seed of the charge basis — [`ddp_seed`](@ref) of its
pairing matrix (`B = -P`, ledger item 5). Vertex `j` of the seed carries the cycle
`cb.contours[j]`; its Voros symbol models the y-variable `y_j`.
"""
bridge_seed(cb::ChargeBasis) = ddp_seed(cb.pairing)
