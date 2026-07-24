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
    problem::SchrodingerProblem
    triangulation::IdealTriangulation
    contours::Vector{Vector{Complex{F}}}
    central_charges::Vector{Complex{F}}
    pairing::Matrix{Int}
    signs::Vector{Int}
end

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

# Decay-oriented contour + central charge for a turning-point pair: the cycle is
# reversed unless Z lies in the closed lower half-plane minus the positive real axis,
# so that Re(Z e^{-iθ_c}) < 0 along the cycle's own wall (ledger item 4) - this is
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
