# DDP layer: the Delabaere–Dillinger–Pham connection formula as the seam between the
# summed Voros symbols (Resurgence side) and cluster y-mutation (ClusterAlgebras side).
# This is the convention layer the cluster bridge builds on. The symbols never
# move: a VorosSymbol is computed once from a fixed contour, and crossing the critical
# phase θ_c only rotates the Laplace ray in the Borel plane, so wall-crossing is a
# lateral-sum statement.
#
# In physics language the DDP jump factor (1 + V_{γ'})^{⟨γ,γ'⟩} is the
# Kontsevich–Soibelman symplectomorphism attached to the BPS ray of γ', and the
# integer Stokes constant is the Donaldson–Thomas statement Ω(γ')·⟨γ,γ'⟩.
#
# ── Convention ledger (each item pinned by an oracle in test/test_ddp.jl) ──────────
#
# 1. Lateral sums / Stokes automorphism are Resurgence's, with its Euler-pinned sign
#    disc_θ = s_θ⁺ − s_θ⁻ and s₋ = s₊ ∘ 𝔖.
# 2. DDP formula. For a wall at phase θ_c with (decaying) wall symbol V_{γ'} and any
#    cycle γ:   s_θ⁺(log V_γ) − s_θ⁻(log V_γ) = ⟨γ, γ'⟩ · log(1 + V_{γ'}),
#    with V_{γ'} evaluated by the median sum. As a factor this reads
#    s₊(V_γ) = s₋(V_γ) · (1 + V_{γ'})^{⟨γ,γ'⟩} - ddp_transform maps the minus-side
#    lateral value to the plus-side one. The Stokes constant IS the integer
#    intersection number Ω(γ')·⟨γ,γ'⟩ with Ω = 1 for a saddle between two simple
#    turning points. Pinned by the cubic assertion round(jump/log1p(V')) == ⟨γ,γ'⟩
#    and by the direction oracle on the exponentiated values.
# 3. ⟨γ,γ'⟩ = intersection_pairing(prob, cγ, cγ') - the signed same-sheet crossing
#    count of the canonical lifts, a crossing of directions (d, d') counting
#    sign(Im(conj(d)·d')). The canonical lift takes the principal √Q branch at the
#    contour's first vertex - the same rule period_integral uses, so the pairing and
#    the classical periods can never disagree about orientation.
# 4. The wall symbol must DECAY along the wall: Re(Z_{γ'} e^{-iθ_c}) < 0. Enforced by
#    a typed error and never flipped silently - reversing the contour flips both
#    Z_{γ'} and every ⟨γ,γ'⟩, which shifts the DDP formula by the tropical monomial;
#    an automatic flip is exactly the bug class this layer exists to kill.
# 5. Cluster dictionary: y_j = V_{γ_j} with exchange matrix B = −P, P[i,j] = ⟨γ_i,γ_j⟩
#    (ddp_seed) - i.e. b_kj = ⟨γ_j, γ_k⟩: an exchange-matrix entry is the pairing of
#    the affected charge with the mutating one. Crossing the wall of γ_k (the initial
#    seed lives on the plus side) is the y-mutation μ_k in the Fomin–Zelevinsky
#    convention of ClusterAlgebras (y'_k = y_k⁻¹,
#    y'_j = y_j y_k^{[b_kj]₊} (1+y_k)^{−b_kj}): the mutated y-variables evaluated at
#    the plus-side values equal the minus-side values of the μ_k-mutated charges
#    V_{γ'_j} = V_{γ_j} V_{γ_k}^{[b_kj]₊}. Pinned by verify_ddp_mutation on the cubic
#    (the opposite B-sign, i.e. the opposite mutation direction, fails by ≥ 10×;
#    B = −P was measured, B = +P off by ~0.14 where −P matched at ~2e-5).
# 6. Chamber ↔ MGS: with these conventions the charges of a maximal green sequence
#    (ClusterAlgebras.ordered_c_vectors) list the chamber's BPS states in the order
#    their walls are crossed as θ DECREASES from the chamber's own θ₀ - i.e. by
#    increasing `mod(θ₀ − θ_c, π)`. Pinned by the cubic's two-state chamber
#    (the reverse order matches no MGS of the seed).
#
#    The θ₀ → π⁻ specialization of this is the plain "θ-decreasing order of the
#    saddle phases" in which the item was originally pinned, and it is what the top
#    (reference) chamber realizes. The original bridge layer used that absolute form everywhere, which is
#    why it enumerated correctly in the top chamber and misordered elsewhere; the
#    chamber-relative form (src/bps.jl `_sweep_chamber`) is the general statement.

using Resurgence: borel, pade, lateral_sum

# Borel–Padé approximant of the quantum tail. reduce = true is mandatory: the odd-only
# Voros series has alternating zero coefficients, so the equal-degree Padé system is
# always degenerate (documented at quantum_series).
function _voros_approximant(vs::VorosSymbol; order = nothing)
    B = borel(quantum_series(vs))
    B, pade(B; order, reduce = true)
end

function _tail_sum(B, r, ħ; theta, side, tilt, rtol)
    side === :median &&
        return (lateral_sum(B, r, ħ; θ = theta, side = :plus, tilt, rtol) +
                lateral_sum(B, r, ħ; θ = theta, side = :minus, tilt, rtol)) / 2
    lateral_sum(B, r, ħ; θ = theta, side, tilt, rtol)
end

_value_type(ħ) = Complex{typeof(float(real(ħ)))}

# log V = v₋₁/ħ + s(tail): the exponent of the Voros symbol. All jump arithmetic
# happens here in log space (no overflow, and jumps compose additively).
function _log_voros(vs::VorosSymbol, ħ::Number; theta::Real = 0, side::Symbol = :median,
                    order = nothing, tilt::Real = 1 // 100, rtol = nothing)
    side in (:plus, :minus, :median) || throw(Resurgence.InvalidArgument(
        "side must be :plus, :minus or :median, got :$side"))
    B, r = _voros_approximant(vs; order)
    _value_type(ħ)(classical_period(vs) / ħ + _tail_sum(B, r, ħ; theta, side, tilt, rtol))
end

"""
    voros_value(vs::VorosSymbol, ħ; theta = 0, side = :median, order = nothing,
                tilt = 1//100, rtol = nothing) -> Complex

The Borel-summed value ``V = \\exp(v_{-1}/ħ + s_θ(\\text{tail}))`` of the Voros
symbol at `ħ`: the quantum tail is summed by Resurgence's `borel`/`pade`
(`reduce = true`, mandatory for the odd-only series) and the lateral Laplace sum
along the ray `theta` on the given `side` (`:plus`, `:minus`, or `:median` - the
average of the two laterals, the canonical value *on* a Stokes ray). `order`
truncates the Padé; `tilt`/`rtol` pass through to `Resurgence.lateral_sum`.
Precision follows the type of `ħ`.
"""
function voros_value(vs::VorosSymbol, ħ::Number; theta::Real = 0, side::Symbol = :median,
                     order = nothing, tilt::Real = 1 // 100, rtol = nothing)
    exp(_log_voros(vs, ħ; theta, side, order, tilt, rtol))
end

"""
    ddp_transform(V, Vp, pairing::Integer) -> Number

The DDP transform on a summed Voros value across the wall of the (decaying) wall
symbol with value `Vp`: ``V ↦ V\\,(1 + V_{γ'})^{⟨γ,γ'⟩}``. Applied to the
**minus-side** lateral value it yields the **plus-side** one (item 2 of the
convention ledger, direction pinned by the cubic oracle). The integer `pairing` is
`⟨γ,γ'⟩ =` [`intersection_pairing`](@ref)`(prob, cγ, cγ')`.
"""
ddp_transform(V::Number, Vp::Number, pairing::Integer) = V * (1 + Vp)^pairing

function _require_decaying(vs::VorosSymbol, theta)
    Z = classical_period(vs)
    re = real(Z * cis(-theta))
    # strictly decaying, with a roundoff margin: a marginal Re ≈ 0 means the charge's
    # wall is elsewhere - passing it as the wall symbol of this ray is a usage error
    re < -sqrt(eps(typeof(re))) * (1 + abs(Z)) || throw(Resurgence.InvalidArgument(
        "the wall Voros symbol must decay along θ = $theta, but Re(Z e^{-iθ}) = $re ≥ 0; " *
        "reverse the contour orientation (this flips Z and the pairing - it is never " *
        "flipped automatically, since that shifts the DDP formula by the tropical " *
        "monomial)"))
end

"""
    verify_ddp(vγ::VorosSymbol, vγp::VorosSymbol, pairing::Integer, ħ; theta = 0,
               order = nothing, tilt = 1//100, rtol = nothing) -> NamedTuple

Numerically verify the DDP connection formula at the wall `theta`: the measured jump
``s_θ^+(\\log V_γ) − s_θ^-(\\log V_γ)`` of the lateral sums of `vγ`'s tail against the
prediction ``⟨γ,γ'⟩·\\log(1 + V_{γ'})`` with `vγp` the (decaying - enforced) wall
symbol evaluated by its median sum. Returns

  - `jump_measured`, `jump_predicted` - the two sides of the identity (log space);
  - `kappa_measured = jump_measured / log1p(V')` - the Stokes-constant readout: an
    **exact integer** (`= pairing`, the DT-invariant statement) up to Padé error;
  - `residual`, `relative_residual` (absolute when `pairing == 0`), `v_prime`.
"""
function verify_ddp(vγ::VorosSymbol, vγp::VorosSymbol, pairing::Integer, ħ::Number;
                    theta::Real = 0, order = nothing, tilt::Real = 1 // 100,
                    rtol = nothing)
    _require_decaying(vγp, theta)
    C = _value_type(ħ)
    B, r = _voros_approximant(vγ; order)
    jump_measured = C(_tail_sum(B, r, ħ; theta, side = :plus, tilt, rtol) -
                      _tail_sum(B, r, ħ; theta, side = :minus, tilt, rtol))
    v_prime = voros_value(vγp, ħ; theta, side = :median, order, tilt, rtol)
    logfac = log1p(v_prime)
    jump_predicted = pairing * logfac
    residual = abs(jump_measured - jump_predicted)
    relative_residual = pairing == 0 ? residual : residual / abs(jump_predicted)
    (; jump_measured, jump_predicted, kappa_measured = jump_measured / logfac,
       residual, relative_residual, v_prime)
end

# ── the pairing ────────────────────────────────────────────────────────────────────

"""
    intersection_pairing(prob::AbstractSchrodingerProblem, c1, c2) -> Int

The homological intersection number ``⟨γ_1, γ_2⟩`` of the canonical lifts of the two
closed contours to the spectral double cover ``w² = Q(z)``: the signed count of
same-sheet crossings, a crossing of segment directions `(d1, d2)` counting
`sign(Im(conj(d1)·d2))`. The canonical lift takes the principal `√Q` branch at each
contour's first vertex - identical to the branch rule of [`period_integral`](@ref),
so the pairing is orientation-coherent with the classical periods (reversing a
contour flips both). Antisymmetric; `0` for disjoint or homologous cycles.

Throws [`ContourError`](@ref) for a non-transversal crossing (touching or nearly
parallel segments, a crossing at a polygon vertex - re-mesh with a different `n` or
`margin`) or a crossing too close to a turning point.
"""
function intersection_pairing(prob::AbstractSchrodingerProblem, c1::AbstractVector,
                              c2::AbstractVector)
    F = promote_type(_contour_float(_as_complex.(c1)), _contour_float(_as_complex.(c2)))
    verts1 = Complex{F}[_as_complex(z) for z in c1]
    verts2 = Complex{F}[_as_complex(z) for z in c2]
    (length(verts1) ≥ 3 && length(verts2) ≥ 3) ||
        throw(ContourError("intersection_pairing needs closed contours with ≥ 3 vertices"))
    push!(verts1, verts1[1])
    push!(verts2, verts2[1])
    K1, K2 = length(verts1) - 1, length(verts2) - 1
    atol = max(_period_atol(prob, verts1, F), _period_atol(prob, verts2, F))
    S1, PA1, UW1, _ = _branch_table(prob, verts1, K1, atol)
    S2, PA2, UW2, _ = _branch_table(prob, verts2, K2, atol)
    tol = sqrt(eps(F))
    total = 0
    for i in 1:K1, j in 1:K2
        a1, d1 = verts1[i], verts1[i + 1] - verts1[i]
        a2, d2 = verts2[j], verts2[j + 1] - verts2[j]
        det = imag(conj(d1) * d2)
        rhs = a2 - a1
        if abs(det) ≤ tol * abs(d1) * abs(d2)
            # (near-)parallel: transversality fails only if the lines also (nearly)
            # coincide and the segments overlap - reject that; otherwise no crossing.
            dist = abs(imag(conj(d1) * rhs)) / abs(d1)
            if dist ≤ tol * (abs(d1) + abs(d2)) &&
               min(abs(rhs), abs(a2 + d2 - a1), abs(a2 - a1 - d1)) ≤ abs(d1) + abs(d2)
                throw(ContourError("contours have (nearly) collinear overlapping " *
                                   "segments near z ≈ $a1; crossings must be " *
                                   "transversal - re-mesh (change n or margin)"))
            end
            continue
        end
        t = imag(conj(rhs) * d2) / det
        u = imag(conj(rhs) * d1) / det
        (-tol < t < 1 + tol && -tol < u < 1 + tol) || continue
        if !(tol < t < 1 - tol && tol < u < 1 - tol)
            throw(ContourError("contours cross at a polygon vertex near " *
                               "z ≈ $(a1 + t * d1); re-mesh (change n or margin)"))
        end
        z = a1 + t * d1
        abs(prob(z)) ≤ atol &&
            throw(ContourError("contours cross at (or too close to) a turning point " *
                               "at z ≈ $z"))
        u1, _, _, _ = _sqrt_Q(prob, (F(i - 1) + t) / K1, S1, PA1, UW1, verts1, K1)
        u2, _, _, _ = _sqrt_Q(prob, (F(j - 1) + u) / K2, S2, PA2, UW2, verts2, K2)
        if abs(u1 - u2) < abs(u1 + u2)          # same sheet of the double cover
            total += Int(sign(det))
        end
    end
    total
end

# ── the cluster seam ────────────────────────────────────────────────────────────────

"""
    ddp_seed(P::AbstractMatrix{<:Integer}) -> ClusterAlgebras.Seed

The principal-coefficient cluster seed of a family of cycles with intersection-pairing
matrix `P[i,j] = ⟨γ_i, γ_j⟩` (antisymmetric - validated): exchange matrix `B = −P`,
i.e. `b_kj = ⟨γ_j, γ_k⟩` (item 5 of the convention ledger, pinned by the
wall-crossing oracle), extended to principal coefficients so its `y_variables` model
the Voros symbols. This is where the ClusterAlgebras dependency of the bridge comes
alive.
"""
function ddp_seed(P::AbstractMatrix{<:Integer})
    size(P, 1) == size(P, 2) || throw(Resurgence.InvalidArgument(
        "pairing matrix must be square, got $(size(P))"))
    B = -Matrix{Int}(P)
    B == -transpose(B) || throw(Resurgence.InvalidArgument(
        "pairing matrix must be antisymmetric (P[i,j] = ⟨γ_i,γ_j⟩ = -P[j,i])"))
    ClusterAlgebras.extend(ClusterAlgebras.Seed(ClusterAlgebras.Quiver(B)))
end

# Evaluate an AbstractAlgebra Frac(ZZ[y…]) element at complex values by manual Horner
# over its terms - AbstractAlgebra never leaks into (or out of) the public API.
function _eval_mpoly(p, vals::AbstractVector{C}) where {C<:Complex}
    acc = zero(C)
    for (c, e) in zip(AA.coefficients(p), AA.exponent_vectors(p))
        term = C(c)
        for (v, m) in zip(vals, e)
            m == 0 && continue
            term *= v^m
        end
        acc += term
    end
    acc
end

_eval_frac(y, vals) = _eval_mpoly(AA.numerator(y), vals) /
                      _eval_mpoly(AA.denominator(y), vals)

"""
    verify_ddp_mutation(vsyms, seed, k, ħ; theta = 0, order = nothing, tilt = 1//100,
                        rtol = nothing) -> NamedTuple

Numerically verify "crossing the wall of ``γ_k`` = the y-mutation ``μ_k``" - the
Iwaki–Nakanishi dictionary statement the cluster bridge consumes. `vsyms[i]` is the Voros
symbol of the cycle attached to vertex `i` of the initial principal-coefficient
`seed` (from [`ddp_seed`](@ref)); `vsyms[k]` is the wall symbol at `theta`
(decaying - enforced). Compares, at `ħ`:

  - `predicted[j]` - the mutated y-variable ``y'_j`` of `mutate(seed, k)` evaluated
    at the plus-side values ``y_i = s_θ^+(V_{γ_i})``;
  - `measured[j]` - the minus-side value of the **mutated-basis** symbol
    ``V_{γ'_j} = V_{γ_j} V_{γ_k}^{[b_{kj}]_+}`` (and ``V_{γ'_k} = V_{γ_k}^{-1}`` by
    its median sum - a cycle is continuous across its own wall, ``⟨γ_k,γ_k⟩ = 0``).

Returns `(predicted, measured, residuals, max_residual)` with
`residuals = abs.(predicted ./ measured .- 1)`.
"""
function verify_ddp_mutation(vsyms::AbstractVector{<:VorosSymbol},
                             seed::ClusterAlgebras.Seed, k::Integer, ħ::Number;
                             theta::Real = 0, order = nothing, tilt::Real = 1 // 100,
                             rtol = nothing)
    n = seed.quiver.n_mutable
    length(vsyms) == n || throw(Resurgence.InvalidArgument(
        "need one Voros symbol per mutable vertex: got $(length(vsyms)) for n = $n"))
    1 ≤ k ≤ n || throw(Resurgence.InvalidArgument("vertex k = $k out of range 1:$n"))
    all(ClusterAlgebras.cmatrix(seed)[i, j] == (i == j ? 1 : 0) for i in 1:n, j in 1:n) ||
        throw(Resurgence.InvalidArgument(
            "seed must be the initial principal seed (C = I): the y-variables of the " *
            "mutated seed are rational functions in the *initial* y's"))
    _require_decaying(vsyms[k], theta)
    C = _value_type(ħ)
    yplus = C[voros_value(v, ħ; theta, side = :plus, order, tilt, rtol) for v in vsyms]
    yminus = C[voros_value(v, ħ; theta, side = :minus, order, tilt, rtol) for v in vsyms]
    ymed_k = voros_value(vsyms[k], ħ; theta, side = :median, order, tilt, rtol)
    ys = ClusterAlgebras.y_variables(ClusterAlgebras.mutate(seed, k))
    predicted = C[_eval_frac(y, yplus) for y in ys]
    Bmat = seed.quiver.B
    measured = C[j == k ? inv(ymed_k) : yminus[j] * ymed_k^max(Bmat[k, j], 0)
                 for j in 1:n]
    residuals = abs.(predicted ./ measured .- 1)
    (; predicted, measured, residuals, max_residual = maximum(residuals))
end
