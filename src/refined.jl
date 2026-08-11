# The refined (motivic) bridge: the BPS chamber of a Schrödinger problem as a quantum
# dilogarithm word, and its wall-crossing invariance.
#
# The classical bridge (src/ddp.jl, src/bps.jl) reads a chamber as a maximal green
# sequence and verifies "crossing the wall of γ_k = the y-mutation μ_k" numerically on
# Borel-summed Voros symbols. This layer is the same chamber read in the quantum torus
# ŷ^α ŷ^β = v^{Λ(α,β)} ŷ^{α+β}: each BPS state contributes a quantum dilogarithm factor,
# and the ordered product is the refined DT invariant of the chamber ([KS08], [DG09]).
# Evaluation is
# ClusterAlgebras' truncated KS algebra (`ks_dilog_product`); nothing about it is
# recomputed here - only the charge lattice, the skew form and the ordering are ours.
#
# ── Convention ledger (each item pinned by an oracle in test/test_refined.jl) ────────
#
# 1. THE SKEW FORM IS THE INTERSECTION PAIRING. Λ = signed_pairing(cb), i.e.
#    Λ(γ_i, γ_j) = ⟨γ_i, γ_j⟩ of the physical cycles, so Λ = −B_seed with B the
#    bridge seed of ledger item 5 (src/ddp.jl). This is exactly ClusterAlgebras'
#    own Λ = −D·B for the same seed (D = I for a quiver from a triangulation), so the
#    physical word and the seed-side word live in one torus by construction. Pinned by:
#    the physical A₂ word reproducing `quantum_dilog_word(bridge_seed(cb), sequence)`
#    factor for factor, and by the pentagon below (the opposite sign of Λ inverts every
#    v-power, and the two sides of the pentagon then disagree at order v).
#
# 2. FACTOR ORDER = SWEEP ORDER. Leftmost factor = first wall crossed as θ decreases
#    from the chamber's θ₀ (ledger item 6), which is `ClusterAlgebras`'s "first mutation
#    leftmost". Pinned by the same reproduction test.
#
# 3. TWO IDENTIFICATIONS, TWO DIFFERENT STATEMENTS - and they disagree, which is the
#    finding of this layer rather than a defect in it.
#      • GEOMETRIC transport (`refined_frame_map`, same problem, same turning points):
#        the cycle of a diagonal continues to the cycle of the diagonal with the same
#        turning-point pair, oriented by the signed frame. Under it two θ-chambers of
#        the cubic differ by ε on the second cycle, and their products are NOT equal.
#        Correct: rotating the half-plane past a BPS ray replaces γ by −γ (a wall of
#        the second kind), and a half-plane product is not invariant under that.
#      • QUIVER isomorphism (`refined_quiver_iso`, the unique permutation matching the
#        skew forms): the only identification available across a *moduli* deformation,
#        where the turning points are relabelled by a continuation we do not track.
#        Under it the cubic z³ − z (2 states) and the cubic z³ − 1.2 (3 states, with
#        the bound state γ₁+γ₂) give the SAME refined product - the quantum pentagon,
#        reached from two potentials. Pinned by `verify_refined_wall_crossing`; the
#        identity itself is [Rei10] / [Kel11].
#    Do not substitute one map for the other: the permutation that matches the skew
#    forms of two θ-chambers is a basis SWAP, and it makes the products agree for a
#    reason that has nothing to do with the physical cycles.
#
# 4. THE REFINED DDP JUMP IS CONJUGATION BY 𝔼(ŷ^{γ_k})^{-1}. `refined_ddp_jump` uses
#    exponent −1 so that its classical limit is literally the DDP factor of ledger
#    item 2, ŷ^{γ_j} ↦ ŷ^{γ_j}(1 + ŷ^{γ_k})^{⟨γ_j,γ_k⟩}; the opposite order gives the
#    inverse exponent. Pinned by comparing the v → 1 limit to `signed_pairing[j,k]`,
#    and (numerically) against `verify_ddp_mutation` on the cubic.
#
# ── what this layer does NOT claim ──────────────────────────────────────────────────
# The refinement is algebraic. A Voros symbol is a Borel-summed *number* and its
# lateral jump is a number, so no `v` can be measured on the WKB side: the q-grading is
# a spin fugacity that the Borel plane does not carry. The v → 1 specialization is the
# only comparison the two pictures admit, and `verify_refined_ddp` is where it is made.
# A genuinely q-deformed WKB side would need non-commuting (quantum) Voros symbols,
# i.e. quantum Fock-Goncharov coordinates - a different package's charter.

using ClusterAlgebras: QuantumDilogWord, ks_dilog_product, ks_dilog_adjoint,
                        ks_classical_limit

"""
    refined_skew(cb::ChargeBasis) -> Matrix{Int}

The skew form of the chamber's quantum torus: `Λ[i,j] = ⟨γ_i, γ_j⟩ =`
[`signed_pairing`](@ref)`(cb)`, so that `ŷ^α ŷ^β = v^{Λ(α,β)} ŷ^{α+β}` with
`v = q^{1/2}`. Equal to `−B` of [`bridge_seed`](@ref) - the same form
`ClusterAlgebras` attaches to that seed (ledger item 1).
"""
refined_skew(cb::ChargeBasis) = signed_pairing(cb)

"""
    refined_dilog_word(sp::BPSSpectrum) -> ClusterAlgebras.QuantumDilogWord

The chamber's quantum dilogarithm word `𝔼(ŷ^{γ₁}) ⋯ 𝔼(ŷ^{γ_ℓ})`: one factor per BPS
state, in sweep order ([`charges`](@ref)), on the lattice of the chamber's
[`ChargeBasis`](@ref) with the skew form [`refined_skew`](@ref).

Every state must carry `Ω(γ) = 1` (true throughout finite type). A state of
nontrivial refined index - a vector multiplet - is not admitted here and never
silently: see [`su2_refined_rays`](@ref) for the affine case, where the index is
supplied explicitly.
"""
function refined_dilog_word(sp::BPSSpectrum)
    all(omega(s) == 1 for s in sp.states) || throw(Resurgence.InvalidArgument(
        "refined_dilog_word needs Ω(γ) = 1 on every state: a state with a nontrivial " *
        "refined index must supply it explicitly (see su2_refined_rays)"))
    QuantumDilogWord(charges(sp), refined_skew(sp.basis), copy(sp.sequence))
end

"""
    refined_dt(sp::BPSSpectrum; truncation_degree = 6) -> Dict

The chamber's **refined DT invariant**: the ordered quantum dilogarithm product of
[`refined_dilog_word`](@ref), evaluated in the quantum torus truncated at total
charge degree `truncation_degree`. Keys are charge vectors in the chamber's basis,
values coefficients in `Frac(ℤ[v])`, `v = q^{1/2}`.

The product itself has no `v → 1` limit (its coefficients have poles there); the
transformation it generates does - see [`refined_ddp_jump`](@ref).
"""
refined_dt(sp::BPSSpectrum; truncation_degree::Int = 6) =
    ks_dilog_product(refined_dilog_word(sp); truncation_degree)

# ── comparing two chambers ──────────────────────────────────────────────────────────

# All permutation matrices of size n, as row-index vectors.
_permutations(n::Int) = n == 0 ? [Int[]] :
    [insert!(copy(p), i, n) for p in _permutations(n - 1) for i in 1:n]

"""
    refined_frame_map(sp1::BPSSpectrum, sp2::BPSSpectrum) -> Matrix{Int}

The **geometric** transport between two θ-chambers of one problem: the basis cycle of
a diagonal continues to the cycle of the diagonal with the same turning-point pair,
and its orientation is the signed frame's, so

    γ_j^{(2)} = ε_j^{(2)} ε_{σ(j)}^{(1)} · γ_{σ(j)}^{(1)},

`σ` the matching of turning-point pairs. A charge `c` of `sp2` becomes `M c` in
`sp1`'s basis. Verified, not assumed: `M Λ₂ Mᵀ == Λ₁` is checked and a mismatch
throws.

Only defined when the two spectra come from the *same* problem: across a deformation
of the potential the turning points are relabelled by continuation, which this
package does not track, so a cross-moduli comparison goes through
[`refined_quiver_iso`](@ref) instead. Throws [`ChamberError`](@ref) otherwise.
"""
function refined_frame_map(sp1::BPSSpectrum, sp2::BPSSpectrum)
    cb1, cb2 = sp1.basis, sp2.basis
    Λ1, Λ2 = refined_skew(cb1), refined_skew(cb2)
    size(Λ1) == size(Λ2) || throw(ChamberError(
        "charge lattices have different ranks: $(size(Λ1, 1)) and $(size(Λ2, 1))"))
    tp1 = [location(t) for t in turning_points(cb1.problem)]
    tp2 = [location(t) for t in turning_points(cb2.problem)]
    (length(tp1) == length(tp2) &&
     all(abs(a - b) ≤ sqrt(eps(Float64)) * (1 + abs(a)) for (a, b) in zip(tp1, tp2))) ||
        throw(ChamberError(
            "geometric transport needs the same turning points on both sides: these " *
            "spectra come from different potentials, where the cycles are related by " *
            "continuation in moduli - use refined_quiver_iso"))
    pairs1, pairs2 = cb1.triangulation.diagonal_tp_pair, cb2.triangulation.diagonal_tp_pair
    ε1, ε2 = signs(cb1), signs(cb2)
    n = size(Λ1, 1)
    M = zeros(Int, n, n)
    for j in 1:n
        i = findfirst(==(pairs2[j]), pairs1)
        i === nothing && throw(ChamberError(
            "diagonal $j of the second chamber has turning-point pair $(pairs2[j]), " *
            "which no diagonal of the first chamber carries: the two triangulations " *
            "are not related by continuation of the diagonals"))
        M[i, j] = ε2[j] * ε1[i]
    end
    M * Λ2 * transpose(M) == Λ1 || throw(ChamberError(
        "the geometric transport $M does not carry the pairing $Λ2 onto $Λ1: the " *
        "cycles named by the same turning-point pair are not the same homology class"))
    M
end

"""
    refined_quiver_iso(sp1::BPSSpectrum, sp2::BPSSpectrum) -> Matrix{Int}

The identification of two chambers' charge lattices as **BPS quivers**: the unique
permutation matrix `M` with `M Λ₂ Mᵀ = Λ₁`, carrying `sp2`'s basis onto `sp1`'s (a
charge `c` becomes `M c`).

Only permutations are searched, deliberately: the charges of a chamber are positive
in its own basis (they are `c`-vectors of a green sequence), so an identification
that keeps both spectra in one positive cone is a relabelling of the basis. This is
the identification available across a *moduli* deformation, where the geometric
transport of [`refined_frame_map`](@ref) is not. Throws [`ChamberError`](@ref) when
no such `M` exists, or when more than one does (a quiver automorphism - then the
comparison would not be a statement).
"""
function refined_quiver_iso(sp1::BPSSpectrum, sp2::BPSSpectrum)
    Λ1, Λ2 = refined_skew(sp1.basis), refined_skew(sp2.basis)
    size(Λ1) == size(Λ2) || throw(ChamberError(
        "charge lattices have different ranks: $(size(Λ1, 1)) and $(size(Λ2, 1))"))
    n = size(Λ1, 1)
    found = Matrix{Int}[]
    for p in _permutations(n)
        M = zeros(Int, n, n)
        for (i, j) in enumerate(p)
            M[i, j] = 1
        end
        M * Λ2 * transpose(M) == Λ1 && push!(found, M)
    end
    isempty(found) && throw(ChamberError(
        "no relabelling of the charge basis identifies the two chambers: the skew " *
        "forms $Λ2 and $Λ1 are not permutation-equivalent, so the two BPS quivers " *
        "are not isomorphic and their refined products are not comparable"))
    length(found) == 1 || throw(ChamberError(
        "$(length(found)) relabellings identify the two chambers (the BPS quiver has " *
        "a nontrivial automorphism), so comparing the products is ambiguous"))
    only(found)
end

"""
    verify_refined_wall_crossing(sp1::BPSSpectrum, sp2::BPSSpectrum;
                                 truncation_degree = 6) -> NamedTuple

Verify **refined (motivic) wall-crossing** between two BPS chambers: their quantum
dilogarithm products agree once `sp2`'s charges are relabelled into `sp1`'s basis by
[`refined_quiver_iso`](@ref). Returns `(; equal, map, n_states, product)`.

The physical instance is a wall of marginal stability in *moduli*: on one side a
2-state chamber, on the other a 3-state chamber whose extra state is the bound state
`γ₁+γ₂`, and the identity is the quantum pentagon - the refined DT invariant does not
see the decay. The identity itself is Reineke/Keller's; what is verified here is that
two *potentials* land on it.

It is **not** a statement about two `θ`-chambers of one problem. Under the geometric
transport of [`refined_frame_map`](@ref) those products differ, and must: rotating
the half-plane past a BPS ray replaces `γ` by `−γ` (a wall of the second kind, ledger
item 3).
"""
function verify_refined_wall_crossing(sp1::BPSSpectrum, sp2::BPSSpectrum;
                                      truncation_degree::Int = 6)
    M = refined_quiver_iso(sp1, sp2)
    w2 = QuantumDilogWord([M * c for c in charges(sp2)],
                          refined_skew(sp1.basis), copy(sp2.sequence))
    p1 = refined_dt(sp1; truncation_degree)
    p2 = ks_dilog_product(w2; truncation_degree)
    (; equal = p1 == p2, map = M, n_states = (n_states(sp1), n_states(sp2)),
       product = p1)
end

# ── the refined DDP jump ────────────────────────────────────────────────────────────

"""
    refined_ddp_jump(cb::ChargeBasis, k::Integer, j::Integer;
                     truncation_degree = 6) -> Dict

The image of `ŷ^{γ_j}` under the quantum wall-crossing transformation at the wall of
`γ_k`: conjugation by `𝔼(ŷ^{γ_k})^{-1}` in the chamber's quantum torus, truncated at
`truncation_degree` extra charge degrees.

This is the `q`-deformation of the DDP jump. Its classical limit
([`ClusterAlgebras.ks_classical_limit`](@ref)) is exactly
`ŷ^{γ_j} (1 + ŷ^{γ_k})^{⟨γ_j,γ_k⟩}` - the factor `ddp_transform` applies to a
Borel-summed Voros value, and the `y`-mutation `μ_k` of the bridge seed (ledger
item 4).
"""
function refined_ddp_jump(cb::ChargeBasis, k::Integer, j::Integer;
                          truncation_degree::Int = 6)
    n = n_charges(cb)
    (1 ≤ k ≤ n && 1 ≤ j ≤ n) ||
        throw(Resurgence.InvalidArgument("vertices k = $k, j = $j out of range 1:$n"))
    γk = [i == k ? 1 : 0 for i in 1:n]
    γj = [i == j ? 1 : 0 for i in 1:n]
    w = QuantumDilogWord([γk], refined_skew(cb), [Int(k)], [1], [0], [-1])
    ks_dilog_adjoint(w, γj; truncation_degree)
end

"""
    verify_refined_ddp(cb::ChargeBasis, k::Integer; truncation_degree = 6)
        -> NamedTuple

Verify that the refined jump at the wall of `γ_k` collapses onto the classical DDP
formula: for every `j`, the `v → 1` limit of [`refined_ddp_jump`](@ref) equals
`ŷ^{γ_j}(1 + ŷ^{γ_k})^{p}` with `p = ⟨γ_j, γ_k⟩ =` [`signed_pairing`](@ref)`[j,k]`,
the exact-integer Stokes constant of [`verify_ddp`](@ref).

Returns `(; matches, exponents, classical)`: `exponents[j]` is the pairing the
classical limit exhibits (read off the binomial coefficients), and `matches` is
`true` when all of them agree with the pairing matrix.
"""
function verify_refined_ddp(cb::ChargeBasis, k::Integer; truncation_degree::Int = 6)
    n = n_charges(cb)
    P = signed_pairing(cb)
    exps = zeros(Int, n)
    classical = Vector{Any}(undef, n)
    ok = true
    for j in 1:n
        cl = ks_classical_limit(refined_ddp_jump(cb, k, j; truncation_degree))
        classical[j] = cl
        expected = _binomial_series(P[j, k], k, j, n, truncation_degree)
        exps[j] = P[j, k]
        cl == expected || (ok = false)
    end
    (; matches = ok, exponents = exps, classical)
end

# ŷ^{γ_j}(1 + ŷ^{γ_k})^p as a truncated classical (commutative) series - the DDP factor.
function _binomial_series(p::Int, k::Integer, j::Integer, n::Int, D::Int)
    out = Dict{Vector{Int}, Rational{BigInt}}()
    coeff = one(Rational{BigInt})
    for m in 0:D
        α = zeros(Int, n)
        α[j] += 1
        α[k] += m
        iszero(coeff) || (out[α] = coeff)
        coeff *= Rational{BigInt}(p - m, m + 1)
    end
    out
end
