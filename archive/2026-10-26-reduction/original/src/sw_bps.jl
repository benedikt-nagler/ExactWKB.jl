# SU(2) BPS spectrum + Kontsevich–Soibelman wall-crossing - rung 2 of the Seiberg–Witten
# SU(2) ascent. Unlike `sw_curve.jl` this file is a **consumer**, not firewalled: it bridges
# the rung-1 SW periods (`central_charge(sw, u, (n_m,n_e))`) to the affine-Kronecker cluster
# machinery of `ClusterAlgebras`. It computes the BPS particle content of pure SU(2) in both
# chambers and verifies the wall-crossing that relates them.
#
# ── Physics ledger (charges, pinned) ──────────────────────────────────────────────────
# Physical charge basis (n_m, n_e), symplectic pairing ⟨(m,e),(m',e')⟩ = m e' − e m'. The BPS
# quiver is the Kronecker quiver (B = [0 2; -2 0], ⟨γ₁,γ₂⟩ = 2 - affine Ã₁, cluster-infinite),
# abstract basis γ₁, γ₂:
#   γ₁ = monopole = (1,0),  γ₂ = dyon = (−1,2),  δ = γ₁+γ₂ = (0,2) = W-boson (Ω = −2).
#   preprojective dyons  P_n = (n+1)γ₁ + n γ₂ = (1, 2n),      n ≥ 0  (Ω = 1);
#   preinjective  dyons  I_n = n γ₁ + (n+1)γ₂ = (−1, 2n+2),   n ≥ 0  (Ω = 1).
#   Strong chamber = {γ₁, γ₂} (2 states). Weak chamber = {P_n} ∪ {δ, Ω=−2} ∪ {I_n}.
# Central charges Z_γ(u) = n_m a_D + n_e a come from rung 1; both dyon towers' phases
# accumulate onto arg Z_δ = arg(a) from opposite sides (marginal stability).
#
# ── Deformed-wall ledger (ours, not forced) ───────────────────────────────────────────
# 1. Under the Nekrasov–Shatashvili deformation the wall of marginal stability is taken to
#    be the zero set of Im(a_D/a) with **both** periods deformed to the same (ħ, order).
#    The alternative - deform only one period, or mix orders - is not excluded by anything
#    in the physics; this choice is the one that keeps the wall the locus where the BPS
#    rays of the two chamber spectra actually align, since both spectra are built from the
#    same deformed pair. Pinned by: ms_wall(...; ħ = 0) reproducing the classical wall
#    exactly, and by the deformed wall converging to it as ħ → 0.
# 2. The real-axis shortcut in `sw_chamber` (real u is :strong iff |u| < 2Λ²) is retained
#    under deformation. Justification: the NS corrections are real-coefficient
#    Picard–Fuchs operators applied to the classical periods, so they inherit the Schwarz
#    reflection symmetry that puts the wall's only real crossings at ±2Λ².
# 3. The charge content and Ω values are NOT deformed - only the central charges are. The
#    NS deformation is a deformation of the periods, not of the BPS quiver.
#
# The wall-crossing is checked classically and self-validatingly (see `verify_su2_wall_crossing`):
# closure of the ray-automorphism identity to a given degree simultaneously pins the spectrum
# content, the Ω values (including the −2 vector factor), the ordering, and the (here trivial -
# the pairing is even) sign convention. No externally supplied "answer" is used.

# abstract Kronecker basis → physical (n_m, n_e): a·γ₁ + b·γ₂ with γ₁=(1,0), γ₂=(−1,2)
_phys_charge(a::Integer, b::Integer) = (a - b, 2b)

"""
    su2_bps_quiver() -> ClusterAlgebras.Quiver

The BPS quiver of pure ``SU(2)``: the Kronecker quiver (two nodes, two arrows,
``B = [0\\ 2;\\ -2\\ 0]``). It is affine ``\\tilde A_1`` - mutation-finite but **not**
finite-type (cluster-infinite: the monopole/dyon nodes mutate to the full dyon tower).
"""
su2_bps_quiver() = ClusterAlgebras.Quiver([0 2; -2 0])

# ── BPS spectrum ─────────────────────────────────────────────────────────────────────

# (abstract charge (a,b), Ω) list for a chamber, phase-ordered by slope b/a
function _su2_rays(chamber::Symbol, tower::Integer)
    if chamber === :strong
        return [((1, 0), 1), ((0, 1), 1)]
    elseif chamber === :weak
        rays = Tuple{Tuple{Int,Int},Int}[]
        for n in 0:tower
            push!(rays, ((n + 1, n), 1))          # preprojective P_n
        end
        push!(rays, ((1, 1), -2))                 # W-boson δ, vector multiplet
        for n in 0:tower
            push!(rays, ((n, n + 1), 1))          # preinjective I_n
        end
        # sort by slope (angle of (a,b)); a = 0 (pure γ₂) sits at the top
        sort!(rays; by = r -> atan(r[1][2], r[1][1]))
        return rays
    else
        throw(Resurgence.InvalidArgument("chamber must be :strong or :weak, got :$chamber"))
    end
end

# the chamber's states from an already-computed period pair, so a caller that also needs
# (a_D, a) itself - the twistor layer - shares one `_ns_periods` call instead of rebuilding
# the NS operator per charge
function _su2_states(aD::Complex{F}, a::Complex{F}, resolved::Symbol,
                     tower::Integer) where {F}
    states = BPSState{F}[]
    for ((p, q), Ω) in _su2_rays(resolved, tower)
        γ = _phys_charge(p, q)
        push!(states, BPSState([γ[1], γ[2]], γ[1] * aD + γ[2] * a, Ω))
    end
    states
end

"""
    su2_bps_states(sw::SeibergWittenSU2, u::Number; chamber = :auto, tower = 4,
                   ħ = 0, order = 0) -> Vector{BPSState}

The BPS spectrum of pure ``SU(2)`` at Coulomb modulus `u`, phase-ordered. The default
`chamber = :auto` resolves the chamber from the wall of marginal stability via
[`sw_chamber`](@ref); `chamber = :strong` returns the two states (monopole `(1,0)`,
dyon `(−1,2)`); `chamber = :weak` returns the W-boson `(0,2)` (with `omega = −2`) and
the preprojective/preinjective dyon towers truncated at `tower`. Each
[`BPSState`](@ref) carries its physical charge `(n_m,n_e)`, its central charge
`Z_γ = n_m a_D + n_e a` evaluated at `u` anywhere on the `u`-plane
([`central_charge`](@ref)), and its DT invariant `omega`.

With `order = m ≥ 1` the central charges carry their Nekrasov–Shatashvili corrections
through `ħ^{2m}` ([`quantum_sw_periods`](@ref)), and `chamber = :auto` resolves against
the **deformed** wall at the same `ħ`. Charges and `omega` values are unchanged - the
deformation acts on the periods, not on the BPS quiver.
"""
function su2_bps_states(sw::SeibergWittenSU2, u::Number; chamber::Symbol = :auto,
                        tower::Integer = 4, ħ = 0, order::Integer = 0)
    aD, a = _ns_periods(sw, u, ħ, order, _sw_float(sw, u))
    resolved = chamber === :auto ? sw_chamber(sw, u; ħ, order) : chamber
    _su2_states(aD, a, resolved, tower)
end

# ── classical Kontsevich–Soibelman wall-crossing (self-validating) ─────────────────────

# skew form of the abstract Kronecker torus: Λ = [0 2; -2 0], ⟨g,m⟩ = 2(g₁m₂ − g₂m₁) (always
# even, so the quadratic-refinement sign is trivial).
_ks_skew(g::NTuple{2,Int}, m::Vector{Int}) = 2 * (g[1] * m[2] - g[2] * m[1])

# (1 − X_γ)^e as a truncated series {exponent vector => coeff}, total degree ≤ D. e any integer;
# generalized binomial ∏(e−i)/k! makes negative e (the vector factor) a geometric-type series.
function _one_minus_pow(γ::NTuple{2,Int}, e::Int, D::Int)
    out = Dict{Vector{Int},Rational{BigInt}}(zeros(Int, 2) => one(Rational{BigInt}))
    dγ = γ[1] + γ[2]
    dγ == 0 && return out
    coeff = one(Rational{BigInt})
    for k in 1:(D ÷ dγ)
        coeff *= Rational{BigInt}(e - (k - 1), k)            # binomial(e,k)/binomial(e,k-1)
        out[[k * γ[1], k * γ[2]]] = coeff * (-1)^k           # (−X_γ)^k
    end
    out
end

# apply the ray automorphism K_γ^Ω to a truncated series: X_μ ↦ X_μ (1 − X_γ)^{Ω⟨γ,μ⟩}
function _apply_ray(series::Dict{Vector{Int},Rational{BigInt}}, γ::NTuple{2,Int}, Ω::Int, D::Int)
    out = Dict{Vector{Int},Rational{BigInt}}()
    for (ν, c) in series
        factor = _one_minus_pow(γ, Ω * _ks_skew(γ, ν), D)
        for (kv, fc) in factor
            μ = ν + kv
            (μ[1] + μ[2]) ≤ D || continue
            newc = get(out, μ, zero(Rational{BigInt})) + c * fc
            iszero(newc) ? delete!(out, μ) : (out[μ] = newc)
        end
    end
    out
end

# the composite ray-automorphism applied to the generator X_{seed}, rays acting in phase
# (slope) order - the orientation in which the strong product K_{γ₂}K_{γ₁} matches the weak
# tower product (pinned by requiring closure of the identity).
function _compose(rays, seed::NTuple{2,Int}, D::Int)
    v = Dict{Vector{Int},Rational{BigInt}}([seed[1], seed[2]] => one(Rational{BigInt}))
    for ((a, b), Ω) in rays
        v = _apply_ray(v, (a, b), Ω, D)
    end
    v
end

"""
    verify_su2_wall_crossing(; degree = 4) -> NamedTuple

Verify the pure-``SU(2)`` Kontsevich–Soibelman wall-crossing identity to total ``ŷ``-degree
`degree`: the strong-chamber ray product ``K_{γ₂}K_{γ₁}`` equals the phase-ordered weak-chamber
product over the dyon towers and the W-boson ``K_δ^{-2}``, as automorphisms of the Kronecker
torus. Both sides are applied to the generators and compared as exact truncated power series.

Returns `(; closed, max_residual)`: `closed` is `true` iff the two sides agree exactly to
`degree` (they do - this is the ``m=2`` analog of the ``A_2`` pentagon). The check is
self-validating: closure pins the spectrum, the Ω values (including the vector multiplet's
`−2`), and the ordering.
"""
function verify_su2_wall_crossing(; degree::Integer = 4)
    degree ≥ 1 || throw(Resurgence.InvalidArgument("degree must be ≥ 1, got $degree"))
    D = Int(degree)
    tower = D ÷ 2 + 1                              # enough dyons to saturate degree D
    strong = [((0, 1), 1), ((1, 0), 1)]           # K_{γ₂} K_{γ₁}
    weak = _su2_rays(:weak, tower)
    max_res = 0.0
    for seed in ((1, 0), (0, 1))
        ls = _compose(strong, seed, D)
        rs = _compose(weak, seed, D)
        for μ in union(keys(ls), keys(rs))
            d = get(ls, μ, zero(Rational{BigInt})) - get(rs, μ, zero(Rational{BigInt}))
            max_res = max(max_res, abs(float(d)))
        end
    end
    (; closed = iszero(max_res), max_residual = max_res)
end

# ── refined (motivic) Kontsevich–Soibelman wall-crossing ──────────────────────────────
#
# The same identity in the quantum torus, where each state contributes quantum
# dilogarithm factors rather than a Poisson ray automorphism. Everything about the
# torus is `src/refined.jl`'s ledger; what is new here is the *affine* case, i.e. a
# state whose refined index is not 1.
#
# ── Ledger (pinned by `verify_su2_refined_wall_crossing`, self-validating) ────────────
# R1. THE TORUS. Λ = −B of `su2_bps_quiver()` = [0 −2; 2 0], i.e. refined.jl's ledger
#     item 1 applied to the Kronecker quiver - no new convention. Factor order is the
#     ray order in the abstract basis: the strong side is 𝔼(ŷ^{γ₂})𝔼(ŷ^{γ₁}), the weak
#     side runs (1,0), (2,1), (3,2), … , δ , … , (2,3), (1,2), (0,1) - preprojective
#     dyons by increasing n, then the W-boson, then preinjective dyons by decreasing n.
# R2. THE W-BOSON FACTOR WAS DERIVED, NOT ASSERTED. Dividing the strong product by the
#     two dyon towers leaves a residue supported on multiples of δ alone, and that
#     residue is
#
#         𝔼(−v ŷ^δ)^{-1} · 𝔼(−v^{-1} ŷ^δ)^{-1}
#
#     verified to ŷ-degree 8. This is Ω(δ, y) = −y − y^{-1} (a vector multiplet, [DG09])
#     *together with the quadratic-refinement sign* σ(δ) = −1: the argument carries a
#     minus, exactly as the classical layer's ray automorphism uses (1 − X_γ). Without
#     the sign nothing closes - the candidates 𝔼(v^{±1}ŷ^δ)^{-1}, 𝔼(ŷ^δ)^{-2} and
#     𝔼(v^{±2}ŷ^δ)^{-1} were each measured and each fails at order ŷ^δ. The two
#     sub-factors commute (δ is isotropic), so their order is not a convention.
# R3. CLASSICAL LIMIT. Summing the exponents of the W-boson's factors gives −2 = the
#     classical Ω(δ) of `_su2_rays`, i.e. the refined index at y → 1.

# One BPS state's refined contribution: the charge and a list of (shift, sign, exponent)
# triples, the factor being ∏ 𝔼(σ v^{shift} ŷ^γ)^{exponent}.
const _RefinedRay = Tuple{Tuple{Int,Int},Vector{NTuple{3,Int}}}

_hypermultiplet() = [(0, 1, 1)]                       # 𝔼(ŷ^γ), Ω = 1
_vector_multiplet() = [(1, -1, -1), (-1, -1, -1)]     # 𝔼(−vŷ^δ)⁻¹𝔼(−v⁻¹ŷ^δ)⁻¹, Ω = −y−y⁻¹

"""
    su2_refined_rays(chamber::Symbol, tower::Integer) -> Vector

The refined BPS rays of pure ``SU(2)`` in the abstract Kronecker basis: one entry
`(γ, factors)` per state, in ray order, where `factors` is a list of
`(shift, sign, exponent)` triples and the state contributes
``∏ 𝔼(σ v^{shift} ŷ^γ)^{exponent}``.

A hypermultiplet is the single factor `𝔼(ŷ^γ)` (refined index `1`); the W-boson
`δ = (1,1)` is `𝔼(−v ŷ^δ)^{-1} 𝔼(−v^{-1} ŷ^δ)^{-1}`, the refined index
`Ω(δ, y) = −y − y^{-1}` with the quadratic-refinement sign (ledger R2). The dyon
towers are truncated at `tower`.
"""
function su2_refined_rays(chamber::Symbol, tower::Integer)
    tower ≥ 0 || throw(Resurgence.InvalidArgument("tower must be ≥ 0, got $tower"))
    if chamber === :strong
        return _RefinedRay[((0, 1), _hypermultiplet()), ((1, 0), _hypermultiplet())]
    elseif chamber === :weak
        rays = _RefinedRay[]
        for n in 0:tower
            push!(rays, ((n + 1, n), _hypermultiplet()))       # preprojective, slope ↗ 1
        end
        push!(rays, ((1, 1), _vector_multiplet()))             # the W-boson
        for n in tower:-1:0
            push!(rays, ((n, n + 1), _hypermultiplet()))       # preinjective, slope 1 ↗ ∞
        end
        return rays
    end
    throw(Resurgence.InvalidArgument("chamber must be :strong or :weak, got :$chamber"))
end

# The ordered dilogarithm product of a refined ray list, in the Kronecker torus.
function _su2_refined_product(rays::Vector{_RefinedRay}, D::Int)
    charges, shifts, signs, exps = Vector{Int}[], Int[], Int[], Int[]
    for (γ, factors) in rays, (s, σ, e) in factors
        push!(charges, [γ[1], γ[2]])
        push!(shifts, s)
        push!(signs, σ)
        push!(exps, e)
    end
    n = length(charges)
    w = ClusterAlgebras.QuantumDilogWord(charges, _su2_refined_skew(), ones(Int, n),
                                         ones(Int, n), shifts, exps, signs)
    ClusterAlgebras.ks_dilog_product(w; truncation_degree = D)
end

# Λ = −B of the Kronecker BPS quiver (refined.jl ledger item 1, ledger R1 here).
_su2_refined_skew() = -su2_bps_quiver().B[1:2, 1:2]

"""
    verify_su2_refined_wall_crossing(; degree = 4) -> NamedTuple

Verify the **refined (motivic)** pure-``SU(2)`` wall-crossing identity to total
``ŷ``-degree `degree`: the strong-chamber quantum dilogarithm product
``𝔼(ŷ^{γ₂}) 𝔼(ŷ^{γ₁})`` equals the ray-ordered weak-chamber product over the two dyon
towers and the W-boson, in the Kronecker quantum torus.

Returns `(; closed, vector_factor, residue)`: `closed` says the two sides agree
exactly as truncated series over `Frac(ℤ[v])`, `residue` is the strong product
divided by the two dyon towers - the W-boson factor **as measured**, which
`vector_factor` reproduces from `Ω(δ, y) = −y − y^{-1}` with the quadratic-refinement
sign (ledger R2). At `v → 1` the exponents sum to the classical `Ω(δ) = −2` of
[`verify_su2_wall_crossing`](@ref).
"""
function verify_su2_refined_wall_crossing(; degree::Integer = 4)
    degree ≥ 1 || throw(Resurgence.InvalidArgument("degree must be ≥ 1, got $degree"))
    D = Int(degree)
    tower = D                                    # saturates: P_n has ŷ-degree 2n+1
    strong = _su2_refined_product(su2_refined_rays(:strong, tower), D)
    weak = _su2_refined_product(su2_refined_rays(:weak, tower), D)

    # the measured residue: strong ÷ (both towers), which must be supported on δ alone
    Λ = _su2_refined_skew()
    _, v = ClusterAlgebras._ks_ring(D)
    left = _su2_refined_product(_RefinedRay[((n + 1, n), _hypermultiplet())
                                            for n in 0:tower], D)
    right = _su2_refined_product(_RefinedRay[((n, n + 1), _hypermultiplet())
                                             for n in tower:-1:0], D)
    residue = ClusterAlgebras._qt_mul(
        ClusterAlgebras._qt_mul(ClusterAlgebras._qt_inv(left, Λ, v, D), strong, Λ, v, D),
        ClusterAlgebras._qt_inv(right, Λ, v, D), Λ, v, D)
    vector = _su2_refined_product(_RefinedRay[((1, 1), _vector_multiplet())], D)

    (; closed = strong == weak, vector_factor = vector, residue,
       delta_only = all(γ[1] == γ[2] for γ in keys(residue)))
end
