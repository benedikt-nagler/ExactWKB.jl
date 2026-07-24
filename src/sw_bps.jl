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

"""
    su2_bps_states(sw::SeibergWittenSU2, u::Number; chamber = :auto, tower = 4)
        -> Vector{BPSState}

The BPS spectrum of pure ``SU(2)`` at Coulomb modulus `u`, phase-ordered. The default
`chamber = :auto` resolves the chamber from the wall of marginal stability via
[`sw_chamber`](@ref); `chamber = :strong` returns the two states (monopole `(1,0)`,
dyon `(−1,2)`); `chamber = :weak` returns the W-boson `(0,2)` (with `omega = −2`) and
the preprojective/preinjective dyon towers truncated at `tower`. Each
[`BPSState`](@ref) carries its physical charge `(n_m,n_e)`, its central charge
`Z_γ = n_m a_D + n_e a` evaluated at `u` anywhere on the `u`-plane
([`central_charge`](@ref)), and its DT invariant `omega`.
"""
function su2_bps_states(sw::SeibergWittenSU2, u::Number; chamber::Symbol = :auto,
                        tower::Integer = 4)
    resolved = chamber === :auto ? sw_chamber(sw, u) : chamber
    rays = _su2_rays(resolved, tower)
    states = BPSState{typeof(float(real(u)))}[]
    for ((a, b), Ω) in rays
        γ = _phys_charge(a, b)
        Z = central_charge(sw, u, γ)
        push!(states, BPSState([γ[1], γ[2]], Z, Ω))
    end
    states
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
