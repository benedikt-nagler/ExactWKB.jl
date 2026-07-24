# Algebraic-geometry periods of the pure SU(2) Seiberg–Witten curve - the first rung of
# the SU(2) north-star ascent. FIREWALLED: this file imports only `Resurgence`
# (FormalSeries, for the ħ-series quantum periods) and `QuadGK`; it must NOT touch the
# polynomial WKB engine (`SchrodingerProblem`, the Riccati ring, `turning_points`) or
# `ClusterAlgebras`. That isolation is deliberate - the module is extraction-ready as a
# future standalone `Periods.jl` foundation package once a second consumer appears.
#
# ── Convention ledger (Λ = the dynamical scale; hand-pinned, oracle-verified) ──────────
# The Nekrasov–Shatashvili limit of pure SU(2) is the Mathieu problem, which in the
# package's own Iwaki–Nakanishi normalization ħ²ψ″ = Q(x)ψ, Q = V − E reads
#     V(x) = 2Λ² cos 2x,     E = u    (the Coulomb modulus).
# The two SW periods are the classical action integrals of that problem:
#   • electric  a(u)   = (1/2π) ∮ √(u − 2Λ² cos 2x) dx   over one real cycle x ∈ [0,2π];
#   • magnetic  a_D(u) = 4 · (i/2π) ∫_{−y_t}^{y_t} √(u − 2Λ² cosh 2y) dy,
#                        y_t = ½ arccosh(u/2Λ²)  (the tunnelling cycle through the barrier).
# a's normalization is fixed by the exact weak-coupling identity (verified in the tests)
#     u(a) = a² + Λ⁴/(2a²) + 5Λ⁸/(32a⁶) + 9Λ¹²/(64a¹⁰) + …
# a_D's factor 4 relative to the naive tunnelling action (rescaled 2026-07-24, rung 3 -
# nothing in rungs 1–2 was sensitive to it) is pinned by three strong-coupling oracles:
# the massless state at u = −2Λ² must carry a *lattice* dyon charge (a_D/a → −2 there,
# approached through Im u > 0, i.e. ±(1,2) massless - the (1,−2)-dyon label transported
# through the upper half-plane); the Picard–Fuchs monodromy must be integral in the
# (a_D, a) basis; and Im τ = Im(da_D/da) > 0 (metric positivity, which kills the
# opposite sign).
# The curve degenerates (a BPS state becomes massless) at u = ±2Λ²: monopole at +2Λ²
# (a_D → 0), dyon at −2Λ². This rung implements the weak-coupling / electric region
# u > 2Λ², where both period integrals are real one-dimensional quadratures.
#
# The NS quantum periods are the all-orders exact-WKB corrections. For ħ²ψ″ = Qψ the
# ħ²-order period integrand is  S₂ = Q″/(8Q^{3/2}) − 5Q′²/(32Q^{5/2})  (the standard
# Dunham term); on the electric cycle (no turning points) it integrates to a finite
# correction  a₂(u), leading behaviour a₂ ≈ −Λ⁴/(4u^{5/2}) at large u.
#
# The magnetic cycle DOES enclose turning points, where S₂ ~ Q^{−5/2} is not integrable
# pointwise - the naive quadrature diverges. But it needs no Voros regularization: split
#   S₂ = (1/8) d/dx[Q′/Q^{3/2}] + (1/32) Q′²/Q^{5/2},
# the total derivative drops over any cycle, and Q′² = 16Λ⁴ − 4(Q+u)² reduces the rest to
# ∮Q^{−1/2}, ∮Q^{−3/2}, ∮Q^{−5/2} = (−2∂_u, −4∂_u², −(8/3)∂_u³) of the classical action
# ∮Q^{1/2}. Collecting and eliminating the third derivative with the Picard–Fuchs relation
# (u²−4Λ⁴)Π″ + Π/4 = 0 gives a universal, cycle-independent first-order operator on the
# classical period Π(u):
#   a₂(Π) = (1/6)·Π′(u) − u·Π(u)/(12(u²−4Λ⁴)).
# Applied to the electric period this reproduces the quadrature above and its −Λ⁴/(4u^{5/2})
# tail exactly; applied to a_D it IS the magnetic quantum period. Valid on the whole
# u-plane (both periods carry their ħ² correction), singular only at u = ±2Λ².
#
# ── Rung 3: complex u (closed forms), Picard–Fuchs continuation, the MS wall ─────────
# The classical periods have closed forms in complete elliptic integrals (parameter
# convention m, i.e. K(m) = ∫₀^{π/2} dθ/√(1 − m sin²θ)):
#     a(u)   = (2/π)·√(u+2Λ²)·E(m̃),                      m̃ = 4Λ²/(u+2Λ²),
#     a_D(u) = −(8iΛ/π)·[E(m) − (1−m)K(m)],              m  = (2Λ²−u)/(4Λ²),
#     da/du  = K(m̃)/(π√(u+2Λ²)),      da_D/du = i·K(m)/(πΛ),
# pinned against the direct quadratures on the real weak-coupling axis (tests). K, E
# are evaluated for complex parameter via in-file Carlson symmetric forms R_F, R_D
# (duplication algorithm, DLMF §19.36, principal square roots). Build-vs-reuse
# decision (2026-07-24): built in-file - Elliptic.jl is real-only, and
# EllipticFunctions.jl (MIT) would pull SpecialFunctions and break this file's
# firewall / Periods.jl extraction path.
#
# Branches: with principal roots, a(u) is analytic on ℂ ∖ (−∞, 2Λ²] and a_D(u) on
# ℂ ∖ (−∞, −2Λ²] (cuts on the real axis; on a cut the value is a one-sided limit).
# The dyon point u = −2Λ² is a branch point of a and is guarded. Genuine analytic
# continuation across the cuts is `continue_periods`: RK4 transport of the
# Picard–Fuchs equation  (u² − 4Λ⁴)Π″ + Π/4 = 0  along an arbitrary u-path. Loops
# around u = ±2Λ², ∞ reproduce the exact integer `sw_monodromy` matrices - the oracle
# that pins every branch and sign above. The wall of marginal stability
# {u : a_D/a ∈ ℝ} is traced by `ms_wall` (bisection along rays from the origin;
# endpoints are the singularities ±2Λ²), and `sw_chamber` classifies a point as
# :strong (inside) or :weak (outside; points on the wall or the real axis beyond
# ±2Λ² count as :weak).

using QuadGK: quadgk

"""
    SeibergWittenSU2{T}

The pure ``SU(2)`` ``\\mathcal N = 2`` Seiberg–Witten geometry, parametrized by its
dynamical scale `Λ` (field `Λ::T`, default `1`). Everything hangs off the `u`-plane:
periods [`sw_periods`](@ref) / [`quantum_sw_periods`](@ref), the singular moduli
[`sw_singularities`](@ref), the [`central_charge`](@ref), and the monodromy
[`sw_monodromy`](@ref). Immutable.

Convention: the Nekrasov–Shatashvili limit is the Mathieu problem
``ħ²ψ″ = (2Λ²\\cos 2x − u)ψ``; see the ledger at the top of `src/sw_curve.jl`.

# Examples
```julia
sw = SeibergWittenSU2()            # Λ = 1
sw = SeibergWittenSU2(Λ = 2)
sw_periods(sw, 10.0)               # (a, a_D) at u = 10
```
"""
struct SeibergWittenSU2{T}
    Λ::T
end

SeibergWittenSU2(; Λ = 1) = SeibergWittenSU2(float(Λ))

"""
    dynamical_scale(sw::SeibergWittenSU2) -> Number

The dynamical scale ``Λ``.
"""
dynamical_scale(sw::SeibergWittenSU2) = sw.Λ

"""
    sw_singularities(sw::SeibergWittenSU2) -> NamedTuple

The singular moduli of the ``u``-plane where a BPS state becomes massless: `monopole`
at ``u = +2Λ²`` (where ``a_D = 0``) and `dyon` at ``u = −2Λ²`` (plus the weak-coupling
puncture at ``∞``).
"""
function sw_singularities(sw::SeibergWittenSU2)
    (monopole = 2 * sw.Λ^2, dyon = -2 * sw.Λ^2)
end

# working float type for the numerics
_sw_float(sw::SeibergWittenSU2{T}, u) where {T} = float(promote_type(T, real(typeof(u))))

# ── Carlson symmetric elliptic integrals (complex-capable, in-file) ────────────────

# R_F(x,y,z) by the duplication theorem (DLMF §19.36); principal square roots fix the
# branches. Converges for all arguments off the negative real axis.
function _carlson_rf(x::Number, y::Number, z::Number)
    xt, yt, zt = promote(complex(float(x)), complex(float(y)), complex(float(z)))
    F = real(typeof(xt))
    errtol = (eps(F) / 3)^(1 // 6)
    for _ in 1:500
        sx, sy, sz = sqrt(xt), sqrt(yt), sqrt(zt)
        λ = sx * sy + sy * sz + sz * sx
        xt, yt, zt = (xt + λ) / 4, (yt + λ) / 4, (zt + λ) / 4
        ave = (xt + yt + zt) / 3
        iszero(ave) && throw(PeriodError("Carlson R_F: degenerate arguments"))
        dx, dy, dz = (ave - xt) / ave, (ave - yt) / ave, (ave - zt) / ave
        if max(abs(dx), abs(dy), abs(dz)) < errtol
            e2 = dx * dy - dz^2
            e3 = dx * dy * dz
            return (1 + (e2 / 24 - F(3) / 44 * e3 - F(1) / 10) * e2 + e3 / 14) / sqrt(ave)
        end
    end
    throw(PeriodError("Carlson R_F did not converge (argument at/near the K(m=1) divergence)"))
end

# R_D(x,y,z) = R_J(x,y,z,z), same duplication scheme with an accumulated sum
function _carlson_rd(x::Number, y::Number, z::Number)
    xt, yt, zt = promote(complex(float(x)), complex(float(y)), complex(float(z)))
    F = real(typeof(xt))
    errtol = (eps(F) / 3)^(1 // 6)
    acc = zero(xt)
    fac = one(F)
    for _ in 1:500
        sx, sy, sz = sqrt(xt), sqrt(yt), sqrt(zt)
        λ = sx * sy + sy * sz + sz * sx
        acc += fac / (sz * (zt + λ))
        fac /= 4
        xt, yt, zt = (xt + λ) / 4, (yt + λ) / 4, (zt + λ) / 4
        ave = (xt + yt + 3zt) / 5
        iszero(ave) && throw(PeriodError("Carlson R_D: degenerate arguments"))
        dx, dy, dz = (ave - xt) / ave, (ave - yt) / ave, (ave - zt) / ave
        if max(abs(dx), abs(dy), abs(dz)) < errtol
            ea = dx * dy
            eb = dz^2
            ec = ea - eb
            ed = ea - 6eb
            ee = ed + 2ec
            return 3acc + fac * (1 + ed * (-F(3) / 14 + F(9) / 88 * ed - F(9) / 52 * dz * ee) +
                                 dz * (ee / 6 + dz * (-F(9) / 22 * ec + dz * F(3) / 26 * ea))) /
                          (ave * sqrt(ave))
        end
    end
    throw(PeriodError("Carlson R_D did not converge"))
end

# complete elliptic integrals in the parameter convention, K(m) = ∫₀^{π/2} dθ/√(1−m sin²θ)
_ellipK(m::Number) = _carlson_rf(zero(m), 1 - m, one(m))
function _ellipE(m::Number)
    m == 1 && return one(complex(float(m)))     # E(1) = 1; the R_F − R_D split is ∞ − ∞ there
    _carlson_rf(zero(m), 1 - m, one(m)) - m / 3 * _carlson_rd(zero(m), 1 - m, one(m))
end

# ── classical periods: closed forms on the whole u-plane ───────────────────────────

# guard the branch point u = −2Λ² (where m̃ = 4Λ²/(u+2Λ²) blows up)
function _check_dyon_point(sw::SeibergWittenSU2, U, Λ)
    iszero(U + 2Λ^2) && throw(PeriodError(
        "u = −2Λ² = $(-2*sw.Λ^2) is the dyon branch point; evaluate nearby, or use " *
        "continue_periods along a path"))
    nothing
end

# a(u) = (2/π)·√(u+2Λ²)·E(m̃), m̃ = 4Λ²/(u+2Λ²); principal branch, cut (−∞, 2Λ²] on ℝ
function _electric_a_closed(sw::SeibergWittenSU2, u, F)
    Λ = F(sw.Λ)
    U = Complex{F}(u)
    _check_dyon_point(sw, U, Λ)
    2 / F(π) * sqrt(U + 2Λ^2) * _ellipE(4Λ^2 / (U + 2Λ^2))
end

# a_D(u) = −(8iΛ/π)·[E(m) − (1−m)K(m)], m = (2Λ²−u)/(4Λ²); cut (−∞, −2Λ²] on ℝ
function _magnetic_aD_closed(sw::SeibergWittenSU2, u, F)
    Λ = F(sw.Λ)
    m = (2Λ^2 - Complex{F}(u)) / (4Λ^2)
    -(8im * Λ / F(π)) * (_ellipE(m) - (1 - m) * _ellipK(m))
end

# derivatives (initial data for the Picard–Fuchs transport):
# da/du = K(m̃)/(π√(u+2Λ²)),  da_D/du = i·K(m)/(4πΛ)
function _da_du_closed(sw::SeibergWittenSU2, u, F)
    Λ = F(sw.Λ)
    U = Complex{F}(u)
    _check_dyon_point(sw, U, Λ)
    _ellipK(4Λ^2 / (U + 2Λ^2)) / (F(π) * sqrt(U + 2Λ^2))
end

function _daD_du_closed(sw::SeibergWittenSU2, u, F)
    Λ = F(sw.Λ)
    im * _ellipK((2Λ^2 - Complex{F}(u)) / (4Λ^2)) / (F(π) * Λ)
end

# ── classical periods: direct quadratures (weak real axis; kept as the pinning oracle
#    for the closed forms and for the ħ² correction) ─────────────────────────────────

# electric period a(u) = (1/2π) ∫₀^{2π} √(u − 2Λ² cos 2x) dx
function _electric_a(sw::SeibergWittenSU2, u, F)
    Λ2 = F(sw.Λ)^2
    U = F(u)
    integrand(x) = sqrt(U - 2Λ2 * cos(2x))
    val, _ = quadgk(integrand, F(0), F(2π); rtol = sqrt(eps(F)))
    val / (2 * F(π))
end

# magnetic period a_D(u) = 4·(i/2π) ∫_{−y_t}^{y_t} √(u − 2Λ² cosh 2y) dy (the lattice
# normalization factor 4 is pinned in rung 3 - see the ledger)
function _magnetic_aD(sw::SeibergWittenSU2, u, F)
    Λ2 = F(sw.Λ)^2
    U = F(u)
    yt = acosh(U / (2Λ2)) / 2
    integrand(y) = sqrt(U - 2Λ2 * cosh(2y))
    # integrable √-singularities at ±y_t; split at 0 so quadgk brackets both endpoints
    val, _ = quadgk(integrand, -yt, zero(F), yt; rtol = sqrt(eps(F)))
    complex(zero(F), 2 * val / F(π))
end

"""
    sw_periods(sw::SeibergWittenSU2, u::Number) -> NamedTuple

The classical Seiberg–Witten periods at Coulomb modulus `u`, valid on the **whole**
`u`-plane: the electric period `a` (real for real `u > 2Λ²`) and the magnetic/dual
period `a_D` (imaginary there). Evaluated in closed form via complete elliptic
integrals (in-file Carlson symmetric forms), principal branch - cuts on the real axis,
`(−∞, 2Λ²]` for `a` and `(−∞, −2Λ²]` for `a_D`; on a cut the value is a one-sided
limit, and genuine analytic continuation across the cuts is [`continue_periods`](@ref).
The normalization is pinned by the weak-coupling identity `u = a² + Λ⁴/(2a²) + …` and
by direct Mathieu-action quadrature (see the ledger in `src/sw_curve.jl`). Throws
[`PeriodError`](@ref) at the dyon branch point `u = −2Λ²` exactly.
"""
function sw_periods(sw::SeibergWittenSU2, u::Number)
    F = _sw_float(sw, u)
    (a = _electric_a_closed(sw, u, F), a_D = _magnetic_aD_closed(sw, u, F))
end

# ── quantum (Nekrasov–Shatashvili) periods ─────────────────────────────────────────

# the ħ²-order correction to a classical period Π(u) with u-derivative dΠ, via the
# universal Picard–Fuchs-reduced operator (see the ledger at the top of the file):
#   a₂ = (1/6)·Π′ − u·Π/(12(u²−4Λ⁴)).
# Cycle-independent (electric ⟵ Π = a, magnetic ⟵ Π = a_D), whole u-plane.
function _period_a2(sw::SeibergWittenSU2, u, Π, dΠ, F)
    U = Complex{F}(u)
    dΠ / 6 - U * Π / (12 * (U^2 - 4 * F(sw.Λ)^4))
end

_electric_a2(sw::SeibergWittenSU2, u, F) =
    _period_a2(sw, u, _electric_a_closed(sw, u, F), _da_du_closed(sw, u, F), F)

_magnetic_aD2(sw::SeibergWittenSU2, u, F) =
    _period_a2(sw, u, _magnetic_aD_closed(sw, u, F), _daD_du_closed(sw, u, F), F)

# independent oracle for the electric correction: a₂ = −(1/2π) ∫₀^{2π} P₂ dx with
# P₂ (in the real barrier variable W = u − 2Λ²cos2x = −Q) equal to
#   W″/(8W^{3/2}) − 5W′²/(32W^{5/2}),   W′ = 4Λ² sin 2x,  W″ = 8Λ² cos 2x.
# The electric cycle has no turning points, so this quadrature converges; it pins the
# operator form above on real u > 2Λ². (The magnetic cycle has no such quadrature.)
function _electric_a2_quad(sw::SeibergWittenSU2, u, F)
    Λ2 = F(sw.Λ)^2
    U = F(u)
    function p2(x)
        W = U - 2Λ2 * cos(2x)
        Wp = 4Λ2 * sin(2x)
        Wpp = 8Λ2 * cos(2x)
        Wpp / (8 * W^F(1.5)) - 5 * Wp^2 / (32 * W^F(2.5))
    end
    val, _ = quadgk(p2, F(0), F(2π); rtol = sqrt(eps(F)))
    -val / (2 * F(π))
end

# the ħ² correction is singular where the Picard–Fuchs coefficient u²−4Λ⁴ vanishes
function _check_h2_singular(sw::SeibergWittenSU2, u)
    Λ2 = sw.Λ^2
    (u == 2Λ2 || u == -2Λ2) && throw(PeriodError(
        "u = $u is a branch point u = ±2Λ² = ±$(2Λ2); the ħ² quantum correction is " *
        "singular there (a period pinches). Evaluate nearby instead."))
    nothing
end

"""
    quantum_sw_periods(sw::SeibergWittenSU2, u::Number; order::Integer = 1) -> NamedTuple

The Nekrasov–Shatashvili quantum periods at modulus `u`, as `Resurgence.FormalSeries` in
the symbol `ħ`. `order` counts powers of `ħ²`: `order = 0` returns the classical periods
as constant series, `order = 1` adds the leading exact-WKB (`ħ²`) correction to **both**
periods, `a = a₀ + a₂·ħ² + O(ħ⁴)` (with `a₂ ≈ −Λ⁴/(4u^{5/2})` at large `u`) and
`a_D = a_{D,0} + a_{D,2}·ħ² + O(ħ⁴)`. Both corrections come from the same
Picard–Fuchs-reduced differential operator `a₂ = (1/6)Π′ − uΠ/(12(u²−4Λ⁴))` - no
turning-point (Voros) regularization is required for the magnetic cycle - so the result
is valid on the **whole** `u`-plane (principal branch, see [`sw_periods`](@ref)).

Throws [`PeriodError`](@ref) for `order < 0` or `order > 1`, and for `order = 1` at the
branch points `u = ±2Λ²`, where the `ħ²` correction is singular.
"""
function quantum_sw_periods(sw::SeibergWittenSU2, u::Number; order::Integer = 1)
    0 ≤ order ≤ 1 || throw(PeriodError(
        "order must be 0 or 1 in this rung (got $order); higher exact-WKB orders are a " *
        "later rung"))
    order == 1 && _check_h2_singular(sw, u)
    F = _sw_float(sw, u)
    a0 = _electric_a_closed(sw, u, F)
    aD0 = _magnetic_aD_closed(sw, u, F)
    if order == 0
        a_coeffs = [complex(a0)]
        aD_coeffs = [complex(aD0)]
    else
        z = zero(complex(a0))
        a_coeffs = [complex(a0), z, complex(_electric_a2(sw, u, F))]
        aD_coeffs = [complex(aD0), z, complex(_magnetic_aD2(sw, u, F))]
    end
    a_series = Resurgence.FormalSeries(a_coeffs, :ħ; power_offset = 0 // 1)
    aD_series = Resurgence.FormalSeries(aD_coeffs, :ħ; power_offset = 0 // 1)
    (a = a_series, a_D = aD_series)
end

"""
    central_charge(sw::SeibergWittenSU2, u, charge::NTuple{2,Int}; ħ = 0, order = 0)
        -> Number

The central charge ``Z_{(n_m,n_e)} = n_m\\,a_D + n_e\\,a`` of the BPS charge
`charge = (n_m, n_e)` at modulus `u`, valid on the whole `u`-plane (closed-form
periods, principal branch - see [`sw_periods`](@ref)). With `order = 0` (default) this
is the classical value; with `order = 1` **both** periods carry their `ħ²` quantum
correction (the Picard–Fuchs-reduced operator, valid on the whole `u`-plane), evaluated
at the given `ħ`. The BPS mass is `abs(central_charge(...))`. This is the seam the
spectrum/wall layer consumes.
"""
function central_charge(sw::SeibergWittenSU2, u, charge::NTuple{2,Int}; ħ = 0, order::Integer = 0)
    0 ≤ order ≤ 1 || throw(PeriodError(
        "order must be 0 or 1 (got $order)"))
    F = _sw_float(sw, u)
    nm, ne = charge
    a = _electric_a_closed(sw, u, F)
    aD = _magnetic_aD_closed(sw, u, F)
    if order == 1
        _check_h2_singular(sw, u)
        a += complex(_electric_a2(sw, u, F)) * F(ħ)^2
        aD += complex(_magnetic_aD2(sw, u, F)) * F(ħ)^2
    end
    nm * aD + ne * a
end

# ── Picard–Fuchs continuation ──────────────────────────────────────────────────────

# right-hand side of the first-order system for Y = [Π_aD Π_a; Π_aD′ Π_a′], from
# (u² − 4Λ⁴)Π″ + Π/4 = 0
function _pf_rhs(Λ::F, u, Y) where {F}
    d = u^2 - 4Λ^4
    abs(d) < sqrt(eps(F)) * (1 + 4Λ^4) && throw(PeriodError(
        "Picard–Fuchs path passes through a singular point u = ±2Λ² (at u = $u)"))
    c = -1 / (4d)
    [Y[2, 1] Y[2, 2]; c * Y[1, 1] c * Y[1, 2]]
end

# transport Y along the straight segment u1 → u2: fixed-step RK4, step-doubling until
# the result is rtol-converged
function _pf_segment(Λ::F, u1, u2, Y, rtol) where {F}
    u1 == u2 && return Y
    prev = nothing
    n = 16
    while true
        Yc = Y
        h = (u2 - u1) / n
        for k in 1:n
            u = u1 + (k - 1) * h
            k1 = _pf_rhs(Λ, u, Yc)
            k2 = _pf_rhs(Λ, u + h / 2, Yc + h / 2 * k1)
            k3 = _pf_rhs(Λ, u + h / 2, Yc + h / 2 * k2)
            k4 = _pf_rhs(Λ, u + h, Yc + h * k3)
            Yc = Yc + h / 6 * (k1 + 2k2 + 2k3 + k4)
        end
        if prev !== nothing && maximum(abs.(Yc - prev)) ≤ rtol * maximum(abs.(Yc))
            return Yc
        end
        prev = Yc
        n *= 2
        n > 2^18 && throw(PeriodError(
            "Picard–Fuchs continuation did not converge on the segment $u1 → $u2 " *
            "(path too close to a singular point?)"))
    end
end

"""
    continue_periods(sw::SeibergWittenSU2, path::AbstractVector{<:Number};
                     rtol = 1e-10) -> NamedTuple

Analytic continuation of the classical periods along the piecewise-linear `u`-path
`path` (waypoints; a curved path - e.g. a monodromy circle - is passed as a dense
polyline). Initial values and derivatives at `path[1]` come from the closed forms;
the fundamental 2×2 system of the Picard–Fuchs equation `(u² − 4Λ⁴)Π″ + Π/4 = 0` is
then transported by adaptive RK4. Unlike [`sw_periods`](@ref) this is genuine
continuation - crossing the principal-branch cuts picks up the corresponding
[`sw_monodromy`](@ref) action (loops around `±2Λ²`/`∞` reproduce those exact integer
matrices; that oracle pins every branch convention in this file).

Returns `(a, a_D, da, da_D)` - the continued periods and their `u`-derivatives at
`path[end]`. A single-point path returns the closed-form data at that point. Throws
[`PeriodError`](@ref) if the path runs through a singular point.
"""
function continue_periods(sw::SeibergWittenSU2, path::AbstractVector{<:Number};
                          rtol::Real = 1e-10)
    isempty(path) && throw(Resurgence.InvalidArgument(
        "path must contain at least one point"))
    F = _sw_float(sw, first(path))
    Λ = F(sw.Λ)
    u0 = Complex{F}(first(path))
    Y = [_magnetic_aD_closed(sw, u0, F) _electric_a_closed(sw, u0, F);
         _daD_du_closed(sw, u0, F) _da_du_closed(sw, u0, F)]
    for i in 2:length(path)
        Y = _pf_segment(Λ, Complex{F}(path[i-1]), Complex{F}(path[i]), Y, F(rtol))
    end
    (a = Y[1, 2], a_D = Y[1, 1], da = Y[2, 2], da_D = Y[2, 1])
end

# ── the wall of marginal stability ─────────────────────────────────────────────────

# Im(a_D/a) at u - the wall is its zero set
function _wall_ratio_im(sw::SeibergWittenSU2, u, F)
    imag(_magnetic_aD_closed(sw, u, F) / _electric_a_closed(sw, u, F))
end

# wall radius along the ray arg(u) = φ (φ strictly off the real axis): grid scan for a
# sign change of Im(a_D/a) in r/2Λ² ∈ (0, 1.5], then bisection
function _wall_radius(sw::SeibergWittenSU2, φ, F)
    R = 2 * F(sw.Λ)^2
    dir = cis(F(φ))
    f(r) = _wall_ratio_im(sw, r * R * dir, F)
    rs = range(F(1) / 20, F(3) / 2; length = 60)
    lo = hi = zero(F)
    flo = zero(F)
    found = false
    rprev, fprev = first(rs), f(first(rs))
    for r in Iterators.drop(rs, 1)
        fr = f(r)
        if fprev * fr ≤ 0
            lo, hi, flo = rprev, r, fprev
            found = true
            break
        end
        rprev, fprev = r, fr
    end
    found || throw(PeriodError(
        "no marginal-stability wall crossing found along arg(u) = $φ"))
    for _ in 1:200
        hi - lo < eps(F) * 100 && break
        mid = (lo + hi) / 2
        fm = f(mid)
        if flo * fm ≤ 0
            hi = mid
        else
            lo, flo = mid, fm
        end
    end
    (lo + hi) / 2 * R
end

"""
    ms_wall(sw::SeibergWittenSU2; n::Integer = 64) -> Vector{<:Complex}

The wall of marginal stability ``\\{u : a_D/a ∈ ℝ\\}`` - the closed curve through the
singular points ``±2Λ²`` separating the strong-coupling chamber (2 BPS states) from
the weak-coupling chamber (W-boson + dyon towers). Traced by bisection of
``\\mathrm{Im}(a_D/a)`` along `n` rays in each half-plane; returns the curve as an
ordered list of `2n + 2` points starting at the monopole point `2Λ²`, through the
upper half-plane to the dyon point `−2Λ²`, and back through the lower half-plane
(conjugation symmetry). See also [`sw_chamber`](@ref).
"""
function ms_wall(sw::SeibergWittenSU2; n::Integer = 64)
    n ≥ 2 || throw(Resurgence.InvalidArgument("n must be ≥ 2, got $n"))
    F = typeof(float(sw.Λ))
    upper = [_wall_radius(sw, k * F(π) / (n + 1), F) * cis(k * F(π) / (n + 1)) for k in 1:n]
    monopole = complex(2 * F(sw.Λ)^2)
    vcat(monopole, upper, -monopole, conj.(reverse(upper)))
end

"""
    sw_chamber(sw::SeibergWittenSU2, u::Number) -> Symbol

Which side of the wall of marginal stability the modulus `u` lies on: `:strong`
(inside - BPS spectrum is monopole + dyon) or `:weak` (outside - W-boson + dyon
towers). Points on the wall itself count as `:weak`; on the real axis the wall meets
`ℝ` only at `±2Λ²`, so real `u` is `:strong` iff `|u| < 2Λ²`. This is the chamber
seam consumed by `su2_bps_states(...; chamber = :auto)`.
"""
function sw_chamber(sw::SeibergWittenSU2, u::Number)
    F = _sw_float(sw, u)
    R = 2 * F(sw.Λ)^2
    if imag(complex(u)) == 0
        return abs(real(u)) < R ? (:strong) : (:weak)
    end
    φ = abs(angle(Complex{F}(u)))
    abs(u) < _wall_radius(sw, φ, F) ? (:strong) : (:weak)
end

# ── monodromy ──────────────────────────────────────────────────────────────────────

"""
    sw_monodromy(sw::SeibergWittenSU2, point::Symbol) -> Matrix{Int}

The ``Sp(2,ℤ)`` monodromy matrix acting on the period vector ``(a_D, a)`` as the modulus
`u` encircles a singular point counterclockwise: `:infinity` (weak coupling),
`:monopole` (`u = 2Λ²`), or `:dyon` (`u = −2Λ²`; lasso through the upper half-plane).
Loops are based at a real weak-coupling point; the matrices are exact integers with
`det = 1`, satisfy ``M_∞ = M_{dyon}·M_{monopole}``, and are **measured** - the
[`continue_periods`](@ref) transport around each loop reproduces them to rounding
(the oracle that pins the `a_D` normalization and every branch choice; see the ledger).
Their invariant charges are the massless states: `(1,0)` for ``M_{monopole}``, `(1,2)`
for ``M_{dyon}``. Throws [`PeriodError`](@ref) on an unknown point.
"""
function sw_monodromy(::SeibergWittenSU2, point::Symbol)
    M_inf = [-1 4; 0 -1]
    M_mon = [1 0; -1 1]
    if point === :infinity
        M_inf
    elseif point === :monopole
        M_mon
    elseif point === :dyon
        # fixed by M_∞ = M_dyon · M_monopole
        Integer.(M_inf * inv(M_mon // 1))
    else
        throw(PeriodError("unknown singular point :$point (expected :infinity, " *
                          ":monopole, or :dyon)"))
    end
end
