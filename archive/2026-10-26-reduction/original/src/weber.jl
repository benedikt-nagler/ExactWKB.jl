# The Weber (parabolic-cylinder) local model of a merging pair of turning points, and
# its Voros coefficient.
#
# Two simple turning points of `Q` that collide as a parameter is tuned form a *merging
# pair*; at the collision `Q` has a double zero and the Stokes graph is degenerate (four
# rays crossing, see the header of `stokes_graph.jl`). The universal local model is the
# Weber equation `ħ²ψ″ = (z²/4 − E)ψ`, whose turning points are `±2√E`.
#
# The object that makes this layer worth having is the **Voros coefficient**: the
# regularized `S_odd` integral from a turning point out to the irregular singular point
# at infinity, with the classical `ħ^{-1}` term removed. For the Weber equation it is,
# order by order, the correction part of the Stirling expansion of `log Γ(½ + ν)`:
#
#   W(ν) = Σ_{n≥1} (2^{1−2n} − 1) B_{2n} / (2n(2n−1)) · ν^{1−2n}
#   log Γ(½ + ν) = ν log ν − ν + ½ log(2π) + (Borel sum of W)
#
# ([Takei08], [KoikeTakei11]; the `(2^{1−2n} − 1)` factor is what distinguishes
# `log Γ(½+x)` from the plain Stirling series for `log Γ(x)`, [DLMF-Gamma] §5.11.)
# So the package's own WKB recursion, integrated numerically, reproduces a closed-form
# special-function expansion - and Borel-summing it computes Γ itself.
#
# ── Weber convention ledger (each item pinned in test_weber.jl) ─────────────────────
#
# 1. Weber index. `ν = Z / (2πi ħ)` with `Z = ∮√Q` the classical period of the merging
#    cycle, oriented so that `ν > 0` for the normal form at real `E > 0` (then `ν =
#    E/ħ`, and Bohr-Sommerfeld puts the spectrum at `ν = n + ½` - the parabolic-cylinder
#    index). The sign is NOT free: `W` has odd powers of `1/ν`, so flipping it flips
#    every coefficient. PINNED by the `weber_voros_coefficient` ↔ `weber_voros_series`
#    comparison, which fails on the whole vector under `ν → −ν`.
# 2. Factor 2. Our contour computes the integral out to infinity on ONE side, which is
#    exactly HALF of `W`. `weber_voros_coefficient` returns the doubled value, so it is
#    `W` itself with no further convention. (The natural reading: the Weber equation has
#    two turning points and two relevant directions at infinity; one side is half the
#    symmetric object.) PINNED by the same comparison - a missing factor 2 shows up
#    uniformly at every order, the tell-tale of a normalization rather than an error.
# 3. Contour orientation. The path runs `R → z₀ → (loop) → z₀ → R`, inward on one sheet
#    and back out on the other; since `S_odd` is odd across the cut this is `−2∫_{z₀}^R`,
#    hence the overall `−1` in `_voros_path_value`. The turning-point end needs no
#    subtraction: on the double cover `S_m dz` has only EVEN powers of the local
#    coordinate `t` (`z − z₀ = t²`), so the small loop carries no residue and the value
#    is independent of the loop radius. PINNED by the loop-radius-independence test.
# 4. Infinity end. Only `R → ∞` is a genuine limit. The tail is a power series in `1/z`
#    at the far endpoint (leading `(3/8)z^{-2}` at order 1 for the normal form); measured
#    in the path parameter `R` it acquires ODD powers too, because the endpoint is
#    `z₀ + R·dir` rather than `R`. We therefore do not assume an even tail: the value is
#    Neville-extrapolated to `1/R = 0` through four radii, killing every power through
#    `R^{-3}`. This, not the quadrature, sets the accuracy floor - the residual is
#    `rtol`-independent and scales like a power of `R`. PINNED by the coefficient
#    comparison (a single radius is off by ~2e-3 at order 1; an even-only two-step
#    Richardson stalls at ~2e-7, which is the odd term this replaces).
# 5. Which side. The two turning points of a pair give coefficients that differ by an
#    overall SIGN, not by nothing: `which` selects a genuinely different *oriented* path
#    on the double cover, and reversing the outward direction reverses `∫…dz` while the
#    principal branch of `√Q` at the far endpoint does not follow. This is a property of
#    the object, not a defect - a Voros coefficient belongs to an oriented path - so it
#    is recorded rather than normalized away. `which = 2` (the turning point listed
#    second, the one further along `+dir`) is the branch that reproduces
#    `weber_voros_series` with a plus sign. PINNED by the two-sided test.
# 6. Which branch of `log ν` in the elementary part. `weber_voros_sum` is defined as
#    `log Γ(½+ν) − (ν log ν − ν + ½log 2π)` with the PRINCIPAL log, and `ν = 0` is the
#    removable case `elementary(0) = ½log 2π`. This is not cosmetic: on the imaginary
#    axis `ν = iε` the principal branch flips at `ε = 0`, and that flip is exactly what
#    makes `Re` of the sum equal `−½log(1+e^{−2π|ε|})` on BOTH sides - the identity the
#    uniform quantization condition of `quantization.jl` is built on. PINNED by the
#    closed-form real-part test at ε of both signs.
# 7. The imaginary axis needs the recurrence, not a different ray. The Borel
#    singularities of the Stirling series sit at `2πik`, i.e. exactly along the Laplace
#    ray `arg(1/ν) = ∓π/2` that `ν = iε` asks for. Rather than choose a lateral, we
#    shift with `log Γ(½+ν) = log Γ(½+ν+k) − Σ_{j<k} log(½+ν+j)` until the ray is well
#    inside a Borel-regular sector (`Re ν ≥ 8` and `|Im ν| ≤ Re ν`, so `|arg x| ≤ π/4`).
#    The recurrence is exact, so this costs accuracy nowhere. PINNED by the functional
#    equation and by `Γ(½) = √π` at `ν = 0`, which only the shift can reach.

using Resurgence: FormalSeries, borel, pade, laplace_sum

# -- Bernoulli numbers (in-file, exact) ----------------------------------------------
#
# The standard recurrence Σ_{k=0}^{m} binom(m+1, k) B_k = 0, in the B₁ = −½ convention
# (irrelevant here - only the even-index values are used). Kept in-file rather than
# taking a dependency, following the in-file Carlson K/E precedent in `sw_curve.jl`.
function _bernoulli(n::Integer)
    B = Vector{Rational{BigInt}}(undef, n + 1)
    B[1] = one(Rational{BigInt})
    for m in 1:n
        s = zero(Rational{BigInt})
        for k in 0:(m - 1)
            s += binomial(BigInt(m + 1), BigInt(k)) * B[k + 1]
        end
        B[m + 1] = -s // (m + 1)
    end
    B
end

# -- the model -----------------------------------------------------------------------

"""
    WeberModel{F}

The Weber (parabolic-cylinder) local model of a merging pair of turning points: their
locations, the classical period `merging_period` ``Z = ∮√Q`` of the cycle around them,
and the Weber index scale ``Z/(2πi)`` (so ``ν = `` [`weber_index`](@ref)``(m, ħ)``).
Build with [`weber_model`](@ref); the reference problem is [`weber_problem`](@ref).
"""
struct WeberModel{F}
    turning_points::Tuple{Complex{F},Complex{F}}
    merging_period::Complex{F}
    index_scale::Complex{F}
end

"""
    merging_period(m::WeberModel) -> Complex

The classical period ``Z = ∮√Q`` of the cycle around the merging pair.
"""
merging_period(m::WeberModel) = m.merging_period

"""
    weber_index(m::WeberModel, ħ) -> Complex

The Weber index ``ν = Z/(2πi ħ)`` (ledger item 1). For the normal form
[`weber_problem`](@ref)`(E)` this is ``E/ħ``, and Bohr-Sommerfeld puts the spectrum at
``ν = n + ½`` - the parabolic-cylinder index of ``Γ(½ + ν)``.
"""
weber_index(m::WeberModel, ħ) = m.index_scale / ħ

"""
    weber_problem(E) -> SchrodingerProblem

The Weber normal form ``ħ²ψ″ = (z²/4 − E)ψ``, whose two simple turning points
``±2√E`` merge at `E = 0`. Pass an exact `E` (`Rational`, `Integer`) to keep the WKB
recursion exact.
"""
weber_problem(E) = SchrodingerProblem([zero(E), zero(E), one(E) / 4]; energy = E)

"""
    weber_model(prob::SchrodingerProblem; pair = nothing, margin = nothing, n = 64,
                rtol = nothing) -> WeberModel

The Weber model of a merging pair of `prob`. `pair` selects the two turning points by
index into [`turning_points`](@ref); by default the closest pair of simple turning
points is used. The classical period is computed with [`period_integral`](@ref) over an
[`encircling_contour`](@ref), and the sign of ``ν`` is fixed by ledger item 1
(`ν > 0` for the normal form at real `E > 0`).

Throws [`CoalescentTurningPoints`](@ref) if the selected points have already merged
(the model is the *approach* to the degenerate graph, not the degenerate graph itself -
at the merge point `Z = 0` and `ν` is not defined).
"""
function weber_model(prob::SchrodingerProblem; pair = nothing, margin = nothing,
                     n::Integer = 64, rtol = nothing)
    tps = turning_points(prob)
    F = _wkb_float(prob)
    simple = [t for t in tps if is_simple(t)]
    length(simple) ≥ 2 || throw(CoalescentTurningPoints(
        energy(prob), isempty(tps) ? zero(Complex{F}) : location(first(tps)), zero(F)))

    if pair === nothing
        best, bi, bj = F(Inf), 0, 0
        for i in eachindex(simple), j in (i + 1):length(simple)
            d = abs(location(simple[i]) - location(simple[j]))
            d < best && ((best, bi, bj) = (d, i, j))
        end
        z1, z2 = location(simple[bi]), location(simple[bj])
    else
        i, j = pair
        z1, z2 = location(tps[i]), location(tps[j])
        (is_simple(tps[i]) && is_simple(tps[j])) || throw(CoalescentTurningPoints(
            energy(prob), (z1 + z2) / 2, abs(z2 - z1)))
    end
    sep = abs(z2 - z1)
    sep ≤ eps(F)^(3 // 4) * (1 + abs(z1)) &&
        throw(CoalescentTurningPoints(energy(prob), (z1 + z2) / 2, sep))

    ct = encircling_contour(z1, z2; margin, n)
    Z = period_integral(prob, ct; rtol)
    ν = Z / (2im * F(π))
    # Ledger item 1: orient so ν > 0 on the normal form's real-E branch. The cycle's
    # orientation is a free choice; the *series* is not (odd powers of 1/ν).
    if real(ν) < 0 || (iszero(real(ν)) && imag(ν) < 0)
        Z = -Z
        ν = -ν
    end
    WeberModel{F}((Complex{F}(z1), Complex{F}(z2)), Complex{F}(Z), Complex{F}(ν))
end

# -- the exact series ----------------------------------------------------------------

"""
    weber_voros_series(; order = 9) -> FormalSeries{Rational{BigInt}}

The **exact** Weber Voros coefficient as a series in ``x = 1/ν``:

``W = Σ_{n≥1} (2^{1−2n} − 1) B_{2n} / (2n(2n−1)) · x^{2n−1}``

(power offset 1, even coefficients exactly zero, coefficients through `x^order`). This
is the correction part of the Stirling expansion of ``\\log Γ(½ + ν)`` - Takei 2008,
Koike-Takei 2011; DLMF §5.11 - and therefore a closed form that the numerically
computed [`weber_voros_coefficient`](@ref) must reproduce.

The series is Gevrey-1 divergent, with Borel singularities at ``2πi k``: hand it to
`Resurgence.borel` / `pade` / `laplace_sum`, or use [`weber_log_gamma`](@ref).
"""
function weber_voros_series(; order::Integer = 9)
    order ≥ 1 || throw(Resurgence.InvalidArgument("order must be ≥ 1, got $order"))
    B = _bernoulli(order + 1)
    c = zeros(Rational{BigInt}, order)
    for n in 1:((order + 1) ÷ 2)
        m = 2n - 1
        m ≤ order || break
        c[m] = (1 // (BigInt(2)^(2n - 1)) - 1) * B[2n + 1] // (2n * (2n - 1))
    end
    FormalSeries(c, :x; power_offset = 1 // 1)
end

"""
    weber_voros_series(m::WeberModel; order = 9) -> FormalSeries

The same series re-expanded in `ħ` for a concrete model: ``x = 1/ν = 2πi ħ/Z``, so the
coefficient of ``ħ^{2n−1}`` is ``c_n (2πi/Z)^{2n−1}``. This is the prediction that
[`weber_voros_coefficient`](@ref) is compared against.
"""
function weber_voros_series(m::WeberModel{F}; order::Integer = 9) where {F}
    S = weber_voros_series(; order)
    u = 1 / m.index_scale                       # x = u·ħ
    c = Complex{F}[Complex{F}(a) * u^k for (k, a) in enumerate(Resurgence.coefficients(S))]
    FormalSeries(c, :ħ; power_offset = 1 // 1)
end

# -- the computed coefficient --------------------------------------------------------

# The contour of ledger item 3: out at radius R on the ray through z₀, in to a small
# loop about z₀, around it once (so √Q flips sheet), and back out to R. Intermediate
# vertices keep the adaptive quadrature (and the branch tracker's grazing test) happy on
# the long outward leg.
function _voros_path(z₀::Complex{F}, dir::Complex{F}, R::F, δ::F, n::Integer) where {F}
    far = z₀ + R * dir
    mid = [z₀ + s * dir for s in (10 * δ, 3 * δ) if s < R]
    v = Complex{F}[far]
    append!(v, mid)
    push!(v, z₀ + δ * dir)
    for k in 1:(n - 1)
        push!(v, z₀ + δ * dir * cis(2 * F(π) * k / n))
    end
    push!(v, z₀ + δ * dir)
    append!(v, reverse(mid))
    push!(v, far)
    v
end

# −∮ over that path = 2∫_{z₀}^{R} S_m dz (ledger items 2 and 3).
function _voros_path_value(w::WKBExpansion, path, m::Integer; rtol)
    -wkb_period(w, path, m; closed = false, rtol)
end

# Neville extrapolation of `vals` sampled at nodes `hs` to h = 0 (ledger item 4). No
# assumption on which powers of h appear - it fits the interpolating polynomial.
function _neville_at_zero(hs::AbstractVector, vals::AbstractVector)
    t = collect(vals)
    for k in 1:(length(t) - 1), i in 1:(length(t) - k)
        t[i] = ((0 - hs[i + k]) * t[i] - (0 - hs[i]) * t[i + 1]) / (hs[i] - hs[i + k])
    end
    t[1]
end

"""
    weber_voros_coefficient(prob::SchrodingerProblem, model::WeberModel; order = 9,
                            which = 2, radius = nothing, loop_radius = nothing,
                            n = 24, rtol = nothing) -> FormalSeries

The Voros coefficient of `model`'s merging pair, **computed** from `prob`'s own WKB
recursion: the regularized ``S_odd`` integral from the turning point out to the
irregular singular point at infinity, doubled (ledger item 2), returned as a
`Resurgence.FormalSeries` in `:ħ` with power offset 1 (even coefficients exactly zero,
like every Voros series here).

`which ∈ (1, 2)` selects the turning point of the pair to start from, and the outward
direction is the one pointing away from its partner. The two sides give coefficients
differing by an overall **sign** - a Voros coefficient belongs to an *oriented* path,
and reversing the outward direction reverses it (ledger item 5); `which = 2` (the
default) is the branch that matches [`weber_voros_series`](@ref). `radius` is the outer
cut-off `R`
(the value is Richardson-extrapolated from `R` and `2R` in ``1/R²``, ledger item 4),
`loop_radius` the small-loop radius `δ` about the turning point (the value is
independent of it - ledger item 3 - but the quadrature is not: the default, a quarter
of the pair separation, keeps the high-order integrands away from the turning point),
and `n` the loop resolution.

For the Weber normal form [`weber_problem`](@ref) this reproduces
[`weber_voros_series`](@ref)`(model)` - the headline oracle of this layer.
"""
function weber_voros_coefficient(prob::SchrodingerProblem, model::WeberModel{F};
                                 order::Integer = 9, which::Integer = 2,
                                 radius = nothing, loop_radius = nothing,
                                 n::Integer = 24, rtol = nothing) where {F}
    order ≥ 1 || throw(Resurgence.InvalidArgument("order must be ≥ 1, got $order"))
    which in (1, 2) || throw(Resurgence.InvalidArgument("which must be 1 or 2, got $which"))
    z₀ = model.turning_points[which]
    other = model.turning_points[3 - which]
    sep = abs(z₀ - other)
    dir = (z₀ - other) / sep                    # outward: away from the partner
    δ = F(something(loop_radius, sep / 4))
    R = F(something(radius, 50 * (1 + abs(z₀))))
    # The R-tail, not the quadrature, is the accuracy floor (ledger item 4), so the
    # default tolerance is capped rather than following `eps(F)` down into BigFloat.
    rt = rtol === nothing ? max(eps(F)^(3 // 4), F(1e-13)) : F(rtol)

    w = wkb_expansion(prob; order)
    radii = F[s * R for s in (1, 2, 3, 4)]
    paths = [_voros_path(z₀, dir, r, δ, n) for r in radii]
    hs = F[1 / r for r in radii]
    v = zeros(Complex{F}, order)
    for m in 1:2:order
        vals = [_voros_path_value(w, p, m; rtol = rt) for p in paths]
        v[m] = _neville_at_zero(hs, vals)       # extrapolate R → ∞ (ledger item 4)
    end
    FormalSeries(v, :ħ; power_offset = 1 // 1)
end

# -- the connection constant ----------------------------------------------------------

# The elementary (non-resurgent) part of `log Γ(½+ν)`: `ν log ν − ν + ½log 2π`, with
# the principal log and the removable value at ν = 0 (ledger item 6).
function _weber_elementary(ν::Number)
    T = Complex{typeof(float(real(ν)))}
    iszero(ν) && return T(log(2 * T(π)) / 2)
    T(ν * log(T(ν)) - ν + log(2 * T(π)) / 2)
end

# How far to walk the recurrence so `x = 1/ν` lands in a Borel-regular sector
# (ledger item 7): `Re ν ≥ 8` and `|Im ν| ≤ Re ν`.
function _weber_shift(ν::Number)
    re, im = float(real(ν)), abs(float(imag(ν)))
    max(0, ceil(Int, max(8 - re, im - re)))
end

"""
    weber_voros_sum(ν; order = 9, shift = nothing, pade_order = nothing,
                    rtol = nothing) -> Complex

The **Borel sum** of the Weber Voros coefficient at `ν` - i.e.
``\\log Γ(½ + ν) − (ν \\log ν − ν + ½\\log 2π)``, the resurgent remainder that
[`weber_voros_series`](@ref) is the asymptotic expansion of. Computed by
Borel-Padé-Laplace on that series at ``x = 1/ν``; no external Gamma function is
involved.

`ν` may sit anywhere off the negative real axis, **including on the imaginary axis**,
where the Borel singularities at ``2πi k`` lie exactly on the Laplace ray: the exact
recurrence ``\\log Γ(½+ν) = \\log Γ(½+ν+k) − Σ_{j<k} \\log(½+ν+j)`` walks `ν` into a
regular sector first (`shift` overrides the automatic `k`). A purely imaginary
``ν = iε`` is the barrier of index `ε` (see [`weber_barrier_amplitude`](@ref)), and
there the real part is the elementary ``−½\\log(1 + e^{−2π|ε|})``.
"""
function weber_voros_sum(ν::Number; order::Integer = 9, shift = nothing,
                         pade_order = nothing, rtol = nothing)
    T = Complex{typeof(float(real(ν)))}
    k = shift === nothing ? _weber_shift(ν) : Int(shift)
    k ≥ 0 || throw(Resurgence.InvalidArgument("shift must be ≥ 0, got $k"))
    νs = T(ν) + k
    x = 1 / νs
    S = weber_voros_series(; order)
    Sf = FormalSeries(T[T(c) for c in Resurgence.coefficients(S)], :x;
                      power_offset = 1 // 1)
    B = borel(Sf)
    r = pade(B; order = pade_order, reduce = true)
    G = laplace_sum(B, r, x; θ = angle(x), rtol)
    # undo the shift: log Γ(½+ν) = log Γ(½+ν+k) − Σ log(½+ν+j), then strip the
    # elementary part at the ORIGINAL ν (ledger items 6, 7)
    for j in 0:(k - 1)
        G -= log(T(ν) + T(1 // 2) + j)
    end
    T(G + _weber_elementary(νs) - _weber_elementary(ν))
end

"""
    weber_log_gamma(ν; order = 9, shift = nothing, pade_order = nothing,
                    rtol = nothing) -> Complex

``\\log Γ(½ + ν)``, computed **by resurgence**: the elementary part
``ν \\log ν − ν + ½\\log(2π)`` plus [`weber_voros_sum`](@ref), the Borel-Padé-Laplace
sum of [`weber_voros_series`](@ref) at ``x = 1/ν``. No external Gamma function is
involved - this is the Weber Voros coefficient doing the work.

Valid for `ν` away from the negative real axis (where ``Γ`` has its poles); small and
imaginary `ν` are reached through the exact recurrence, so `ν = 0` gives
``\\log Γ(½) = ½\\log π``. `pade_order` truncates the Padé approximant
(`reduce = true` is mandatory and applied - the series is odd-only).
"""
weber_log_gamma(ν::Number; kwargs...) =
    _weber_elementary(ν) + weber_voros_sum(ν; kwargs...)

"""
    weber_connection(ν; kwargs...) -> Complex

The Weber connection constant ``\\sqrt{2π}/Γ(½ + ν)`` - the factor that a degenerate
(double) turning point contributes to a connection formula, where a simple turning
point contributes the Airy constant. Built from [`weber_log_gamma`](@ref), so it too
computes ``Γ`` by Borel summation of the Voros coefficient. `kwargs` pass through.
"""
weber_connection(ν::Number; kwargs...) =
    sqrt(2 * oftype(float(real(ν)), π)) * exp(-weber_log_gamma(ν; kwargs...))

"""
    weber_barrier_amplitude(ε; kwargs...) -> Real

The transmission amplitude ``|T|`` of a parabolic barrier of index `ε`, computed **by
resurgence** from [`weber_connection`](@ref) at ``ν = iε``:

``|T| = e^{-πε/2}\\,/\\,|\\sqrt{2π}/Γ(½+iε)| = (1 + e^{2πε})^{-1/2}``,

the closed form on the right being the oracle. `ε > 0` is a barrier the particle
tunnels through (``|T| ≈ e^{-πε}``, the one-instanton weight); `ε = 0` is the barrier
top (``|T| = 1/\\sqrt2``); `ε < 0` is over-barrier (``|T| → 1``). This is the
connection constant a *merging pair* of turning points contributes where a single
simple turning point contributes the Airy ``½`` - the content of the uniform
quantization condition of [`quantization_condition`](@ref).
"""
function weber_barrier_amplitude(ε::Real; kwargs...)
    F = typeof(float(ε))
    exp(-F(π) * ε / 2) / abs(weber_connection(complex(zero(F), F(ε)); kwargs...))
end
