# [Takei08], [DLMF-Gamma] 5.11.8. Weber normal form only.
# Principal log(nu); recurrence shifts to Re(nu) >= 8 before summation.
# Oracles: Gamma(1/2), the gamma recurrence and independent gamma values.
using Resurgence: FormalSeries, borel, pade, laplace_sum

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
    weber_problem(E) -> SchrodingerProblem

The Weber normal form ``ħ²ψ″ = (z²/4 − E)ψ``, whose two simple turning points
``±2√E`` merge at `E = 0`. Pass an exact `E` (`Rational`, `Integer`) to keep the WKB
recursion exact.
"""
weber_problem(E) = SchrodingerProblem([zero(E), zero(E), one(E) / 4]; energy = E)

"""
    weber_voros_series(; order = 9) -> FormalSeries{Rational{BigInt}}

The **exact** Weber Voros coefficient as a series in ``x = 1/ν``:

``W = Σ_{n≥1} (2^{1−2n} − 1) B_{2n} / (2n(2n−1)) · x^{2n−1}``

(power offset 1, even coefficients exactly zero, coefficients through `x^order`). This
is the correction part of the Stirling expansion of ``\\log Γ(½ + ν)`` - Takei 2008,
Koike-Takei 2011; DLMF §5.11 - and is retained only for the Weber normal form.

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
``ν = iε`` is the barrier of index `ε`, and
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
