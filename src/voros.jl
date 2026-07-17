# Voros symbols / quantum periods. The Voros symbol of a cycle is ∮ S_odd dz, i.e. the
# integral of the *odd* part of the WKB solution around the cycle - the even periods
# vanish (S_even is a total z-derivative of a single-valued function, so ∮S_even = 0),
# which is exactly why one works with ∮S_odd alone. The result is
#
#   a(ħ) = v₋₁ ħ⁻¹ + v₁ ħ + v₃ ħ³ + …,   v_m = ∮ S_m dz,
#
# an odd-only Laurent series. The leading v₋₁ = ∮√Q is the classical period (the
# future Transseries action) and is stored apart; the ħ¹-and-up tail is a FormalSeries
# with power offset 1 and alternating-zero coefficients - exactly the shape of
# Resurgence's Euler series, so it feeds `borel`/`pade`/`laplace_sum` unchanged.

using Resurgence: FormalSeries, coefficients

"""
    VorosSymbol{F}

The Voros symbol (quantum period) of a cycle: `classical_period` ``v_{-1} = ∮√Q``
stored separately, and `series` the quantum tail ``v₁ħ + v₃ħ³ + …`` as a
`Resurgence.FormalSeries` in `:ħ` with `power_offset = 1` and zero even coefficients.
Build with [`voros_symbol`](@ref).
"""
struct VorosSymbol{F}
    classical_period::Complex{F}
    series::FormalSeries{Complex{F}}
end

"""
    classical_period(vs::VorosSymbol) -> Complex

The classical period ``v_{-1} = ∮√Q`` (the coefficient of ``ħ^{-1}``).
"""
classical_period(vs::VorosSymbol) = vs.classical_period

"""
    quantum_series(vs::VorosSymbol) -> FormalSeries

The quantum tail ``v₁ħ + v₃ħ³ + …`` (power offset 1, even coefficients zero). This is
the object to hand to Resurgence's `borel`/`pade`/`laplace_sum`. Because the odd-only
coefficients alternate with exact zeros, an equal-degree `[n/n]` Padé is degenerate -
pass `order` or `reduce = true` to `pade`.
"""
quantum_series(vs::VorosSymbol) = vs.series

"""
    full_series(vs::VorosSymbol) -> FormalSeries

The complete Voros symbol ``∮S_odd = v_{-1}ħ^{-1} + v₁ħ + v₃ħ³ + …`` as one
`FormalSeries` with `power_offset = -1` (the ``ħ^0`` Maslov slot is zero - that index
is carried by the quantization convention, not the summable series).
"""
function full_series(vs::VorosSymbol)
    FormalSeries(vcat(vs.classical_period, zero(vs.classical_period),
                      coefficients(vs.series)), :ħ; power_offset = -1 // 1)
end

"""
    voros_symbol(w::WKBExpansion, contour; closed = true, atol = nothing,
                 rtol = nothing) -> VorosSymbol

The Voros symbol of the cycle `contour`: computes the quantum periods
`v_m = ∮ S_m dz` for `m = -1, 1, …, order` with the shared √Q branch tracking of
[`wkb_period`](@ref). The even periods `v₂, v₄, …` are asserted `≈ 0` (tolerance
`atol`, default scaled to the period magnitudes) and stored as exact zero. `rtol`
passes through to the quadrature of each period (default `√eps` of the working
float - override to trade precision for speed in energy sweeps).
"""
function voros_symbol(w::WKBExpansion, contour::AbstractVector; closed::Bool = true,
                      atol = nothing, rtol = nothing)
    w.order ≥ 1 ||
        throw(Resurgence.InvalidArgument("voros_symbol needs a WKBExpansion of order ≥ 1"))
    classical = wkb_period(w, contour, -1; closed, rtol)
    F = typeof(real(classical))
    rt = rtol === nothing ? sqrt(eps(F)) : F(rtol)
    v = Vector{Complex{F}}(undef, w.order)
    # Odd periods first, tracking their scale. A pure relative quadrature target is
    # pathological on a cancelling integral (it exhausts maxevals chasing digits of
    # zero - as the odd periods do whenever the series truncates, e.g. the harmonic
    # oscillator), so each also gets an absolute floor at the classical-period scale.
    oddscale = abs(classical)
    for m in 1:2:w.order
        v[m] = wkb_period(w, contour, m; closed, rtol = rt, atol = rt * (1 + abs(classical)))
        oddscale = max(oddscale, abs(v[m]))
    end
    # Even periods vanish identically; they are only sanity-checked (below) at ~1e-3
    # relative to the period scale. Their integrand magnitude is that of the odd
    # periods - up to O(1e6) at high order - so their atol must track `oddscale`, not
    # the (often tiny) classical period: an atol tied to |v₋₁| alone makes the
    # quadrature burn `maxevals` resolving zero to a hopeless precision (the barrier
    # cycle of a shallow double well, 25 s → 0.1 s once this floor is right).
    evtol = F(1e-3) * (1 + oddscale)
    for m in 2:2:w.order
        # Capped eval budget: the high-order 1/Qᵏ integrand has enormous peaks near
        # the turning points, and adaptive quadrature keeps refining them long past
        # the point where the (cancelling, ≈ 0) integral is resolved to `evtol`. A
        # tight cap returns the good rough estimate without the runaway - enough for
        # the vanishing sanity check, never used as a value.
        v[m] = wkb_period(w, contour, m; closed, rtol = F(1e-2), atol = evtol,
                          maxevals = 10^4)
    end
    # The even periods vanish identically (S_even is a total z-derivative of a
    # single-valued function). We store them as exact zero and assert only that the
    # computed value is small relative to the period scale - a generous sanity check
    # (a branch error yields an O(|v₋₁|) even period), not a precision test, since the
    # high-order 1/Qᵏ integrands are numerically stiff near the turning points.
    tol = atol === nothing ?
          F(1e-3) * (1 + abs(classical) + maximum(abs, v; init = zero(F))) : F(atol)
    for m in 2:2:w.order
        abs(v[m]) ≤ tol || throw(Resurgence.InvalidArgument(
            "even quantum period v_$m ≈ $(v[m]) is not ≈ 0 (|·| > $tol); check the " *
            "contour / branch"))
        v[m] = 0
    end
    VorosSymbol{F}(classical, FormalSeries(v, :ħ; power_offset = 1 // 1))
end
