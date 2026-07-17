# The ZJJ / Dunne–Ünsal layer for the symmetric double well: the two quantum periods
# as ħ-series with E-dependent coefficients, the exact perturbative/nonperturbative
# (P/NP) relation between them, and the median-summed energy splitting.
#
# ── ZJJ/DU convention ledger (pinned in test_double_well.jl) ────────────────────────
#
# 1. The two characteristic functions, in Zinn-Justin–Jentschura's normalization
#    (Jentschura–Zinn-Justin, "Multi-instantons and exact results I", eqs. 2.23–2.25):
#      B(E, ħ) = −(1/2πi)·∮_well S_odd   (perturbative condition B = n + ½),
#      A(E, ħ) = −∮_barrier S_odd        (splitting weight e^{−A/2}; A ≈ 2S_I/ħ > 0),
#    both odd Laurent series in ħ with offset −1 (`perturbative_b`, `instanton_a`).
#    The harmonic pin: for the ZJJ potential V = q²(1−q)²/2 at ħ = g/√2 these match
#    ZJJ's B(E_zjj, g) and A(E_zjj, g) under E_ours = g·E_zjj.
# 2. Dunne–Ünsal relation (DU 2014; Gahramanov–Tezgin 1512.08466 eq. 4): for the
#    ZJJ functions at fixed B,  ∂E/∂B = −(g/2S)·(2B + g ∂A_zjj/∂g),  S = the
#    instanton action coefficient; genus-one potentials only (NOT quintic+).
#    **Period language**: ZJJ's A is not the barrier quantum action - equating the
#    Γ-form condition with the parity condition (item 3 of the quantization ledger)
#    through the reflection formula gives the Weber dictionary
#      A_zjj = A_per + 2B·ln(2/g) − 2·lnΓ(½+B) + const,
#    and `g ∂/∂g|_B` kills the Γ and constant terms while the log contributes −2B:
#    the 2B in the DU bracket cancels exactly. For the quantum periods themselves,
#      ∂E/∂B|_ħ = −c·ħ²·[ħ·(∂A/∂ħ − (∂A/∂E)(∂B/∂ħ)/(∂B/∂E))]   (all at fixed E),
#    one potential-specific constant c > 0. Exact values from mapping coordinates
#    onto ZJJ's potential: c = 6 for V = q²(1−q)²/2, c = 3/2 for V = (z²−1)²
#    (z = 2q−1). `verify_zjj` fits c at leading ħ-order and verifies every higher
#    order - each order is a theorem check ("every theorem is a test").
# 3. Splitting. The parity eigenvalues from `wkb_eigenvalue(…; parity = ±1)` obey
#    the ZJJ one-instanton law: for the ZJJ potential at N = 0,
#      E₋ − E₊ = g·[2ξ(g)·(1 − (71/12)g − (6299/288)g² + O(g³))],
#      ξ(g) = e^{−1/6g}/√(πg)   (ZJJ-I eqs. 2.19–2.20, coefficients e_{0,10l}),
#    pinned together with the exact g = 0.08 eigenvalues quoted in ZJJ-I fig. 2.1.
# 4. Two-instanton log sector (the M6b consumer contract, tested in
#    test_double_well.jl): expanding the exact ZJJ quantization condition
#      1/Γ(½−B) + (εi/√2π)·(−2/g)^B e^{−A/2} = 0    (ZJJ-I eq. 2.22)
#    to second order in ξ produces exactly one log: E^{(2)} = ξ²·(ln(−2/g) + γ),
#    i.e. e_{0,210} = 1, e_{0,200} = γ (ZJJ-I eq. 2.19c). The log is minted by the
#    resonant (1,1) channel of the A/−A lattice - `Resurgence.resonant_solve` with
#    an explicitly passed lattice, since the numeric actions are inexact reals.

# ── the two quantum periods as ħ-series ─────────────────────────────────────────

# the validated (well, barrier) cycles of the double-well layout at energy E
function _dw_cycles(prob::SchrodingerProblem, E)
    cycles = spectral_cycles(prob, E)
    wells = _wells(cycles)
    barriers = _barriers(cycles)
    (length(wells) == 2 && length(barriers) == 1) || throw(QuantizationError(
        "the double-well layer needs a symmetric double-well layout (two wells, one " *
        "barrier) at E = $E, found $(length(wells)) well(s), $(length(barriers)) " *
        "barrier(s)"))
    wells[1], barriers[1]
end

# B from the well symbol / A from the barrier symbol (ledger item 1)
function _b_from(vs::VorosSymbol)
    T = typeof(classical_period(vs))
    (T(im) / (2 * oftype(real(one(T)), π))) * full_series(vs)
end
_a_from(vs::VorosSymbol) = -1 * full_series(vs)

# both series from ONE WKB expansion (the shared stencil workhorse of verify_zjj)
function _dw_ab(prob::SchrodingerProblem, E; order, quad_rtol = nothing)
    well, barrier = _dw_cycles(prob, E)
    w = wkb_expansion(with_energy(prob, E); order)
    (_b_from(voros_symbol(w, well.contour; rtol = quad_rtol)),
     _a_from(voros_symbol(w, barrier.contour; rtol = quad_rtol)))
end

"""
    perturbative_b(prob::SchrodingerProblem, E; order = 12, quad_rtol = nothing)
        -> FormalSeries

The ZJJ perturbative function ``B(E, ħ) = −(1/2πi)\\,∮_{well} S_{odd}`` of a
symmetric double well as an odd ħ-Laurent series (`power_offset = -1`), normalized
so the perturbative quantization condition is ``B = n + ½`` (ledger item 1). The
coefficients are E-dependent complex numbers with vanishing imaginary part.
"""
function perturbative_b(prob::SchrodingerProblem, E; order::Integer = 12,
                        quad_rtol = nothing)
    well, _ = _dw_cycles(prob, E)
    w = wkb_expansion(with_energy(prob, E); order)
    _b_from(voros_symbol(w, well.contour; rtol = quad_rtol))
end

"""
    instanton_a(prob::SchrodingerProblem, E; order = 12, quad_rtol = nothing)
        -> FormalSeries

The ZJJ instanton function ``A(E, ħ) = −∮_{barrier} S_{odd}`` of a symmetric double
well as an odd ħ-Laurent series (`power_offset = -1`), normalized so the splitting
weight is ``e^{−A/2}`` with ``A = 2S_I/ħ + O(ħ)`` positive (ledger item 1).
"""
function instanton_a(prob::SchrodingerProblem, E; order::Integer = 12,
                     quad_rtol = nothing)
    _, barrier = _dw_cycles(prob, E)
    w = wkb_expansion(with_energy(prob, E); order)
    _a_from(voros_symbol(w, barrier.contour; rtol = quad_rtol))
end

# ── series helpers (local; FormalSeries has no division) ────────────────────────

# multiply by ħ^s: shift the power offset, keep every coefficient (a monomial
# Cauchy product would truncate to one term - Resurgence `*` keeps min(n_terms))
_dw_shift(Φ::Resurgence.FormalSeries, s::Integer) =
    Resurgence.FormalSeries(Resurgence.coefficients(Φ), Φ.var;
                            power_offset = Φ.power_offset + s)

# multiplicative inverse of Φ (leading coefficient nonzero), to `n` coefficients
function _series_inv(Φ::Resurgence.FormalSeries{T}, n::Integer) where {T}
    c = Resurgence.coefficients(Φ)
    iszero(c[1]) && throw(Resurgence.InvalidArgument(
        "cannot invert a series with vanishing leading coefficient"))
    d = Vector{T}(undef, n)
    d[1] = inv(c[1])
    for k in 2:n
        acc = zero(T)
        for j in 1:(k - 1)
            acc += d[j] * (k - j + 1 ≤ length(c) ? c[k - j + 1] : zero(T))
        end
        d[k] = -acc / c[1]
    end
    Resurgence.FormalSeries(d, Φ.var; power_offset = -Φ.power_offset)
end

"""
    verify_zjj(prob::SchrodingerProblem, E; order = 13, delta = nothing, c = nothing,
               quad_rtol = nothing) -> NamedTuple

Verify the Dunne–Ünsal P/NP relation (ledger item 2) for the symmetric double well
`prob` **order by order in ħ as an identity at energy `E`** - no small-ħ limit, no
well-bottom coalescence. Computes ``B(E,ħ)``/``A(E,ħ)`` on a 5-point E-stencil
(spacing `delta`), forms both sides of the period-language relation

``1/(∂B/∂E) = −c\\,ħ³\\,(∂A/∂ħ|_B)``

as ħ-series, fits the constant `c` from the leading order (or checks a supplied
`c`), and returns `(c, orders, residuals)` - the relative residual at each higher
odd order, every one of which is an independent theorem check. Use exact/BigFloat
potentials with an explicit `quad_rtol` for clean residuals.
"""
function verify_zjj(prob::SchrodingerProblem, E; order::Integer = 13,
                    delta = nothing, c = nothing, quad_rtol = nothing)
    F = promote_type(_wkb_float(prob), typeof(float(real(E))))
    EF = F(real(E))
    δ = delta === nothing ?
        (F === Float64 ? F(1 // 500) : F(1 // 10)^6) * (1 + abs(EF)) : F(delta)
    # B and A on the 5-point stencil (one shared WKB expansion per stencil energy)
    Bs = Vector{Any}(undef, 5)
    As = Vector{Any}(undef, 5)
    for (i, k) in enumerate(-2:2)
        Bs[i], As[i] = _dw_ab(prob, EF + k * δ; order, quad_rtol)
    end
    B0, A0 = Bs[3], As[3]
    fd(f) = (f[1] - 8 * f[2] + 8 * f[4] - f[5]) * (1 / (12 * δ))   # 4th-order ∂/∂E
    dBdE = fd(Bs)
    dAdE = fd(As)
    dBdh = Resurgence.derivative(B0)
    dAdh = Resurgence.derivative(A0)
    T = eltype(typeof(B0))
    L = Resurgence.n_terms(B0)
    inv_dBdE = _series_inv(dBdE, L)
    # ∂A/∂ħ at fixed B = ∂A/∂ħ|_E − (∂A/∂E)(∂B/∂ħ|_E)/(∂B/∂E)
    dAdh_B = dAdh - dAdE * dBdh * inv_dBdE
    lhs = inv_dBdE                                    # ∂E/∂B, power offset +1
    rhs = _dw_shift(dAdh_B, 3)                        # ħ³·(∂A/∂ħ|_B)
    cfit = c === nothing ? -lhs[0] / rhs[0] : T(c)
    resid = lhs + cfit * rhs
    pmax = order - 2
    orders = collect(3:2:pmax)
    residuals = map(orders) do p
        k = p - 1                                     # 0-based index at offset +1
        scale = max(abs(lhs[k]), abs(cfit * rhs[k]))
        iszero(scale) ? abs(resid[k]) : abs(resid[k]) / scale
    end
    (c = cfit, orders = orders, residuals = residuals)
end

"""
    energy_splitting(prob::SchrodingerProblem, n::Integer, ħ; kwargs...) -> Real

The median-summed even/odd splitting ``ΔE_n = E_{odd,n} − E_{even,n}`` of doublet
`n` of a symmetric double well: the difference of the two parity roots of the exact
quantization condition ([`wkb_eigenvalue`](@ref) at `parity = ∓1`). Positive - the
even (symmetric) member lies below. Keywords pass to [`wkb_eigenvalue`](@ref).
One-instanton law: ``ΔE ≈ (dE/dB)·(2/π)·e^{−A/2}``-type, pinned quantitatively by
the ZJJ oracle in `test_double_well.jl` (ledger item 3).
"""
function energy_splitting(prob::SchrodingerProblem, n::Integer, ħ; kwargs...)
    E_even = wkb_eigenvalue(prob, n, ħ; parity = +1, kwargs...)
    E_odd = wkb_eigenvalue(prob, n, ħ; parity = -1, kwargs...)
    E_odd - E_even
end
