# Exact (Voros) quantization conditions, spectral determinants, and the eigenvalue
# solver. The conditions are built from median-summed Voros multipliers (`_log_voros`,
# the DDP layer's summation - never re-implemented here) of the spectral cycles at
# energy E; the eigenvalues are the Newton roots in E.
#
# ── Quantization ledger (each item pinned by an oracle in test_quantization /
#    test_double_well) ────────────────────────────────────────────────────────────
#
# 1. Well multiplier. With the spectral cycle orientation (quantum_periods ledger) the
#    well cycle's summed exponent is log V_B = −i·J_qu(E, ħ)/ħ with J_qu > 0 the
#    all-orders quantum action. SINGLE-WELL exact quantization (all-orders
#    Bohr–Sommerfeld): 1 + V_B = 0 in the summed sense, i.e.
#        imag(log V_B) + 2π(n + ½) = 0,   n = 0, 1, ….
#    Pinned by the harmonic oracle, where the WKB series truncates and the condition
#    is EXACT at every order: E_n = 2ħ(n + ½) for Q = z² − E to full precision.
# 2. Laplace ray. Real potentials with real spectra put every Borel singularity of a
#    well symbol on the real axis (barrier saddles) or off it (complex saddles), so
#    θ = 0 is the physical ray; `side = :median` is the canonical value on it and is
#    mandatory in the double-well case (the A/−A resonant ray).
# 3. Symmetric double well (parity-factorized DDP condition). With φ(E) =
#    J_qu,med/(2ħ) = −imag(log V_B^med)/2 (half the well action) and V_A^med the
#    median-summed BARRIER multiplier (real, exponentially small; our barrier
#    orientation gives log V_A ≈ −2S_I/ħ < 0, so V_A = e^{−A} in ZJJ's A):
#        cos φ = (σ/2)·√(V_A / (1 + V_A)),   σ = parity · (−1)ⁿ,
#    parity = +1 the even (symmetric) member of doublet n - always the lower one -
#    and −1 the odd. The factor ½ is the connection-formula tunneling amplitude
#    (the Airy 1/2 that a single barrier contributes) - PINNED by the oracle: with
#    it, the diagonalization eigenvalues converge with WKB order (3.7e-11 at
#    ħ = 0.1, order 10) and the splitting ratio → 1; WITHOUT it the mean is still
#    exact but the splitting is exactly 2× too wide (an order-independent factor,
#    the tell-tale of a missing connection constant). V_A → 0 recovers the per-well
#    condition φ = π(n + ½); the linearized splitting Δφ = √V_A = e^{−S_I/ħ}(1+O(ħ))
#    is the one-instanton factor. The all-orders barrier symbol V_A resums what ZJJ
#    writes as the separate (2/g)^B, Γ(½+B), and prefactor pieces (why V_A alone,
#    with no extra Γ, gives the right magnitude). The (1 + V_A)^{−1/2} tail is a
#    two-instanton (e^{−2S_I/ħ}) correction, below what the oracle resolves.
# 4. Spectral determinant. Thin, from the same ingredients: single well
#    D = 1 + V_B; double well D = cos²φ − ¼·V_A/(1+V_A) = ∏_± D_±, the parity
#    factors D_± = cos φ ∓ ½√(V_A/(1+V_A)). Zeros ⟺ eigenvalues; no Hadamard/ζ
#    regularization (deliberately out of the spectral layer's scope).
# 5. UNIFORM (Weber) condition, `uniform = true`. Item 3's ½ is the AIRY connection
#    constant, so it is the deep-barrier form: it has no continuation past the barrier
#    top, where the barrier's two turning points merge and the connection becomes the
#    parabolic-cylinder (Weber) one. Written as a phase condition on the SAME two
#    all-orders quantities - `φ = −Im log V_B/2` on the continued well cycle and the
#    barrier index `ε` of `uniform_cycles` - the uniform condition is
#        2φ + p·arctan(e^{−πε}) = π(2n + 1),   p = parity (+1 even, −1 odd).
#    It comes from the exact parabolic-barrier amplitudes: `|T| = (1+e^{2πε})^{−1/2}`
#    and `|T|/|R| = e^{−πε}` ([`weber_barrier_amplitude`](@ref) computes `|T|` from the
#    package's own resurgent Γ), with the round trip e^{2iφ}·e^{−iπ/2}·r_p = 1 and the
#    parity channel r_p = e^{G}(p e^{−πε} − i), `G` the Weber Voros coefficient. The
#    real part of that complex equation cancels IDENTICALLY on both sides of the top -
#    `Re G = −½log(1 + e^{−2π|ε|})` against the modulus of `(p e^{−πε} − i)`, the
#    branch of `log ν` flipping at ε = 0 (weber ledger item 6) - which is why the
#    residual is real above the top as well as below. Two limits pin it:
#      * ε → +∞: arctan(e^{−πε}) → e^{−πε} = √V_A, i.e. cos φ = σ√V_A/2, item 3 at
#        one instanton. The two differ at O(V_A^{3/2}) - two instantons - where item 3
#        is the exact (DDP) one and this is the parabolic model's.
#      * ε → −∞: arctan → π/2, giving 2φ = π(2n + 1 − p/2), i.e. the single-well
#        Bohr-Sommerfeld spectrum with the even levels at N = 2n and the odd at
#        N = 2n + 1 - the interleaving of a symmetric well with no barrier left.
#    PINNED end to end by dense diagonalization at ħ = 0.25 for the quartic double
#    well, on levels both below the barrier top (where it reproduces item 3) and above
#    it (where item 3 does not exist).
#    NOT uniform in the last window: `φ` is a period of the cycle PINCHED by the
#    vanishing one, so its ħ-series diverges near the top and no truncation order
#    rescues it (measured for V = (z²−1)² at ħ = 0.25: fine at |E − E_top| = 1e-2,
#    `v_8 ≈ 1e19` at 1e-3, at every order). `ε` - the VANISHING cycle's own index - is
#    finite and order-stable in the same place, to 5%: a measured asymmetry, not an
#    assumption, and the reason the refusal names the well cycle. That window throws
#    `QuantizationError`; resumming it is a summation problem, not a
#    connection-formula one.
# 6. DERIVATIVES of the condition (`quantization_derivative`, and the Newton step).
#    Ours is the SUMMED DERIVATIVE SERIES, not the derivative of the summed condition:
#    the Riccati recursion is differentiated symbolically (`wkb_derivative`), the
#    resulting ∂S_m are integrated over the contour FROZEN at E - exact for a closed
#    cycle, which is what every spectral cycle is - and the derivative series is then
#    Borel-summed by the ordinary `_log_voros` path, since ∂_λ commutes with Borel
#    summation. The choice matters because the two are not the same at finite Padé
#    order: the Padé approximant of the derivative series is not the derivative of the
#    primal series' Padé approximant. Measured on the BARRIER cycle of V = (z²−1)² at
#    ħ = 0.25 - the worst-behaved ħ-series in the package - the gap to a central
#    difference is 4e-10 at order 4, 1.2e-3 at order 8 and 1.2e-5 at order 10, never
#    exceeding the derivative's own order-to-order movement (1.2e-3 from order 8 to
#    10); on a single well the two agree at the finite difference's own floor at every
#    order. So the ambiguity is the resummation's, and both halves of that statement
#    are tests. Newton is therefore formally quasi-Newton on the double well, which
#    costs nothing: the root is where the condition vanishes either way.
#    The alternative - carrying dual numbers through the summation - was rejected, and
#    not only on effort: `_value_type(ħ)` (ddp.jl) converts the summed value on ħ's
#    type and would silently strip a dual carried in E, and the summation itself lives
#    in Resurgence.jl.

# ── internal: cycles + shared WKB expansion + summed multipliers at one energy ──

# The condition data at energy E: spectral cycles and one WKBExpansion, reused by
# every multiplier at this E (the plan's per-(prob, E) cache).
function _condition_data(prob::SchrodingerProblem, E; order, quad_rtol = nothing)
    cycles = spectral_cycles(prob, E)
    isempty(cycles) && throw(QuantizationError(
        "no classically allowed region at E = $E (E below the potential minimum?)"))
    w = wkb_expansion(with_energy(prob, E); order)
    (cycles = cycles, w = w, quad_rtol = quad_rtol)
end

_cycle_symbol(data, c::SpectralCycle) =
    voros_symbol(data.w, c.contour; rtol = data.quad_rtol)

# median/lateral-summed exponent of a cycle's Voros symbol at ħ, θ = 0
function _summed_log(data, c::SpectralCycle, ħ; side, pade_order, tilt, rtol)
    _log_voros(_cycle_symbol(data, c), ħ; theta = 0, side, order = pade_order,
               tilt, rtol)
end

# ∂_λ of the same quantity. The derivative *series* is summed by the ordinary
# `_log_voros` path - ∂_λ commutes with Borel summation, so the summation itself is
# never differentiated and nothing in Resurgence.jl has to know about derivatives.
function _dsummed_log(data, dw::WKBDerivative, c::SpectralCycle, ħ;
                      side, pade_order, tilt, rtol)
    _log_voros(voros_symbol(dw, c.contour; rtol = data.quad_rtol), ħ;
               theta = 0, side, order = pade_order, tilt, rtol)
end

_wells(cycles) = [c for c in cycles if kind(c) === :well]
_barriers(cycles) = [c for c in cycles if kind(c) === :barrier]

function _require_layout(cycles, parity, E)
    nw, nb = length(_wells(cycles)), length(_barriers(cycles))
    if parity === nothing
        (nw == 1 && nb == 0) || throw(QuantizationError(
            "the single-well quantization condition needs exactly one well and no " *
            "barrier at E = $E (found $nw well(s), $nb barrier(s)); for a symmetric " *
            "double well pass parity = ±1"))
    else
        parity in (-1, 1) || throw(Resurgence.InvalidArgument(
            "parity must be +1 (even), -1 (odd), or nothing, got $parity"))
        (nw == 2 && nb == 1) || throw(QuantizationError(
            "the parity-factorized condition needs a double-well layout (two wells, " *
            "one barrier) at E = $E, found $nw well(s), $nb barrier(s); for a single " *
            "well omit parity"))
    end
end

# ── the uniform (Weber) branch: ledger item 5 ───────────────────────────────────

# arctan(e^{−πε}), written so the over-barrier side (ε ≪ 0) cannot overflow.
_barrier_phase(ε::Real) =
    ε ≥ 0 ? atan(exp(-oftype(ε, π) * ε)) :
            oftype(ε, π) / 2 - atan(exp(oftype(ε, π) * ε))

# (φ, ε) at one energy, on the cycles continued through the barrier top.
# (φ, ε) of the uniform condition, and - when `params` is given - their ∂_λ for each
# parameter, off one `uniform_cycles` and one `wkb_expansion`.
function _uniform_data(prob::SchrodingerProblem, E, ħ; order, side, pade_order, tilt,
                       rtol, quad_rtol, params = nothing)
    cyc = uniform_cycles(prob, E)
    w = wkb_expansion(with_energy(prob, E); order)
    sym(c, ex) = try
        voros_symbol(ex, c.contour; rtol = quad_rtol)
    catch e
        e isa Resurgence.InvalidArgument || rethrow()
        # the pinched-cycle divergence of ledger item 5 - name it rather than let the
        # generic even-period guard surface
        throw(QuantizationError(
            "the quantum period of the well cycle diverges at E = $E: it is pinched " *
            "by the vanishing cycle, so its ħ-series is useless for |E − E_top| ≲ ħ " *
            "(barrier index ε ≈ 0). Lower ħ, or step the eigenvalue past this window"))
    end
    lg(c, ex) = _log_voros(sym(c, ex), ħ; theta = 0, side, order = pade_order, tilt, rtol)
    logVB = lg(cyc.well, w)
    logVA = lg(cyc.barrier, w)
    φ = -imag(logVB) / 2
    s = cyc.below ? 1 : -1
    ε = s * abs(real(logVA)) / (2 * oftype(φ, π))                      # ledger item 3
    params === nothing && return (φ, ε)
    dws = [wkb_derivative(w; wrt, index) for (wrt, index) in params]
    dφ = [-imag(lg(cyc.well, dw)) / 2 for dw in dws]
    dε = [s * sign(real(logVA)) * real(lg(cyc.barrier, dw)) / (2 * oftype(φ, π))
          for dw in dws]
    (φ, ε, dφ, dε)
end

function _uniform_condition(prob::SchrodingerProblem, E, ħ; n, parity, side, order,
                            pade_order, tilt, rtol, quad_rtol)
    parity in (-1, 1) || throw(Resurgence.InvalidArgument(
        "the uniform condition is the parity-factorized one: parity must be +1 " *
        "(even) or −1 (odd), got $parity"))
    φ, ε = _uniform_data(prob, E, ħ; order, side, pade_order, tilt, rtol, quad_rtol)
    2φ + parity * _barrier_phase(ε) - oftype(φ, π) * (2n + 1)
end

# ── exact parameter derivatives of the condition ────────────────────────────────
#
# Both branches are functions of the summed multipliers log V_B, log V_A alone, so one
# `_dsummed_log` per cycle differentiates either of them. Nothing here differentiates
# the solver: `wkb_eigenvalue`'s Newton step consumes ∂cond/∂E, and the implicit
# function theorem on cond(E, v) = 0 turns ∂cond/∂v_k into ∂E/∂v_k.

# d/dε arctan(e^{−πε}) = −π/(2 cosh πε), written so it cannot overflow: `cosh` runs to
# Inf on the flat wings and the derivative correctly goes to zero there.
_barrier_phase_deriv(ε) = -oftype(ε, π) / (2 * cosh(oftype(ε, π) * ε))

# The value of the quantization condition together with its ∂_λ for each parameter in
# `params` (a vector of `(wrt, index)` pairs, as taken by `wkb_derivative`). One
# `spectral_cycles`, one `wkb_expansion` and one primal summation serve all of them.
function _condition_and_derivatives(prob::SchrodingerProblem, E, ħ, params;
                                    n::Integer, parity, uniform::Bool, side, order,
                                    pade_order, tilt, rtol, quad_rtol)
    uniform && return _uniform_and_derivatives(prob, E, ħ, params; n, parity, side,
                                               order, pade_order, tilt, rtol, quad_rtol)
    data = _condition_data(prob, E; order, quad_rtol)
    _require_layout(data.cycles, parity, E)
    dws = [wkb_derivative(data.w; wrt, index) for (wrt, index) in params]
    well = _wells(data.cycles)[1]
    lB = _summed_log(data, well, ħ; side, pade_order, tilt, rtol)
    dlB = [_dsummed_log(data, dw, well, ħ; side, pade_order, tilt, rtol) for dw in dws]
    if parity === nothing
        value = imag(lB) + 2 * oftype(imag(lB), π) * (n + 1 // 2)
        return value, [imag(d) for d in dlB]
    end
    barrier = _barriers(data.cycles)[1]
    lA = _summed_log(data, barrier, ħ; side, pade_order, tilt, rtol)
    dlA = [_dsummed_log(data, dw, barrier, ħ; side, pade_order, tilt, rtol) for dw in dws]
    φ = -imag(lB) / 2
    VA = exp(real(lA))
    r = sqrt(VA / (1 + VA))
    σ = parity * (iseven(n) ? 1 : -1)
    # ∂r/r = ∂log V_A / (2(1+V_A)) - the log form, so a deep barrier (V_A ≈ 0, r ≈ 0)
    # does not divide by a vanishing r.
    ders = [-sin(φ) * (-imag(db) / 2) - σ * r * real(da) / (2 * (1 + VA)) / 2
            for (db, da) in zip(dlB, dlA)]
    cos(φ) - σ * r / 2, ders
end

function _uniform_and_derivatives(prob::SchrodingerProblem, E, ħ, params;
                                  n::Integer, parity, side, order, pade_order, tilt,
                                  rtol, quad_rtol)
    parity in (-1, 1) || throw(Resurgence.InvalidArgument(
        "the uniform condition is the parity-factorized one: parity must be +1 " *
        "(even) or −1 (odd), got $parity"))
    φ, ε, dφ, dε = _uniform_data(prob, E, ħ; order, side, pade_order, tilt, rtol,
                                 quad_rtol, params)
    value = 2φ + parity * _barrier_phase(ε) - oftype(φ, π) * (2n + 1)
    value, [2 * a + parity * _barrier_phase_deriv(ε) * b for (a, b) in zip(dφ, dε)]
end

"""
    quantization_condition(prob::SchrodingerProblem, E, ħ; n::Integer,
                           parity = nothing, uniform = false, side = :median,
                           order = 12, pade_order = nothing, tilt = 1//100,
                           rtol = nothing, quad_rtol = nothing) -> Real

The exact (Voros) quantization condition at energy `E`: a real residual whose zero
in `E` is the eigenvalue. Two supported layouts (see the ledger above):

- **Single well** (`parity = nothing`): the all-orders Bohr–Sommerfeld residual
  ``\\mathrm{Im}\\log V_B^{med} + 2π(n + ½)``, exact at every order for the
  harmonic oscillator.
- **Symmetric double well** (`parity = ±1`, doublet index `n`): the median-summed
  parity condition ``\\cos φ − σ\\sqrt{V_A/(1+V_A)}``, ``σ = parity·(−1)^n``
  (`+1` = even = lower member).

With `uniform = true` (double well only) the Airy connection constant ``½`` is
replaced by the parabolic-cylinder (Weber) one and the condition becomes the phase
residual ``2φ + p\\arctan(e^{-πε}) − π(2n+1)``, `p = parity`, on the cycles
[`uniform_cycles`](@ref) continues through the barrier top - so it is defined, and
pinned against dense diagonalization, for levels **above** the barrier top, where the
non-uniform condition has no classical layout to stand on. It throws
[`QuantizationError`](@ref) in the window `|E − E_top| ≲ ħ`, where the well cycle is
pinched by the vanishing one and its ħ-series diverges (ledger item 5).

`side` picks the lateral/median summation (`:median` is canonical and mandatory on
the double-well resonant ray); `order` is the WKB order; `pade_order`/`tilt`/`rtol`
pass to the Borel–Padé summation; `quad_rtol` to the period quadratures. Precision
follows `ħ` and the potential's float type.
"""
function quantization_condition(prob::SchrodingerProblem, E, ħ; n::Integer,
                                parity = nothing, uniform::Bool = false,
                                side::Symbol = :median,
                                order::Integer = 12, pade_order = nothing,
                                tilt::Real = 1 // 100, rtol = nothing,
                                quad_rtol = nothing)
    n ≥ 0 || throw(Resurgence.InvalidArgument("the quantum number n must be ≥ 0, got $n"))
    uniform && return _uniform_condition(prob, E, ħ; n, parity, side, order, pade_order,
                                         tilt, rtol, quad_rtol)
    data = _condition_data(prob, E; order, quad_rtol)
    _require_layout(data.cycles, parity, E)
    wells = _wells(data.cycles)
    logVB = _summed_log(data, wells[1], ħ; side, pade_order, tilt, rtol)
    if parity === nothing
        return imag(logVB) + 2 * oftype(imag(logVB), π) * (n + 1 // 2)
    end
    logVA = _summed_log(data, _barriers(data.cycles)[1], ħ; side, pade_order,
                        tilt, rtol)
    φ = -imag(logVB) / 2
    VA = exp(real(logVA))
    σ = parity * (iseven(n) ? 1 : -1)
    cos(φ) - σ * sqrt(VA / (1 + VA)) / 2
end

"""
    spectral_determinant(prob::SchrodingerProblem, E, ħ; parity = nothing,
                         side = :median, order = 12, kwargs...) -> Number

The spectral determinant `D(E)` assembled from the same median-summed Voros
multipliers as [`quantization_condition`](@ref) (ledger item 4): single well
``D = 1 + V_B``; symmetric double well ``D = \\cos²φ − ¼V_A/(1+V_A)``, or the parity
factor ``D_± = \\cos φ ∓ ½\\sqrt{V_A/(1+V_A)}`` when `parity = ±1` is given.
Zeros in `E` are the eigenvalues. Thin by design - values only, no entire-function
theory.
"""
function spectral_determinant(prob::SchrodingerProblem, E, ħ; parity = nothing,
                              side::Symbol = :median, order::Integer = 12,
                              pade_order = nothing, tilt::Real = 1 // 100,
                              rtol = nothing, quad_rtol = nothing)
    data = _condition_data(prob, E; order, quad_rtol)
    wells = _wells(data.cycles)
    barriers = _barriers(data.cycles)
    logVB = _summed_log(data, wells[1], ħ; side, pade_order, tilt, rtol)
    if length(wells) == 1 && isempty(barriers)
        return 1 + exp(logVB)
    end
    (length(wells) == 2 && length(barriers) == 1) || throw(QuantizationError(
        "spectral_determinant supports a single well or a symmetric double well; " *
        "found $(length(wells)) well(s), $(length(barriers)) barrier(s) at E = $E"))
    logVA = _summed_log(data, barriers[1], ħ; side, pade_order, tilt, rtol)
    φ = -imag(logVB) / 2
    VA = exp(real(logVA))
    parity === nothing && return cos(φ)^2 - VA / (4 * (1 + VA))
    parity in (-1, 1) || throw(Resurgence.InvalidArgument(
        "parity must be +1, -1, or nothing, got $parity"))
    cos(φ) - parity * sqrt(VA / (1 + VA)) / 2
end

# ── eigenvalue solver ───────────────────────────────────────────────────────────

# minimum of V over its real critical points (the well-bottom floor for seeding)
function _potential_floor(prob::SchrodingerProblem, ::Type{F}) where {F}
    v = [F(c) for c in prob.v_coeffs]
    dv = [k * v[k + 1] for k in 1:(length(v) - 1)]
    crits = length(dv) == 1 ? Complex{F}[] :
            PolynomialRoots.roots(Complex{F}.(dv))
    zs = [real(z) for z in crits if abs(imag(z)) ≤ sqrt(eps(F)) * (1 + abs(z))]
    isempty(zs) && return F(-Inf)
    minimum(_horner(v, z) for z in zs)
end

# classical action J(E) of the seed well (first well when parity mode): cheap, no
# WKB recursion. Returns NaN when the layout does not match (used by the bracketing).
function _seed_action(prob::SchrodingerProblem, E, parity, ::Type{F}) where {F}
    cycles = try
        spectral_cycles(prob, E)
    catch e
        e isa CoalescentTurningPoints && return F(NaN)
        rethrow()
    end
    wells = _wells(cycles)
    if parity === nothing
        (length(wells) == 1 && isempty(_barriers(cycles))) || return F(NaN)
    else
        (length(wells) == 2 && length(_barriers(cycles)) == 1) || return F(NaN)
    end
    # loose quadrature: the seed only needs ~1e-3 relative accuracy
    abs(period_integral(with_energy(prob, E), wells[1].contour; rtol = 1e-8))
end

# leading Bohr–Sommerfeld seed: solve J(E) = 2πħ(n + ½) by bracketed bisection
function _seed_energy(prob::SchrodingerProblem, n, ħ, parity, ::Type{F}) where {F}
    target = 2 * F(π) * F(ħ) * (n + F(1 // 2))
    floor_E = _potential_floor(prob, F)
    isfinite(floor_E) || throw(QuantizationError(
        "cannot seed the eigenvalue: the potential has no real critical point"))
    # find a starting scale: grow δ until the layout carries a well
    δ = F(ħ) * (1 + abs(floor_E))
    lo = floor_E + δ
    J_lo = _seed_action(prob, lo, parity, F)
    tries = 0
    while isnan(J_lo)
        δ /= 8
        lo = floor_E + δ
        J_lo = _seed_action(prob, lo, parity, F)
        (tries += 1) > 60 && throw(QuantizationError(
            "cannot seed the eigenvalue: no valid spectral-cycle layout above the " *
            "potential floor $floor_E (parity = $parity)"))
    end
    J_lo ≥ target && return lo
    hi = lo
    step = max(δ, F(ħ))
    while true
        hi += step
        J_hi = _seed_action(prob, hi, parity, F)
        if isnan(J_hi)
            throw(QuantizationError(
                "seeding left the supported layout at E = $hi before reaching the " *
                "Bohr–Sommerfeld action $target (doublet n = $n above the barrier " *
                "top? - the parity-factorized condition needs the double-well layout)"))
        end
        J_hi ≥ target && break
        step *= 2
        hi > floor_E + F(10)^9 * (1 + abs(floor_E)) && throw(QuantizationError(
            "seeding did not reach the Bohr–Sommerfeld action $target"))
    end
    # bisect J(E) = target (J is increasing in E on the supported layout)
    for _ in 1:80
        mid = (lo + hi) / 2
        J_mid = _seed_action(prob, mid, parity, F)
        isnan(J_mid) && (hi = mid; continue)
        if J_mid < target
            lo = mid
        else
            hi = mid
        end
        (hi - lo) ≤ F(1 // 1000) * (1 + abs(mid)) && break
    end
    (lo + hi) / 2
end

# The barrier top: the largest value of V at a real critical point that is a local
# maximum (`NaN` when the potential has none - then there is no top to avoid).
function _barrier_top(prob::SchrodingerProblem, ::Type{F}) where {F}
    v = [F(c) for c in prob.v_coeffs]
    dv = [k * v[k + 1] for k in 1:(length(v) - 1)]
    length(dv) ≤ 1 && return F(NaN)
    d2v = [k * dv[k + 1] for k in 1:(length(dv) - 1)]
    crits = PolynomialRoots.roots(Complex{F}.(dv))
    zs = [real(z) for z in crits if abs(imag(z)) ≤ sqrt(eps(F)) * (1 + abs(z))]
    maxima = [_horner(v, z) for z in zs if _horner(d2v, z) < 0]
    isempty(maxima) ? F(NaN) : maximum(maxima)
end

# Leading-order seed for the uniform route: the same phase condition with CLASSICAL
# periods only (no WKB recursion, no Borel sum), which is monotone increasing in E and
# - unlike `_seed_action` - defined on both sides of the barrier top.
function _uniform_seed(prob::SchrodingerProblem, n, ħ, parity, ::Type{F}) where {F}
    function f1(E)
        cyc = try
            uniform_cycles(prob, E)
        catch e
            (e isa QuantizationError || e isa CoalescentTurningPoints) &&
                return F(NaN)
            rethrow()
        end
        probE = with_energy(prob, E)
        φ = -imag(period_integral(probE, cyc.well.contour; rtol = 1e-8)) / (2 * F(ħ))
        ε = (cyc.below ? 1 : -1) *
            abs(real(period_integral(probE, cyc.barrier.contour; rtol = 1e-8))) /
            (2 * F(π) * F(ħ))
        2φ + parity * _barrier_phase(ε) - F(π) * (2n + 1)
    end
    # `f` is continuous through the barrier top but its GEOMETRY is not: at E_top the
    # inner pair has merged (three turning points, and a well cycle that would enclose
    # a single branch point) and just off it the well contour has to squeeze through a
    # gap that shrinks with the separation. So the seed simply does not sample a thin
    # window around E_top; `f` is monotone, so clamping to the window's edges keeps
    # every sign and only blurs the root by `w`, far inside the seed's own tolerance.
    # Leaving a NaN there instead would make the bisection treat E_top as an upper
    # bound and converge to it, hiding every level above the barrier.
    Etop = _barrier_top(prob, F)
    w = F(1 // 10)^4 * (1 + abs(Etop))
    function f(E)
        isnan(Etop) && return f1(E)
        abs(E - Etop) ≥ w && return f1(E)
        f1(E ≥ Etop ? Etop + w : Etop - w)
    end
    floor_E = _potential_floor(prob, F)
    isfinite(floor_E) || throw(QuantizationError(
        "cannot seed the eigenvalue: the potential has no real critical point"))
    δ = F(ħ) * (1 + abs(floor_E))
    lo, f_lo, tries = floor_E + δ, F(NaN), 0
    while isnan(f_lo)
        lo = floor_E + δ
        f_lo = f(lo)
        δ /= 8
        (tries += 1) > 60 && throw(QuantizationError(
            "cannot seed the uniform eigenvalue: no four-turning-point double-well " *
            "layout above the potential floor $floor_E"))
    end
    f_lo ≥ 0 && return lo
    hi, step = lo, max(δ, F(ħ))
    while true
        hi += step
        f_hi = f(hi)
        if isnan(f_hi)                          # stepped onto E_top itself
            hi -= step
            step /= 2
            step ≤ eps(F) * (1 + abs(hi)) && throw(QuantizationError(
                "cannot bracket the uniform eigenvalue above E = $hi"))
            continue
        end
        f_hi ≥ 0 && break
        step *= 2
        hi > floor_E + F(10)^9 * (1 + abs(floor_E)) && throw(QuantizationError(
            "seeding did not reach the uniform quantization phase for n = $n"))
    end
    for _ in 1:80
        mid = (lo + hi) / 2
        f_mid = f(mid)
        isnan(f_mid) && (hi = mid; continue)
        f_mid < 0 ? (lo = mid) : (hi = mid)
        (hi - lo) ≤ F(1 // 1000) * (1 + abs(mid)) && break
    end
    (lo + hi) / 2
end

"""
    wkb_eigenvalue(prob::SchrodingerProblem, n::Integer, ħ; parity = nothing,
                   uniform = false, side = :median, order = 12, rtol = nothing,
                   maxiter = 30, kwargs...) -> Real

The `n`-th eigenvalue of ``−ħ²ψ″ + Vψ = Eψ`` from the exact quantization condition:
Newton's method in `E` on [`quantization_condition`](@ref) with the **exact**
derivative of [`quantization_derivative`](@ref) (one evaluation per step, no step
size), seeded by leading-order Bohr–Sommerfeld (bracketed bisection on the classical
action). For a symmetric double well pass
`parity = +1` (even, the lower member of doublet `n`) or `−1` (odd); `n` is then
the doublet index. `rtol` is the relative tolerance on `E` (default `√eps` of the
working float type - set it explicitly for BigFloat runs); remaining keywords pass
to [`quantization_condition`](@ref). Throws [`QuantizationError`](@ref) on
non-convergence.

`uniform = true` switches to the Weber (parabolic-cylinder) condition and to a seed
built from the same phase, so `n` counts doublets **through and above the barrier
top**: level `2n` of the whole well for `parity = +1` and `2n + 1` for `parity = −1`
once the barrier is submerged. The default route cannot reach those levels at all.
"""
function wkb_eigenvalue(prob::SchrodingerProblem, n::Integer, ħ; parity = nothing,
                        uniform::Bool = false,
                        side::Symbol = :median, order::Integer = 12,
                        rtol = nothing, maxiter::Integer = 30, pade_order = nothing,
                        tilt::Real = 1 // 100, sum_rtol = nothing, quad_rtol = nothing)
    n ≥ 0 || throw(Resurgence.InvalidArgument("the quantum number n must be ≥ 0, got $n"))
    F = promote_type(_wkb_float(prob), typeof(float(real(ħ))))
    ħF = F(real(ħ))
    tol = rtol === nothing ? sqrt(eps(F)) : F(rtol)
    floor_E = _potential_floor(prob, F)
    E = uniform ? _uniform_seed(prob, n, ħF, parity, F) :
                  _seed_energy(prob, n, ħF, parity, F)
    # One evaluation per Newton step returns both the residual and its exact ∂/∂E: the
    # derivative series is summed alongside the primal one off a single
    # `spectral_cycles` + `wkb_expansion`, so this costs less than the three
    # `quantization_condition` calls a central difference needed - and there is no `h`.
    function cond(x)
        r, d = _condition_and_derivatives(prob, x, ħF, ((:energy, 0),); n, parity,
                                          uniform, side, order, pade_order, tilt,
                                          rtol = sum_rtol, quad_rtol)
        r, d[1]
    end
    scale(x) = 1 + abs(x)
    # a step must land strictly inside the supported layout: above the potential
    # floor (below it there is no allowed region), and evaluable without hitting a
    # coalescent/critical energy. The double-well condition oscillates in E, so a
    # leading-order seed can sit on the wrong slope - backtrack such steps.
    # The step-validity probe returns the data it computed, so an accepted step needs
    # no second evaluation at the new energy.
    function step_to(E, ΔE)
        for _ in 1:60
            En = E - ΔE
            if isfinite(floor_E) && En ≤ floor_E + sqrt(eps(F)) * (1 + abs(floor_E))
                ΔE /= 2; continue
            end
            got = try
                cond(En)
            catch e
                (e isa CoalescentTurningPoints || e isa QuantizationError) || rethrow()
                nothing
            end
            got === nothing || return En, ΔE, got[1], got[2]
            ΔE /= 2
        end
        throw(QuantizationError(
            "Newton cannot take a valid step from E ≈ $E (parity = $parity); the " *
            "requested level may lie outside the supported layout"))
    end
    r, drdE = cond(E)
    for iter in 1:maxiter
        iszero(drdE) && throw(QuantizationError(
            "flat quantization condition at E = $E (iteration $iter)"))
        ΔE = r / drdE
        # damp absurd steps (bad seed or near-critical derivative)
        maxstep = scale(E) / 4
        abs(ΔE) > maxstep && (ΔE = sign(ΔE) * maxstep)
        E, ΔEtaken, r, drdE = step_to(E, ΔE)
        abs(ΔEtaken) ≤ tol * scale(E) && return E
    end
    throw(QuantizationError(
        "Newton did not converge to rtol = $tol in $maxiter iterations (last E = $E); " *
        "loosen rtol, raise maxiter, or check n/parity against the layout"))
end

"""
    quantization_derivative(prob::SchrodingerProblem, E, ħ; wrt = :energy, index = 0,
                            n, parity = nothing, uniform = false, kwargs...) -> Real

The exact derivative of [`quantization_condition`](@ref) with respect to one parameter
of the potential: `wrt = :energy` gives ``∂/∂E``, and `wrt = :coefficient` with
`index = k` gives ``∂/∂v_k``. Remaining keywords are those of
[`quantization_condition`](@ref).

Exact in the sense that the WKB recursion is differentiated symbolically and the
resulting derivative series is Borel-summed by the ordinary route - no step size, and
no finite difference anywhere. The contour is frozen at `E`, which costs nothing: a
closed cycle's period does not depend on its representative.

What that does **not** mean: the Padé approximant of the derivative series is not the
derivative of the primal series' Padé approximant, so this differs from
differentiating the summed condition by the resummation's own ambiguity. Measured on
the symmetric double well at `ħ = 0.25`, where the barrier cycle's asymptotics is at
its worst, the two routes agree to 4e-10 at `order = 4` and to 1e-5 at `order = 10`,
in both cases within the order-to-order movement of the condition itself. On a single
well they agree to the finite difference's own error at every order.
"""
function quantization_derivative(prob::SchrodingerProblem, E, ħ; wrt::Symbol = :energy,
                                 index::Integer = 0, n::Integer, parity = nothing,
                                 uniform::Bool = false, side::Symbol = :median,
                                 order::Integer = 12, pade_order = nothing,
                                 tilt::Real = 1 // 100, rtol = nothing,
                                 quad_rtol = nothing)
    n ≥ 0 || throw(Resurgence.InvalidArgument("the quantum number n must be ≥ 0, got $n"))
    _, d = _condition_and_derivatives(prob, E, ħ, ((wrt, index),); n, parity, uniform,
                                      side, order, pade_order, tilt, rtol, quad_rtol)
    d[1]
end

"""
    eigenvalue_sensitivity(prob::SchrodingerProblem, n::Integer, ħ; kwargs...)
        -> NamedTuple

The eigenvalue and its gradient in the potential coefficients:
`(E, dE_dv, dcond_dE)`, where `dE_dv[k+1] = ∂E/∂v_k` for the coefficient of ``z^k`` in
`V`. By the implicit function theorem on the quantization condition
``F(E, v) = 0``: ``∂E/∂v_k = −(∂F/∂v_k)/(∂F/∂E)``, with every partial derivative exact
(see [`quantization_derivative`](@ref)). Keywords are those of
[`wkb_eigenvalue`](@ref).

`dcond_dE` is returned because it is the conditioning of the whole answer: a level
approaching a degeneracy sends it to zero and the sensitivities to infinity, and that
is real, not numerical.
"""
function eigenvalue_sensitivity(prob::SchrodingerProblem, n::Integer, ħ;
                                parity = nothing, uniform::Bool = false,
                                side::Symbol = :median, order::Integer = 12,
                                rtol = nothing, maxiter::Integer = 30,
                                pade_order = nothing, tilt::Real = 1 // 100,
                                sum_rtol = nothing, quad_rtol = nothing)
    E = wkb_eigenvalue(prob, n, ħ; parity, uniform, side, order, rtol, maxiter,
                       pade_order, tilt, sum_rtol, quad_rtol)
    F = promote_type(_wkb_float(prob), typeof(float(real(ħ))))
    d = length(q_coefficients(prob)) - 1
    params = vcat([(:energy, 0)], [(:coefficient, k) for k in 0:d])
    _, ders = _condition_and_derivatives(prob, E, F(real(ħ)), params; n, parity,
                                         uniform, side, order, pade_order, tilt,
                                         rtol = sum_rtol, quad_rtol)
    dFdE = ders[1]
    iszero(dFdE) && throw(QuantizationError(
        "the quantization condition is flat at E = $E, so the level's sensitivity to " *
        "the potential is not defined there (a degeneracy, or too loose a summation)"))
    (E = E, dE_dv = [-c / dFdE for c in ders[2:end]], dcond_dE = dFdE)
end
