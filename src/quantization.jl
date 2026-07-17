# Exact (Voros) quantization conditions, spectral determinants, and the eigenvalue
# solver. The conditions are built from median-summed Voros multipliers (`_log_voros`,
# the DDP layer's summation — never re-implemented here) of the spectral cycles at
# energy E; the eigenvalues are the Newton roots in E.
#
# ── M5 quantization ledger (each item pinned by an oracle in test_quantization /
#    test_double_well) ────────────────────────────────────────────────────────────
#
# 1. Well multiplier. With the M5 cycle orientation (quantum_periods ledger) the
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
#    parity = +1 the even (symmetric) member of doublet n — always the lower one —
#    and −1 the odd. The factor ½ is the connection-formula tunneling amplitude
#    (the Airy 1/2 that a single barrier contributes) — PINNED by the oracle: with
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
#    regularization (deliberately out of the M5 scope).

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

"""
    quantization_condition(prob::SchrodingerProblem, E, ħ; n::Integer,
                           parity = nothing, side = :median, order = 12,
                           pade_order = nothing, tilt = 1//100, rtol = nothing,
                           quad_rtol = nothing) -> Real

The exact (Voros) quantization condition at energy `E`: a real residual whose zero
in `E` is the eigenvalue. Two supported layouts (see the ledger above):

- **Single well** (`parity = nothing`): the all-orders Bohr–Sommerfeld residual
  ``\\mathrm{Im}\\log V_B^{med} + 2π(n + ½)``, exact at every order for the
  harmonic oscillator.
- **Symmetric double well** (`parity = ±1`, doublet index `n`): the median-summed
  parity condition ``\\cos φ − σ\\sqrt{V_A/(1+V_A)}``, ``σ = parity·(−1)^n``
  (`+1` = even = lower member).

`side` picks the lateral/median summation (`:median` is canonical and mandatory on
the double-well resonant ray); `order` is the WKB order; `pade_order`/`tilt`/`rtol`
pass to the Borel–Padé summation; `quad_rtol` to the period quadratures. Precision
follows `ħ` and the potential's float type.
"""
function quantization_condition(prob::SchrodingerProblem, E, ħ; n::Integer,
                                parity = nothing, side::Symbol = :median,
                                order::Integer = 12, pade_order = nothing,
                                tilt::Real = 1 // 100, rtol = nothing,
                                quad_rtol = nothing)
    n ≥ 0 || throw(Resurgence.InvalidArgument("the quantum number n must be ≥ 0, got $n"))
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
Zeros in `E` are the eigenvalues. Thin by design — values only, no entire-function
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
                "top? — the parity-factorized condition needs the double-well layout)"))
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

"""
    wkb_eigenvalue(prob::SchrodingerProblem, n::Integer, ħ; parity = nothing,
                   side = :median, order = 12, rtol = nothing, maxiter = 30,
                   kwargs...) -> Real

The `n`-th eigenvalue of ``−ħ²ψ″ + Vψ = Eψ`` from the exact quantization condition:
Newton's method in `E` (finite-difference derivative) on
[`quantization_condition`](@ref), seeded by leading-order Bohr–Sommerfeld
(bracketed bisection on the classical action). For a symmetric double well pass
`parity = +1` (even, the lower member of doublet `n`) or `−1` (odd); `n` is then
the doublet index. `rtol` is the relative tolerance on `E` (default `√eps` of the
working float type — set it explicitly for BigFloat runs); remaining keywords pass
to [`quantization_condition`](@ref). Throws [`QuantizationError`](@ref) on
non-convergence.
"""
function wkb_eigenvalue(prob::SchrodingerProblem, n::Integer, ħ; parity = nothing,
                        side::Symbol = :median, order::Integer = 12,
                        rtol = nothing, maxiter::Integer = 30, pade_order = nothing,
                        tilt::Real = 1 // 100, sum_rtol = nothing, quad_rtol = nothing)
    n ≥ 0 || throw(Resurgence.InvalidArgument("the quantum number n must be ≥ 0, got $n"))
    F = promote_type(_wkb_float(prob), typeof(float(real(ħ))))
    ħF = F(real(ħ))
    tol = rtol === nothing ? sqrt(eps(F)) : F(rtol)
    floor_E = _potential_floor(prob, F)
    E = _seed_energy(prob, n, ħF, parity, F)
    cond(x) = quantization_condition(prob, x, ħF; n, parity, side, order,
                                     pade_order, tilt, rtol = sum_rtol, quad_rtol)
    scale(x) = 1 + abs(x)
    # a step must land strictly inside the supported layout: above the potential
    # floor (below it there is no allowed region), and evaluable without hitting a
    # coalescent/critical energy. The double-well condition oscillates in E, so a
    # leading-order seed can sit on the wrong slope — backtrack such steps.
    function step_to(E, ΔE)
        for _ in 1:60
            En = E - ΔE
            if isfinite(floor_E) && En ≤ floor_E + sqrt(eps(F)) * (1 + abs(floor_E))
                ΔE /= 2; continue
            end
            ok = try
                cond(En); true
            catch e
                (e isa CoalescentTurningPoints || e isa QuantizationError) || rethrow()
                false
            end
            ok && return En, ΔE
            ΔE /= 2
        end
        throw(QuantizationError(
            "Newton cannot take a valid step from E ≈ $E (parity = $parity); the " *
            "requested level may lie outside the supported layout"))
    end
    for iter in 1:maxiter
        h = max(cbrt(eps(F)) * scale(E), 8 * tol * scale(E))
        r, rp, rm = cond(E), cond(E + h), cond(E - h)
        drdE = (rp - rm) / (2h)
        iszero(drdE) && throw(QuantizationError(
            "flat quantization condition at E = $E (iteration $iter)"))
        ΔE = r / drdE
        # damp absurd steps (bad seed or near-critical derivative)
        maxstep = scale(E) / 4
        abs(ΔE) > maxstep && (ΔE = sign(ΔE) * maxstep)
        Enew, ΔEtaken = step_to(E, ΔE)
        E = Enew
        abs(ΔEtaken) ≤ tol * scale(E) && return E
    end
    throw(QuantizationError(
        "Newton did not converge to rtol = $tol in $maxiter iterations (last E = $E); " *
        "loosen rtol, raise maxiter, or check n/parity against the layout"))
end
