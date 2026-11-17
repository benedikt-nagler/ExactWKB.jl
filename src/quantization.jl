# Single-well numerical quantization only.
# Ledger: well orientation gives log(V) = -i J/hbar; harmonic E_n = hbar*(2n+1).
# For anharmonic wells this finite Borel-Pade condition is an approximation.
# Its root residual does not bound the error relative to the Schrodinger spectrum.

_wells(cycles) = [c for c in cycles if kind(c) === :well]
_barriers(cycles) = [c for c in cycles if kind(c) === :barrier]

function _require_single_well(prob::SchrodingerProblem)
    _require_real_potential(prob, "quantization")
    v = real.(q_coefficients(with_energy(prob, 0)))
    d = length(v) - 1
    supported = (d == 2 && v[3] > 0) ||
                (d == 4 && v[5] > 0 && v[3] >= 0 && iszero(v[2]) && iszero(v[4]))
    supported || throw(QuantizationError(
        "quantization supports confining quadratics and V0 + a*z^2 + b*z^4 with a >= 0, b > 0"))
end

function _single_symbol(prob, E; order, quad_rtol)
    _require_single_well(prob)
    cycles = spectral_cycles(prob, E)
    length(cycles) == 1 && kind(only(cycles)) === :well ||
        throw(QuantizationError("energy must have one allowed interval with two simple real endpoints"))
    w = wkb_expansion(with_energy(prob, E); order)
    c = only(cycles)
    w, c, voros_symbol(w, c.contour; rtol = quad_rtol)
end

function _quantization_args(n, hbar)
    n >= 0 || throw(Resurgence.InvalidArgument("n must be nonnegative"))
    isreal(hbar) && isfinite(hbar) && real(hbar) > 0 ||
        throw(Resurgence.InvalidArgument("hbar must be finite, real and positive"))
end

"""
    quantization_condition(prob, E, ħ; n, order = 12, kwargs...)

Evaluate the single-well phase residual `imag(log(V)) + 2pi*(n+1/2)`.
Supports confining quadratics and symmetric quartics with nonnegative quadratic term.
Only the harmonic case is exact. Finite-order anharmonic residuals are approximations.
Throws [`QuantizationError`](@ref) outside this domain.
"""
function quantization_condition(prob::SchrodingerProblem, E, ħ; n::Integer,
                                side::Symbol = :median, order::Integer = 12,
                                pade_order = nothing, tilt::Real = 1//100,
                                rtol = nothing, quad_rtol = nothing)
    _quantization_args(n, ħ)
    _, _, vs = _single_symbol(prob, E; order, quad_rtol)
    z = _log_voros(vs, ħ; side, order = pade_order, tilt, rtol)
    imag(z) + 2 * oftype(imag(z), pi) * (n + 1//2)
end

"""
    spectral_determinant(prob, E, ħ; kwargs...)

Evaluate `1 + V` for the supported single-well approximation.
This is a numerical quantization residual, not a regularized spectral determinant.
Its anharmonic zeros need not equal the exact eigenvalues.
"""
function spectral_determinant(prob::SchrodingerProblem, E, ħ;
                              side::Symbol = :median, order::Integer = 12,
                              pade_order = nothing, tilt::Real = 1//100,
                              rtol = nothing, quad_rtol = nothing)
    _quantization_args(0, ħ)
    _, _, vs = _single_symbol(prob, E; order, quad_rtol)
    1 + voros_value(vs, ħ; side, order = pade_order, tilt, rtol)
end

function _potential_floor(prob::SchrodingerProblem, ::Type{F}) where {F}
    v = [F(c) for c in prob.v_coeffs]
    dv = [k * v[k + 1] for k in 1:(length(v) - 1)]
    crits = length(dv) == 1 ? Complex{F}[] :
            PolynomialRoots.roots(Complex{F}.(dv))
    zs = [real(z) for z in crits if abs(imag(z)) ≤ sqrt(eps(F)) * (1 + abs(z))]
    isempty(zs) && return F(-Inf)
    minimum(_horner(v, z) for z in zs)
end

function _seed_action(prob::SchrodingerProblem, E, ::Type{F}) where {F}
    cycles = try
        spectral_cycles(prob, E)
    catch e
        e isa CoalescentTurningPoints && return F(NaN)
        rethrow()
    end
    wells = _wells(cycles)
    (length(wells) == 1 && isempty(_barriers(cycles))) || return F(NaN)
    # loose quadrature: the seed only needs ~1e-3 relative accuracy
    abs(period_integral(with_energy(prob, E), wells[1].contour; rtol = 1e-8))
end

# leading Bohr–Sommerfeld seed: solve J(E) = 2πħ(n + ½) by bracketed bisection
function _seed_energy(prob::SchrodingerProblem, n, ħ, ::Type{F}) where {F}
    target = 2 * F(π) * F(ħ) * (n + F(1 // 2))
    floor_E = _potential_floor(prob, F)
    isfinite(floor_E) || throw(QuantizationError(
        "cannot seed the eigenvalue: the potential has no real critical point"))
    # find a starting scale: grow δ until the layout carries a well
    δ = F(ħ) * (1 + abs(floor_E))
    lo = floor_E + δ
    J_lo = _seed_action(prob, lo, F)
    tries = 0
    while isnan(J_lo)
        δ /= 8
        lo = floor_E + δ
        J_lo = _seed_action(prob, lo, F)
        (tries += 1) > 60 && throw(QuantizationError(
            "cannot seed the eigenvalue: no valid spectral-cycle layout above the " *
            "potential floor $floor_E"))
    end
    J_lo ≥ target && return lo
    hi = lo
    step = max(δ, F(ħ))
    while true
        hi += step
        J_hi = _seed_action(prob, hi, F)
        if isnan(J_hi)
            throw(QuantizationError(
                "seeding left the supported layout at E = $hi before reaching the " *
                "Bohr–Sommerfeld action $target"))
        end
        J_hi ≥ target && break
        step *= 2
        hi > floor_E + F(10)^9 * (1 + abs(floor_E)) && throw(QuantizationError(
            "seeding did not reach the Bohr–Sommerfeld action $target"))
    end
    # bisect J(E) = target (J is increasing in E on the supported layout)
    for _ in 1:80
        mid = (lo + hi) / 2
        J_mid = _seed_action(prob, mid, F)
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
    wkb_eigenvalue(prob, n, ħ; order = 12, rtol = nothing, maxiter = 30, kwargs...)

Find a zero of the supported single-well phase condition.
`rtol` controls solver termination, not the error relative to the exact spectrum.
The iteration sums a derivative series; it is a quasi-Newton method at finite Padé order.
Throws [`QuantizationError`](@ref) when the residual and step criteria are not both met.
"""
function wkb_eigenvalue(prob::SchrodingerProblem, n::Integer, ħ;
                        side::Symbol = :median, order::Integer = 12,
                        rtol = nothing, maxiter::Integer = 30, pade_order = nothing,
                        tilt::Real = 1//100, sum_rtol = nothing, quad_rtol = nothing)
    _quantization_args(n, ħ)
    _require_single_well(prob)
    F = promote_type(_wkb_float(prob), typeof(float(real(ħ))))
    h = F(real(ħ))
    tol = rtol === nothing ? sqrt(eps(F)) : F(rtol)
    isfinite(tol) && tol > 0 && maxiter > 0 ||
        throw(Resurgence.InvalidArgument("rtol and maxiter must be positive"))
    floor_E = _potential_floor(prob, F)
    E = _seed_energy(prob, n, h, F)
    function cond(x)
        w, c, vs = _single_symbol(prob, x; order, quad_rtol)
        z = _log_voros(vs, h; side, order = pade_order, tilt, rtol = sum_rtol)
        dw = wkb_derivative(w; wrt = :energy)
        dv = voros_symbol(dw, c.contour; rtol = quad_rtol)
        dz = _log_voros(dv, h; side, order = pade_order, tilt, rtol = sum_rtol)
        imag(z) + 2 * F(pi) * (n + 1//2), imag(dz)
    end
    residual, slope = cond(E)
    for _ in 1:maxiter
        isfinite(residual) && isfinite(slope) && !iszero(slope) ||
            throw(QuantizationError("nonfinite or flat phase condition"))
        step = clamp(residual / slope, -(1 + abs(E))/4, (1 + abs(E))/4)
        accepted = false
        for _ in 1:60
            next = E - step
            if next > floor_E
                result = try
                    cond(next)
                catch e
                    e isa CoalescentTurningPoints || e isa QuantizationError || rethrow()
                    nothing
                end
                if result !== nothing && abs(result[1]) <= abs(residual) + eps(F)
                    E, residual, slope = next, result[1], result[2]
                    accepted = true
                    break
                end
            end
            step /= 2
        end
        accepted || throw(QuantizationError("no residual-reducing step within the supported well"))
        abs(step) <= tol*(1 + abs(E)) && abs(residual) <= tol*2*F(pi)*(n+1) && return E
    end
    throw(QuantizationError("single-well phase solver did not converge"))
end
