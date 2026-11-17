# Quantum periods: contour integrals ∮ S_m dz along piecewise-linear paths, with the
# √Q branch tracked by unwrapping arg Q(z) along the contour. The classical period
# ∮√Q dz is the m = -1 case. Branch bookkeeping is the top correctness risk, so it is
# isolated here: an adaptive pre-walk lays down nodes dense enough that arg Q never
# jumps by more than ~π/2 between neighbours, and every quadrature sample unwraps
# against the nearest node - no discrete sign choices inside the integrand.
#
# ── convention ledger ──────────────────────────────────────────────────────────────
#
# 1. Parameter derivatives (`period_derivative`, `wkb_period` of a `WKBDerivative`) are
#    taken **at frozen contour**: the vertices, the branch table and the quadrature
#    breakpoints are those of the base point, and only the integrand is differentiated.
#    For a CLOSED cycle this is exact, not an approximation - the period depends on the
#    contour only through its homology class, so the moving turning points contribute
#    nothing and there is no boundary term. Pinned by harmonic derivatives, coefficient differences and Cauchy-Riemann
#    checks in test_periods.jl. These are numerical comparisons, not error bounds.
# 2. Consequently the derivative is **refused on an open contour**. On the
#    turning-point-to-turning-point form the endpoints are branch points that move with
#    the parameter, so a frozen contour drops a boundary term; `closed = false` is a
#    usage error, not a supported mode.

using QuadGK: quadgk

# principal-branch-safe wrap of an angle into (−π, π] (type-generic: uses BigFloat π
# for BigFloat arguments, so branch unwrapping keeps full precision)
function _wrap(θ)
    τ = 2 * oftype(θ, π)
    θ - τ * round(θ / τ)
end

# the float type carried by a contour / turning-point list
_contour_float(::Type{Complex{F}}) where {F} = F
_contour_float(::Type{F}) where {F<:AbstractFloat} = F
_contour_float(::Type{<:Real}) = BigFloat
_contour_float(v::AbstractVector) = _contour_float(eltype(v))

_as_complex(z::Number) = complex(z)
_as_complex(tp::TurningPoint) = tp.z

# Map global parameter s ∈ [0,1] to a point on the piecewise-linear contour and its
# derivative dz/ds (piecewise constant per segment).
function _contour_point(verts, K, s)
    seg = clamp(ceil(Int, s * K), 1, K)
    t = s * K - (seg - 1)
    z = verts[seg] + t * (verts[seg + 1] - verts[seg])
    dz = K * (verts[seg + 1] - verts[seg])
    z, dz
end

# A disk bound for the polynomial image of a segment about its midpoint.
# If the disk excludes zero, its argument stays within asin(radius/abs(center)).
# For rational Q the denominator contributes the sum of the pole-factor bounds.
# [IN14] 2.3 uses a continued square root; this bound prevents hidden full turns.
function _segment_phase_bound(prob, z0, z1)
    mid = (z0 + z1) / 2
    h = abs(z1 - z0) / 2
    c = _taylor_shift(_turning_polynomial(prob), mid)
    radius = sum(abs(c[k]) * h^(k - 1) for k in 2:length(c); init = zero(h))
    center = abs(c[1])
    radius < center || return oftype(h, Inf)
    bound = asin(radius / center)
    for (pole, m) in zip(poles(prob), pole_orders(prob))
        d = abs(mid - pole)
        h < d || return oftype(h, Inf)
        bound += m * asin(h / d)
    end
    bound
end

function _connect!(prob, s0, z0, pa0, uw0, s1, z1, S, PA, UW, atol, depth)
    Q1 = prob(z1)
    isfinite(Q1) && abs(Q1) > atol || throw(ContourError(
        "contour reaches a zero, pole or unresolved value near z = $z1"))
    if _segment_phase_bound(prob, z0, z1) >= oftype(s0, pi) / 4
        depth > 0 || throw(ContourError("cannot resolve the square-root branch on a segment"))
        sm, zm = (s0 + s1) / 2, (z0 + z1) / 2
        (sm != s0 && sm != s1) || throw(ContourError("contour subdivision reached working precision"))
        pam, uwm = _connect!(prob, s0, z0, pa0, uw0, sm, zm, S, PA, UW, atol, depth - 1)
        return _connect!(prob, sm, zm, pam, uwm, s1, z1, S, PA, UW, atol, depth - 1)
    end
    pa1 = angle(Q1)
    uw1 = uw0 + _wrap(pa1 - pa0)
    push!(S, s1); push!(PA, pa1); push!(UW, uw1)
    pa1, uw1
end

# Build the branch-tracking node table (sorted global s, principal args, unwrapped
# args) over the contour. Also returns the net unwrapped Δ over the whole path.
function _branch_table(prob, verts::Vector{Complex{F}}, K, atol; maxdepth = 40) where {F}
    Q0 = prob(verts[1])
    (!isfinite(Q0) || abs(Q0) ≤ atol) &&
        throw(ContourError("contour starts at a turning point (z ≈ $(verts[1]))"))
    pa0 = angle(Q0)
    S = F[zero(F)]; PA = [pa0]; UW = [pa0]
    pa, uw = pa0, pa0
    for k in 1:K
        pa, uw = _connect!(prob, F(k - 1) / K, verts[k], pa, uw, F(k) / K, verts[k + 1],
                           S, PA, UW, atol, maxdepth)
    end
    S, PA, UW, uw - pa0
end

# Continuous √Q at global parameter s, unwrapping against the nearest table node.
function _sqrt_Q(prob, s, S, PA, UW, verts, K)
    z, dz = _contour_point(verts, K, s)
    Q = prob(z)
    i = clamp(searchsortedlast(S, s), 1, length(S))
    uw = UW[i] + _wrap(angle(Q) - PA[i])
    sqrt(abs(Q)) * cis(uw / 2), z, Q, dz
end

# Geometric exclusion scale; independent of the requested quadrature tolerance.
function _period_atol(prob, verts, F)
    isempty(verts) && return sqrt(eps(F))
    mags = sort!([F(abs(prob(z))) for z in verts])
    eps(F) * mags[cld(length(mags), 2)]
end

# Gauss–Kronrod order for the period quadratures. The integrands are analytic on the
# contour, so raising the order converges exponentially - at BigFloat tolerances this
# beats adaptive bisection at the default order 7 by well over an order of magnitude
# (the high-m 1/Q^k integrands are stiff near the turning points).
function _quad_order(::Type{F}) where {F<:AbstractFloat}
    F === Float64 && return 7
    clamp(Int(cld(precision(F), 8)), 21, 61)
end

"""
    period_integral(prob::AbstractSchrodingerProblem, contour; closed = true) -> Complex

The classical period ``∮ √Q\\,dz`` (``= v_{-1}``) along the piecewise-linear
`contour` (a vector of complex vertices). With `closed = true` the contour is closed
back to its first vertex. The `√Q` branch is tracked continuously; a vertex sitting on
a turning point, or a contour enclosing an odd number of turning points (so `√Q` is
not single-valued), throws [`ContourError`](@ref). Precision follows the contour's
float type. `atol` and `rtol` control quadrature only.
An unmet quadrature error target raises [`ContourError`](@ref).
"""
function period_integral(prob::AbstractSchrodingerProblem, contour::AbstractVector;
                         closed::Bool = true, rtol = nothing, atol = nothing, maxevals::Integer = 10^6)
    F = _contour_float(_as_complex.(contour))
    verts = Complex{F}[_as_complex(z) for z in contour]
    length(verts) >= 2 || throw(ContourError("contour needs at least two vertices"))
    all(isfinite, verts) || throw(ContourError("contour vertices must be finite"))
    closed && push!(verts, verts[1])
    K = length(verts) - 1
    K ≥ 1 || throw(ContourError("contour needs ≥ 2 vertices"))
    geom_tol = _period_atol(prob, verts, F)
    S, PA, UW, Δ = _branch_table(prob, verts, K, geom_tol)
    if closed && isodd(round(Int, Δ / (2π)))
        throw(ContourError("closed contour encloses an odd number of turning points; " *
                           "√Q is not single-valued around it (branch ambiguity)"))
    end
    function f(s)
        u, _, _, dz = _sqrt_Q(prob, s, S, PA, UW, verts, K)
        u * dz
    end
    rt = rtol === nothing ? sqrt(eps(F)) : F(rtol)
    at = atol === nothing ? zero(F) : F(atol)
    isfinite(rt) && isfinite(at) && rt >= 0 && at >= 0 && maxevals > 0 ||
        throw(Resurgence.InvalidArgument("quadrature tolerances must be finite and nonnegative; maxevals must be positive"))
    val, err = quadgk(f, S...;
                    rtol = rt, atol = at, maxevals = maxevals, order = _quad_order(F))
    err <= max(rt * abs(val), at) || throw(ContourError("period quadrature did not reach its requested tolerance"))
    val
end

# Ledger item 2: a frozen contour is exact only around a closed cycle.
_require_closed_for_derivative(::WKBExpansion, closed) = nothing
_require_closed_for_derivative(::WKBDerivative, closed) =
    closed || throw(ContourError(
        "a parameter derivative is taken at frozen contour, which is exact only for a " *
        "closed cycle; on an open contour the turning-point endpoints move with the " *
        "parameter and the dropped boundary term is not small"))

"""
    wkb_period(w::WKBExpansion, contour, m; closed = true) -> Complex
    wkb_period(dw::WKBDerivative, contour, m; closed = true) -> Complex

The order-`m` quantum period ``∮ S_m\\,dz`` along `contour`, using the same √Q branch
tracking as [`period_integral`](@ref). `m = -1` reproduces the classical period.

Given a [`WKBDerivative`](@ref) instead, the same integral of ``∂S_m`` - that is,
``∂_λ`` of the order-`m` period, taken **at frozen contour** (see the ledger above).
"""
function wkb_period(w::Union{WKBExpansion,WKBDerivative}, contour::AbstractVector,
                    m::Integer; closed::Bool = true, rtol = nothing, atol = nothing,
                    maxevals::Integer = 10^6)
    -1 ≤ m ≤ order(w) ||
        throw(Resurgence.InvalidArgument("m must be in -1:$(order(w)), got $m"))
    _require_closed_for_derivative(w, closed)
    prob = _wkb_prob(w)
    F = _contour_float(_as_complex.(contour))
    verts = Complex{F}[_as_complex(z) for z in contour]
    length(verts) >= 2 || throw(ContourError("contour needs at least two vertices"))
    all(isfinite, verts) || throw(ContourError("contour vertices must be finite"))
    closed && push!(verts, verts[1])
    K = length(verts) - 1
    K >= 1 || throw(ContourError("contour needs at least two vertices"))
    geom_tol = _period_atol(prob, verts, F)
    S, PA, UW, Δ = _branch_table(prob, verts, K, geom_tol)
    if closed && isodd(round(Int, Δ / (2π)))
        throw(ContourError("closed contour encloses an odd number of turning points; " *
                           "√Q is not single-valued around it (branch ambiguity)"))
    end
    term = _s_term(w, m)
    function f(s)
        u, z, Q, dz = _sqrt_Q(prob, s, S, PA, UW, verts, K)
        _eval_term(term, Complex{F}(z), Q, u) * dz
    end
    rt = rtol === nothing ? sqrt(eps(F)) : F(rtol)
    at = atol === nothing ? zero(F) : F(atol)
    isfinite(rt) && isfinite(at) && rt >= 0 && at >= 0 && maxevals > 0 ||
        throw(Resurgence.InvalidArgument("quadrature tolerances must be finite and nonnegative; maxevals must be positive"))
    val, err = quadgk(f, S...;
                    rtol = rt, atol = at, maxevals = maxevals,
                    order = _quad_order(F))
    err <= max(rt * abs(val), at) || throw(ContourError("period quadrature did not reach its requested tolerance (m=$m, estimated error=$err, target=$(max(rt*abs(val),at)))"))
    val
end

"""
    period_derivative(prob::AbstractSchrodingerProblem, contour, dQ; closed = true)
        -> Complex

``∂_λ`` of the classical period ``∮√Q\\,dz`` along `contour`, for any parameter `λ`
entering only through `Q`: `dQ` is a callable giving ``∂_λ Q(z)``, and the integral is
``∮ ∂_λQ / (2√Q)\\,dz`` over the **same** contour and branch table as
[`period_integral`](@ref).

Frozen contour, so this is exact only for a **closed** cycle - see the ledger at the
top of this file. `dQ = z -> -1` is ``∂/∂E`` for a [`SchrodingerProblem`](@ref);
`dQ = z -> z^k` is ``∂/∂v_k``. Unlike [`wkb_derivative`](@ref) this needs no Riccati
algebra, so it also covers a [`RationalProblem`](@ref), where `λ` is whichever
parameter the caller differentiates `Q` by.
"""
function period_derivative(prob::AbstractSchrodingerProblem, contour::AbstractVector,
                           dQ; closed::Bool = true, rtol = nothing, atol = nothing,
                           maxevals::Integer = 10^6)
    closed || throw(ContourError(
        "a parameter derivative is taken at frozen contour, which is exact only for a " *
        "closed cycle; on an open contour the turning-point endpoints move with the " *
        "parameter and the dropped boundary term is not small"))
    F = _contour_float(_as_complex.(contour))
    verts = Complex{F}[_as_complex(z) for z in contour]
    length(verts) >= 2 || throw(ContourError("contour needs at least two vertices"))
    all(isfinite, verts) || throw(ContourError("contour vertices must be finite"))
    closed && push!(verts, verts[1])
    K = length(verts) - 1
    K ≥ 1 || throw(ContourError("contour needs ≥ 2 vertices"))
    geom_tol = _period_atol(prob, verts, F)
    S, PA, UW, Δ = _branch_table(prob, verts, K, geom_tol)
    if closed && isodd(round(Int, Δ / (2π)))
        throw(ContourError("closed contour encloses an odd number of turning points; " *
                           "√Q is not single-valued around it (branch ambiguity)"))
    end
    function f(s)
        u, z, _, dz = _sqrt_Q(prob, s, S, PA, UW, verts, K)
        dQ(z) / (2 * u) * dz
    end
    rt = rtol === nothing ? sqrt(eps(F)) : F(rtol)
    at = atol === nothing ? zero(F) : F(atol)
    isfinite(rt) && isfinite(at) && rt >= 0 && at >= 0 && maxevals > 0 ||
        throw(Resurgence.InvalidArgument("quadrature tolerances must be finite and nonnegative; maxevals must be positive"))
    val, err = quadgk(f, S...;
                    rtol = rt, atol = at, maxevals = maxevals, order = _quad_order(F))
    err <= max(rt * abs(val), at) || throw(ContourError("period quadrature did not reach its requested tolerance"))
    val
end

"""
    period_derivative(prob::SchrodingerProblem, contour; wrt = :energy, index = 0)
        -> Complex

``∂/∂E`` (`wrt = :energy`) or ``∂/∂v_k`` (`wrt = :coefficient`, `index = k`) of the
classical period, i.e. the two `dQ` a polynomial potential supplies itself.
"""
function period_derivative(prob::SchrodingerProblem, contour::AbstractVector;
                           wrt::Symbol = :energy, index::Integer = 0, kwargs...)
    dQ = if wrt === :energy
        _ -> -1
    elseif wrt === :coefficient
        dmax = length(q_coefficients(prob)) - 1
        0 ≤ index ≤ dmax || throw(Resurgence.InvalidArgument(
            "index must be in 0:$dmax (the powers of z carried by V), got $index"))
        z -> z^index
    else
        throw(Resurgence.InvalidArgument(
            "wrt must be :energy or :coefficient, got :$wrt"))
    end
    period_derivative(prob, contour, dQ; kwargs...)
end

"""
    encircling_contour(tp1, tp2; margin = nothing, n = 64) -> Vector{Complex}

A closed elliptical `n`-gon enclosing the two turning points `tp1`, `tp2` (either
`TurningPoint`s or complex numbers), oriented along their separation. Suitable as the
`contour` argument to [`period_integral`](@ref) / [`wkb_period`](@ref) for the quantum
period of the cycle around that pair. `margin` sets the padding (default `≈` half the
separation).
"""
function encircling_contour(tp1, tp2; margin = nothing, n::Integer = 32)
    z1, z2 = _as_complex(tp1), _as_complex(tp2)
    F = _contour_float(typeof(z1 + z2 + 0.0im))
    c = (z1 + z2) / 2
    d = z2 - z1
    a = abs(d) / 2
    dhat = a ≈ 0 ? complex(one(F)) : d / abs(d)
    perp = im * dhat
    pad = margin === nothing ? (a > 0 ? a / 2 : one(F)) : F(margin)
    rmaj = a + pad
    rmin = pad
    [Complex{F}(c + rmaj * cos(2π * j / n) * dhat + rmin * sin(2π * j / n) * perp)
     for j in 0:(n - 1)]
end
