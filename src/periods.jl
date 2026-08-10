# Quantum periods: contour integrals ∮ S_m dz along piecewise-linear paths, with the
# √Q branch tracked by unwrapping arg Q(z) along the contour. The classical period
# ∮√Q dz is the m = -1 case. Branch bookkeeping is the top correctness risk, so it is
# isolated here: an adaptive pre-walk lays down nodes dense enough that arg Q never
# jumps by more than ~π/2 between neighbours, and every quadrature sample unwraps
# against the nearest node - no discrete sign choices inside the integrand.

using QuadGK: quadgk

# principal-branch-safe wrap of an angle into (−π, π] (type-generic: uses BigFloat π
# for BigFloat arguments, so branch unwrapping keeps full precision)
function _wrap(θ)
    τ = 2 * oftype(θ, π)
    θ - τ * round(θ / τ)
end

# the float type carried by a contour / turning-point list
_contour_float(::Type{Complex{F}}) where {F} = F
_contour_float(::Type{<:Real}) = Float64
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

# Adaptive connection: append nodes from (s0,z0) to (s1,z1), subdividing until the
# principal arg-Q step is ≤ π/2, carrying the unwrapped angle. Returns (pa1, uw1).
function _connect!(prob, s0, z0, pa0, uw0, s1, z1, S, PA, UW, atol, depth)
    Q1 = prob(z1)
    abs(Q1) ≤ atol &&
        throw(ContourError("contour passes through (or grazes) a turning point at " *
                           "z ≈ $z1 where |Q| ≈ $(abs(Q1))"))
    pa1 = angle(Q1)
    d = _wrap(pa1 - pa0)
    if abs(d) > π / 2 && depth > 0
        sm = (s0 + s1) / 2
        zm, _ = _contour_point_from_ends(z0, z1)
        pam, uwm = _connect!(prob, s0, z0, pa0, uw0, sm, zm, S, PA, UW, atol, depth - 1)
        return _connect!(prob, sm, zm, pam, uwm, s1, z1, S, PA, UW, atol, depth - 1)
    end
    uw1 = uw0 + d
    push!(S, s1); push!(PA, pa1); push!(UW, uw1)
    pa1, uw1
end

_contour_point_from_ends(z0, z1) = ((z0 + z1) / 2, z1 - z0)

# Build the branch-tracking node table (sorted global s, principal args, unwrapped
# args) over the contour. Also returns the net unwrapped Δ over the whole path.
function _branch_table(prob, verts::Vector{Complex{F}}, K, atol; maxdepth = 40) where {F}
    Q0 = prob(verts[1])
    abs(Q0) ≤ atol &&
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

# atol for turning-point / grazing detection, scaled to the potential's size on the path.
#
# `sqrt(eps)` is a RELATIVE threshold, so the scale it multiplies has to be a *typical*
# `|Q|` along the path, not an extremal one. The maximum is not typical once `Q` has
# poles: a contour that runs anywhere near one takes `max|Q|` to infinity, and the
# tolerance then swamps perfectly good values. Measured case (an order-4 pole 0.017 away
# from the contour): `max|Q| = 6.0e6` put `atol` at 0.09 while the honest `|Q|` at the
# contour's first vertex was 0.015, so a legitimate vertex was rejected as a turning
# point and `charge_basis(prob, g)` failed on two of three chambers.
#
# The median is unmoved by those few huge samples, still goes to zero only where `Q`
# genuinely does, and on a polynomial contour differs from the old scale by a constant -
# so this loosens nothing that used to be caught (a vertex *on* a turning point has
# `|Q| = 0` and trips any positive tolerance).
function _period_atol(prob, verts, F)
    isempty(verts) && return sqrt(eps(F))
    mags = sort!([F(abs(prob(z))) for z in verts])
    sqrt(eps(F)) * (1 + mags[cld(length(mags), 2)])
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
float type.
"""
function period_integral(prob::AbstractSchrodingerProblem, contour::AbstractVector;
                         closed::Bool = true, rtol = nothing, maxevals::Integer = 10^6)
    F = _contour_float(_as_complex.(contour))
    verts = Complex{F}[_as_complex(z) for z in contour]
    closed && push!(verts, verts[1])
    K = length(verts) - 1
    K ≥ 1 || throw(ContourError("contour needs ≥ 2 vertices"))
    atol = _period_atol(prob, verts, F)
    S, PA, UW, Δ = _branch_table(prob, verts, K, atol)
    if closed && isodd(round(Int, Δ / (2π)))
        throw(ContourError("closed contour encloses an odd number of turning points; " *
                           "√Q is not single-valued around it (branch ambiguity)"))
    end
    function f(s)
        u, _, _, dz = _sqrt_Q(prob, s, S, PA, UW, verts, K)
        u * dz
    end
    rt = rtol === nothing ? sqrt(eps(F)) : F(rtol)
    val, _ = quadgk(f, range(zero(F), one(F), length = K + 1)...;
                    rtol = rt, maxevals = maxevals, order = _quad_order(F))
    val
end

"""
    wkb_period(w::WKBExpansion, contour, m; closed = true) -> Complex

The order-`m` quantum period ``∮ S_m\\,dz`` along `contour`, using the same √Q branch
tracking as [`period_integral`](@ref). `m = -1` reproduces the classical period.
"""
function wkb_period(w::WKBExpansion, contour::AbstractVector, m::Integer;
                    closed::Bool = true, rtol = nothing, atol = nothing,
                    maxevals::Integer = 10^6)
    -1 ≤ m ≤ w.order ||
        throw(Resurgence.InvalidArgument("m must be in -1:$(w.order), got $m"))
    F = _contour_float(_as_complex.(contour))
    verts = Complex{F}[_as_complex(z) for z in contour]
    closed && push!(verts, verts[1])
    K = length(verts) - 1
    atol = _period_atol(w.prob, verts, F)
    S, PA, UW, Δ = _branch_table(w.prob, verts, K, atol)
    if closed && isodd(round(Int, Δ / (2π)))
        throw(ContourError("closed contour encloses an odd number of turning points; " *
                           "√Q is not single-valued around it (branch ambiguity)"))
    end
    term = _s_term(w, m)
    function f(s)
        u, z, Q, dz = _sqrt_Q(w.prob, s, S, PA, UW, verts, K)
        _eval_term(term, Complex{F}(z), Q, u) * dz
    end
    rt = rtol === nothing ? sqrt(eps(F)) : F(rtol)
    at = atol === nothing ? zero(F) : F(atol)
    val, _ = quadgk(f, range(zero(F), one(F), length = K + 1)...;
                    rtol = rt, atol = at, maxevals = maxevals,
                    order = _quad_order(F))
    val
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
