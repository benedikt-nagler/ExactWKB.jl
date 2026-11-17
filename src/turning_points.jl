# Turning points: the zeros of Q(z). Found numerically with PolynomialRoots, polished
# by Newton against the exact Q, then clustered so a multiple root is reported once
# with its multiplicity (its order). Exact-rational potentials are lifted to BigFloat
# at the caller's `setprecision`; float potentials stay in their own precision.

import PolynomialRoots

# The working float type for numerics derived from a coefficient element type:
# exact coefficients → BigFloat (caller's precision), floats keep their own type.
_wkb_float(::Type{T}) where {T<:AbstractFloat} = T
_wkb_float(::Type{Complex{T}}) where {T} = _wkb_float(T)
_wkb_float(::Type{<:Union{Integer,Rational}}) = BigFloat
_wkb_float(prob::AbstractSchrodingerProblem) = _wkb_float(eltype(_turning_polynomial(prob)))

"""
    TurningPoint{F}

A turning point of a [`SchrodingerProblem`](@ref): a zero `z::Complex{F}` of ``Q``
with multiplicity `order` (`order == 1` is simple). `F` is the working float type.
"""
struct TurningPoint{F}
    z::Complex{F}
    order::Int
end

"""
    location(tp::TurningPoint) -> Complex

The turning point's position ``z``.
"""
location(tp::TurningPoint) = tp.z

"""
    order(tp::TurningPoint) -> Int

The multiplicity of the turning point as a zero of ``Q``.
"""
order(tp::TurningPoint) = tp.order

"""
    is_simple(tp::TurningPoint) -> Bool

`true` for a simple (order-1) turning point.
"""
is_simple(tp::TurningPoint) = tp.order == 1

Base.:(==)(s::TurningPoint, t::TurningPoint) = s.z == t.z && s.order == t.order

# Exact Euclidean arithmetic on the polynomial represented by the input numbers.
# Floating coefficients are converted to their exact binary rational values.
# Multiplicity is algebraic; distances between approximate roots never define it.
_trim_poly(a) = (while length(a) > 1 && iszero(last(a)); pop!(a); end; a)
function _root_divrem(a, b)
    a = copy(a)
    q = zeros(eltype(a), max(1, length(a) - length(b) + 1))
    while length(a) >= length(b) && !all(iszero, a)
        k = length(a) - length(b)
        c = last(a) / last(b)
        q[k + 1] = c
        for j in eachindex(b)
            a[k + j] -= c * b[j]
        end
        _trim_poly(a)
    end
    _trim_poly(q), a
end
function _root_gcd(a, b)
    while !all(iszero, b)
        _, rem = _root_divrem(a, b)
        a, b = b, rem
    end
    a ./ last(a)
end
_root_quotient(a, b) = first(_root_divrem(a, b))
function _squarefree_parts(q)
    T = any(c -> !isreal(c), q) ? Complex{Rational{BigInt}} : Rational{BigInt}
    f = T[isreal(c) ? Rational{BigInt}(real(c)) :
          complex(Rational{BigInt}(real(c)), Rational{BigInt}(imag(c))) for c in q]
    f ./= last(f)
    c = _root_gcd(f, _poly_derivative(f))
    w = _root_quotient(f, c)
    out = Tuple{Vector{T},Int}[]
    m = 1
    while length(w) > 1
        y = _root_gcd(w, c)
        z = _root_quotient(w, y)
        length(z) > 1 && push!(out, (z, m))
        w, c = y, _root_quotient(c, y)
        m += 1
    end
    out
end

"""
    turning_points(prob::AbstractSchrodingerProblem; cluster_tol = nothing)

Locate zeros of the polynomial numerator and report their algebraic multiplicities.
Multiplicity refers to the exact polynomial represented by the supplied coefficients.
Floating coefficients are treated as exact binary rationals for this purpose.
Root positions are numerical and follow the coefficient precision.
`cluster_tol` is a separation guard, never permission to merge roots.
Throws `Resurgence.InvalidArgument` when distinct roots are numerically unresolved.
"""
function turning_points(prob::AbstractSchrodingerProblem; cluster_tol = nothing)
    F = _wkb_float(prob)
    q = _turning_polynomial(prob)
    all(isfinite, q) || throw(InvalidPotential("potential coefficients must be finite"))
    tol = cluster_tol === nothing ? 8 * eps(F) : F(cluster_tol)
    isfinite(tol) && tol >= 0 || throw(Resurgence.InvalidArgument(
        "cluster_tol must be finite and nonnegative"))
    tps = TurningPoint{F}[]
    for (factor, multiplicity) in _squarefree_parts(q)
        c = Complex{F}.(factor)
        dc = _poly_derivative(c)
        roots = length(c) == 2 ? [-c[1] / c[2]] : PolynomialRoots.roots(c)
        for z0 in roots
            z = z0
            for _ in 1:50
                dz = _horner(dc, z)
                iszero(dz) && break
                step = _horner(c, z) / dz
                z -= step
                abs(step) <= eps(F) * (1 + abs(z)) && break
            end
            isfinite(z) || throw(Resurgence.InvalidArgument("root location is not finite"))
            scale = sum(abs(c[k])*abs(z)^(k-1) for k in eachindex(c))
            abs(_horner(c,z)) <= 100*eps(F)*scale || throw(Resurgence.InvalidArgument(
                "root residual is unresolved at the working precision"))
            push!(tps, TurningPoint{F}(z, multiplicity))
        end
    end
    for i in eachindex(tps), j in (i + 1):length(tps)
        a, b = location(tps[i]), location(tps[j])
        abs(a - b) > tol * (1 + max(abs(a), abs(b))) ||
            throw(Resurgence.InvalidArgument(
                "distinct turning points are unresolved; increase coefficient precision"))
    end
    sort!(tps; by = t -> (real(t.z), imag(t.z)))
end

"""
    turning_points(prob; ...) restricted with `simple = true`

The simple turning points only (order 1).
"""
simple_turning_points(prob::AbstractSchrodingerProblem; kwargs...) =
    filter(is_simple, turning_points(prob; kwargs...))
