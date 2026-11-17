# Rational Q = N(z)/prod((z-p)^m), with supplied pole locations and orders.
# Ledger: finite orders >= 2; degree(Q) >= -1 at infinity. These are input limits.
# Only classical periods and numerical trajectories are supported here.

"""
    RationalProblem(n_coeffs, poles, pole_orders; var = :z)

Represent `Q = N(z)/prod((z-p)^m)` with ascending numerator coefficients.
Pole locations and orders are supplied data. Finite orders must be at least two.
The exponent of Q at infinity must be at least -1.
Supports classical periods and numerical tracing, not the Riccati recursion.
Throws [`InvalidPotential`](@ref) for unsupported or cancelling pole data.
"""
struct RationalProblem{T} <: AbstractSchrodingerProblem
    n_coeffs::Vector{T}
    poles::Vector{T}
    pole_orders::Vector{Int}
    var::Symbol

    function RationalProblem{T}(n_coeffs::Vector{T}, poles::Vector{T},
                                pole_orders::Vector{Int}, var::Symbol) where {T}
        isempty(n_coeffs) &&
            throw(InvalidPotential("a RationalProblem needs ≥ 1 numerator coefficient"))
        last_nz = findlast(!iszero, n_coeffs)
        last_nz === nothing &&
            throw(InvalidPotential("the numerator of Q is identically zero"))
        n = n_coeffs[1:last_nz]
        length(poles) == length(pole_orders) || throw(InvalidPotential(
            "$(length(poles)) poles but $(length(pole_orders)) pole orders"))
        allunique(poles) || throw(InvalidPotential("the poles of Q must be distinct"))
        for (j, m) in enumerate(pole_orders)
            m ≥ 2 || throw(InvalidPotential(
                "pole $j (at $(poles[j])) has order $m: a simple pole is an " *
                "orbifold point ([IN16]), which needs orbifold points on the " *
                "surface side and a Bessel local model on the WKB side. Poles of " *
                "order 2 (a puncture) and ≥ 3 (irregular) are supported"))
        end
        # the numerator must not vanish at a pole, or the true order is lower
        scale = maximum(abs, n; init = one(real(T)))
        for (j, p) in enumerate(poles)
            abs(_horner(n, p)) > sqrt(eps(float(real(one(T))))) * scale ||
                throw(InvalidPotential(
                    "the numerator vanishes at pole $j (z = $p), so its order is " *
                    "not the declared $(pole_orders[j])"))
        end
        e = (length(n) - 1) - sum(pole_orders; init = 0)
        e ≥ -1 || throw(InvalidPotential(
            "Q ~ a·z^$e at infinity, so Q dz² has a pole of order $(e + 4) there: " *
            "infinity is a regular singularity, and only *finite* poles may be " *
            "regular. Move it with a Möbius change of variable"))
        new{T}(n, copy(poles), copy(pole_orders), var)
    end
end

function RationalProblem(n_coeffs::AbstractVector, poles::AbstractVector,
                         pole_orders::AbstractVector{<:Integer}; var::Symbol = :z)
    T = isempty(poles) ? eltype(n_coeffs) :
                         promote_type(eltype(n_coeffs), eltype(poles))
    RationalProblem{T}(collect(T, n_coeffs), collect(T, poles),
                       collect(Int, pole_orders), var)
end

# -- accessors ---------------------------------------------------------------------

"""
    q_numerator(prob::RationalProblem) -> Vector

Coefficients of the numerator `N` of `Q = N/∏(z − p_j)^{m_j}`, ascending, with
trailing zeros trimmed.
"""
q_numerator(prob::RationalProblem) = copy(prob.n_coeffs)

"""
    poles(prob) -> Vector

The finite poles of `Q` (empty for a [`SchrodingerProblem`](@ref)).
"""
poles(prob::RationalProblem) = copy(prob.poles)
poles(::SchrodingerProblem{T}) where {T} = T[]

"""
    pole_orders(prob) -> Vector{Int}

The orders of the finite poles of `Q`, in the order of [`poles`](@ref).
"""
pole_orders(prob::RationalProblem) = copy(prob.pole_orders)
pole_orders(::AbstractSchrodingerProblem) = Int[]

"""
    n_finite_poles(prob) -> Int

The number of finite poles of `Q`, counted **without** multiplicity - the length of
[`poles`](@ref), not the sum of [`pole_orders`](@ref). Zero for a problem with no finite
poles, so this is the cheap test for whether a problem is polynomial in disguise.
"""
n_finite_poles(prob::RationalProblem) = length(prob.poles)

"""
    puncture_indices(prob) -> Vector{Int}

The indices (into [`poles`](@ref)) of the **double** poles of `Q` - the punctures of
the surface. An irregular pole (order ≥ 3) carries a boundary circle instead; see
[`RationalProblem`](@ref).
"""
puncture_indices(prob::RationalProblem) = findall(==(2), prob.pole_orders)
puncture_indices(::AbstractSchrodingerProblem) = Int[]

"""
    n_punctures(prob) -> Int

The number of punctures (double poles) of `Q`. Zero for a polynomial problem and for
a purely irregular rational one.
"""
n_punctures(prob::AbstractSchrodingerProblem) = length(puncture_indices(prob))

variable(prob::RationalProblem) = prob.var

"""
    degree(prob::RationalProblem) -> Int

The exponent of `Q` at infinity: `Q ~ a·z^{degree}`, i.e. the numerator degree minus
the total pole order. May be negative. This is the same quantity `degree` returns for
a polynomial problem, so [`asymptotic_directions`](@ref) is one formula for both; the
pole order of `Q dz²` at infinity is `degree + 4`.
"""
degree(prob::RationalProblem) =
    (length(prob.n_coeffs) - 1) - sum(prob.pole_orders; init = 0)

_turning_polynomial(prob::RationalProblem) = prob.n_coeffs
_infinity_exponent(prob::RationalProblem) = degree(prob)
_infinity_coefficient(prob::RationalProblem) = last(prob.n_coeffs)

# `e` and `c` in `Q ~ c·(z − p_j)^e` at the j-th finite pole.
_pole_exponent(prob::RationalProblem, j::Integer) = -prob.pole_orders[j]

function _pole_coefficient(prob::RationalProblem, j::Integer)
    p = prob.poles[j]
    c = _horner(prob.n_coeffs, p)
    for i in eachindex(prob.poles)
        i == j && continue
        c /= (p - prob.poles[i])^prob.pole_orders[i]
    end
    c
end

_pole_exponent(::AbstractSchrodingerProblem, ::Integer) =
    throw(InvalidPotential("this problem has no finite poles"))
_pole_coefficient(::AbstractSchrodingerProblem, ::Integer) =
    throw(InvalidPotential("this problem has no finite poles"))

# -- evaluation --------------------------------------------------------------------

# the denominator ∏(z − p_j)^{m_j}
function _denominator(prob::RationalProblem, z)
    d = one(zero(eltype(prob.poles)) * zero(z) + one(eltype(prob.poles)))
    for (p, m) in zip(prob.poles, prob.pole_orders)
        d *= (z - p)^m
    end
    d
end

"""
    (prob::RationalProblem)(z) -> Number

Evaluate ``Q(z) = N(z)/∏(z − p_j)^{m_j}``.
"""
(prob::RationalProblem)(z) = _horner(prob.n_coeffs, z) / _denominator(prob, z)

"""
    q_derivative_at(prob::RationalProblem, z) -> Number

Evaluate ``Q′(z)``. Computed as `(N′ − N·Σ m_j/(z − p_j))/D` rather than by the
quotient rule on `D²`, so it is stable at a turning point (where `N` vanishes) and
needs no cancellation.
"""
function q_derivative_at(prob::RationalProblem, z)
    n = _horner(prob.n_coeffs, z)
    dn = _horner(_poly_derivative(prob.n_coeffs), z)
    s = zero(n) * zero(z)
    for (p, m) in zip(prob.poles, prob.pole_orders)
        s += m / (z - p)
    end
    (dn - n * s) / _denominator(prob, z)
end

# Truncated series division a/b (both ascending, b[1] ≠ 0), keeping `len` terms.
function _series_divide(a::AbstractVector, b::AbstractVector, len::Integer)
    T = typeof(zero(eltype(a)) / one(eltype(b)))
    c = zeros(T, len)
    for k in 1:len
        acc = k ≤ length(a) ? T(a[k]) : zero(T)
        for i in 1:(k - 1)
            # the guard is on the index INTO b, not on i: with a double pole the
            # numerator outruns the denominator, which no irregular-pole problem does
            k - i + 1 ≤ length(b) && (acc -= c[i] * b[k - i + 1])
        end
        c[k] = acc / b[1]
    end
    c
end

"""
    q_taylor_at(prob::RationalProblem, z0) -> Vector

The Taylor coefficients of ``Q`` about `z0` in **ascending** order (`z0` must not be
a pole). Both numerator and denominator are Taylor-shifted exactly and then divided
as truncated series, so a turning point's leading coefficient `c_m = Q^{(m)}(z0)/m!`
- the datum that sets its Stokes-ray directions - is as exact as the inputs.
"""
function q_taylor_at(prob::RationalProblem, z0)
    any(p -> p == z0, prob.poles) &&
        throw(InvalidPotential("q_taylor_at: z0 = $z0 is a pole of Q"))
    num = _taylor_shift(prob.n_coeffs, z0)
    den = [one(eltype(num))]
    for (p, m) in zip(prob.poles, prob.pole_orders)
        f = _taylor_shift([-p, one(p)], z0)          # (z − p) about z0
        for _ in 1:m
            den = _poly_multiply(den, f)
        end
    end
    _series_divide(num, den, length(num))
end

function _poly_multiply(a::AbstractVector, b::AbstractVector)
    T = typeof(zero(eltype(a)) * zero(eltype(b)))
    c = zeros(T, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i + j - 1] += a[i] * b[j]
    end
    c
end

Base.:(==)(p::RationalProblem, q::RationalProblem) =
    p.var == q.var && p.n_coeffs == q.n_coeffs && p.poles == q.poles &&
    p.pole_orders == q.pole_orders
