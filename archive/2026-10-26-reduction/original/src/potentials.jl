# The one-dimensional Schrödinger problem ħ²ψ″ = Q(z)ψ, Q = V − E (Iwaki–Nakanishi
# normalization). The potential is a polynomial in z given by its coefficients in
# ascending order; the energy is a separate scalar so `with_energy` is cheap and the
# family {Q_E} is one object. Exact coefficient types (Rational, Int) are first-class.

# Horner evaluation of a coefficient vector (ascending) at `z`, promoting to the join
# of the coefficient and argument types. Mirrors Resurgence's FormalSeries.evaluate.
function _horner(coeffs, z)
    foldr((c, acc) -> c + z * acc, coeffs; init = zero(first(coeffs)) * zero(z))
end

"""
    AbstractSchrodingerProblem

Supertype of the one-dimensional problems ``ħ²ψ″ = Q(z)ψ`` the package traces:
[`SchrodingerProblem`](@ref) (polynomial `Q`) and [`RationalProblem`](@ref)
(rational `Q`, i.e. `Q` with poles). Everything the Stokes-geometry layer needs is
reached through four methods - `prob(z)`, [`q_derivative_at`](@ref),
[`q_taylor_at`](@ref) and [`turning_points`](@ref) - so a new potential class is one
file, not a rewrite.
"""
abstract type AbstractSchrodingerProblem end

# Taylor shift of a coefficient vector (ascending) to the point `z0`, by repeated
# synthetic division. Exact whenever the coefficients and `z0` are.
function _taylor_shift(q::AbstractVector, z0)
    T = typeof(zero(eltype(q)) * zero(z0) + zero(eltype(q)))
    a = T[T(c) for c in q]
    out = T[]
    while true
        n = length(a)
        if n == 1
            push!(out, a[1])
            break
        end
        b = Vector{T}(undef, n - 1)
        b[n - 1] = a[n]
        for i in (n - 1):-1:2
            b[i - 1] = a[i] + z0 * b[i]
        end
        push!(out, a[1] + z0 * b[1])
        a = b
    end
    out
end

# Derivative of a coefficient vector (ascending).
function _poly_derivative(q::AbstractVector)
    length(q) < 2 && return [zero(eltype(q))]
    [k * q[k + 1] for k in 1:(length(q) - 1)]
end

"""
    SchrodingerProblem{T}

A one-dimensional Schrödinger problem ``ħ²ψ″ = Q(z)ψ`` with ``Q(z) = V(z) − E``.
The potential `V` is a polynomial in the variable `var` (default `:z`) stored by its
coefficients `v_coeffs` in **ascending** order (`v_coeffs[k]` multiplies ``z^{k-1}``);
`energy` is ``E``. Immutable: [`with_energy`](@ref) returns a new problem.

The coefficient type `T` is generic - exact `Rational`/`Int` and real floats are all
first-class. Construction throws [`InvalidPotential`](@ref) on an empty coefficient
vector or on a `Q` that is identically constant (no turning points).

# Examples
```julia
airy   = SchrodingerProblem([0, 1])                       # Q = z
harm   = SchrodingerProblem([0, 0, 1]; energy = 1)        # Q = z² − 1
dwell  = SchrodingerProblem([3//4, 0, -2, 0, 1])          # Q = (z²−1)² − 1/4
airy(2.0)                                                 # Q(2) = 2.0
```
"""
struct SchrodingerProblem{T} <: AbstractSchrodingerProblem
    v_coeffs::Vector{T}
    energy::T
    var::Symbol

    function SchrodingerProblem{T}(v_coeffs::Vector{T}, energy::T,
                                   var::Symbol) where {T}
        isempty(v_coeffs) &&
            throw(InvalidPotential("a SchrodingerProblem needs ≥ 1 potential coefficient"))
        q = copy(v_coeffs)
        q[1] -= energy
        last_nz = findlast(!iszero, q)
        (last_nz === nothing || last_nz < 2) &&
            throw(InvalidPotential("Q = V − E is constant (degree < 1); no turning " *
                                   "points and no WKB problem"))
        new{T}(copy(v_coeffs), energy, var)
    end
end

function SchrodingerProblem(v_coeffs::AbstractVector; energy = 0, var::Symbol = :z)
    T = promote_type(eltype(v_coeffs), typeof(energy))
    SchrodingerProblem{T}(collect(T, v_coeffs), convert(T, energy), var)
end

# -- accessors ---------------------------------------------------------------------

"""
    q_coefficients(prob::SchrodingerProblem) -> Vector

Coefficients of ``Q = V − E`` in ascending order, with trailing zeros trimmed so the
length is `degree(prob) + 1`.
"""
function q_coefficients(prob::SchrodingerProblem)
    q = copy(prob.v_coeffs)
    q[1] -= prob.energy
    q[1:findlast(!iszero, q)]
end

"""
    energy(prob::SchrodingerProblem) -> Number

The energy ``E``.
"""
energy(prob::SchrodingerProblem) = prob.energy

"""
    degree(prob::SchrodingerProblem) -> Int

Degree of the polynomial ``Q = V − E`` (≥ 1 by construction).
"""
degree(prob::SchrodingerProblem) = length(q_coefficients(prob)) - 1

"""
    variable(prob::SchrodingerProblem) -> Symbol

The coordinate symbol (default `:z`).
"""
variable(prob::SchrodingerProblem) = prob.var

"""
    with_energy(prob::SchrodingerProblem, E) -> SchrodingerProblem

A new problem with the same potential and energy `E`.
"""
function with_energy(prob::SchrodingerProblem, E)
    SchrodingerProblem(prob.v_coeffs; energy = E, var = prob.var)
end

# -- evaluation --------------------------------------------------------------------

"""
    (prob::SchrodingerProblem)(z) -> Number

Evaluate ``Q(z) = V(z) − E`` by Horner, promoting to `typeof(z)`.
"""
(prob::SchrodingerProblem)(z) = _horner(q_coefficients(prob), z)

"""
    q_derivative_at(prob::SchrodingerProblem, z) -> Number

Evaluate ``Q′(z)`` by Horner on the differentiated coefficient vector.
"""
function q_derivative_at(prob::SchrodingerProblem, z)
    q = q_coefficients(prob)
    length(q) < 2 && return zero(first(q)) * zero(z)
    _horner(_poly_derivative(q), z)
end

"""
    q_taylor_at(prob::SchrodingerProblem, z0) -> Vector

The Taylor coefficients of ``Q`` about `z0` in **ascending** order: the returned
`t` satisfies ``Q(z) = Σ_k t[k+1] (z − z_0)^k``, so `t[1] = Q(z₀)`,
`t[2] = Q′(z₀)` and `t[m+1] = Q^{(m)}(z₀)/m!`.

Computed by repeated synthetic division, so the result is **exact** whenever the
coefficient type and `z0` are (`Rational`, `Integer`). This is how the leading
coefficient ``c_m`` of a turning point of order `m` is read off - the datum that
sets its Stokes-ray directions (see [`stokes_graph`](@ref)).
"""
q_taylor_at(prob::SchrodingerProblem, z0) = _taylor_shift(q_coefficients(prob), z0)

# -- the generic interface the Stokes-geometry layer reads --------------------------
#
# `_turning_polynomial` - the polynomial whose roots are the turning points.
# `_infinity_exponent` / `_infinity_coefficient` - `e` and `a` in `Q ~ a·z^e` as
# `z → ∞`. Together they fix the asymptotic Stokes directions and the escape mass.
# `n_finite_poles` is 0 for a polynomial problem; see `src/rational_potentials.jl`.

_turning_polynomial(prob::SchrodingerProblem) = q_coefficients(prob)
_infinity_exponent(prob::SchrodingerProblem) = degree(prob)
_infinity_coefficient(prob::SchrodingerProblem) = last(q_coefficients(prob))
n_finite_poles(prob::AbstractSchrodingerProblem) = 0

"""
    asymptotic_directions(prob, theta; pole = 0) -> Vector

The directions in which Stokes lines run into a singularity of `Q` at phase `theta`.
Locally `Q ~ c·ζ^e` (`ζ = z − p` at a finite pole `p`, `ζ = 1/z` inverted at
infinity), and `Im[e^{−iθ}∫√Q] = 0` puts `|e + 2|` directions at

    arg ζ = (2θ − arg c + 2πk)/(e + 2) ,   k = 0 … |e+2| − 1 .

`pole = 0` selects infinity (`e = degree(prob)`, so a degree-`d` polynomial gets the
familiar `d + 2` asymptotic directions); `pole = j ≥ 1` selects the `j`-th finite pole
of a [`RationalProblem`](@ref), where `e = −m_j` and the count is `m_j − 2`. The
returned angles are the *marked points* of the boundary circle at that singularity.

A **double** pole gives `e + 2 = 0` and the empty vector: it is a puncture, carrying
no marked points, and the trajectories spiral into it instead of approaching a
direction (see [`ring_domain_walls`](@ref)).
"""
function asymptotic_directions(prob::AbstractSchrodingerProblem, theta; pole::Integer = 0)
    e, c = pole == 0 ? (_infinity_exponent(prob), _infinity_coefficient(prob)) :
                       (_pole_exponent(prob, pole), _pole_coefficient(prob, pole))
    F = typeof(float(real(zero(c) + zero(theta))))
    n = abs(e + 2)
    # A double pole (e = −2) is a puncture: it carries no marked points at all, and
    # the Stokes lines spiral in rather than approaching any direction. Returning
    # nothing is the right answer, not an error - it is what makes the marked-point
    # count `Σ|e+2|` one formula across irregular poles and punctures.
    n ≥ 1 || return F[]
    [mod((2F(theta) - angle(Complex{F}(c)) + 2 * F(π) * k) / (e + 2), 2 * F(π))
     for k in 0:(n - 1)]
end

Base.:(==)(p::SchrodingerProblem, q::SchrodingerProblem) =
    p.var == q.var && p.energy == q.energy && p.v_coeffs == q.v_coeffs
