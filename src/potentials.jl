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
struct SchrodingerProblem{T}
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
    dq = [k * q[k + 1] for k in 1:(length(q) - 1)]
    _horner(dq, z)
end

Base.:(==)(p::SchrodingerProblem, q::SchrodingerProblem) =
    p.var == q.var && p.energy == q.energy && p.v_coeffs == q.v_coeffs
