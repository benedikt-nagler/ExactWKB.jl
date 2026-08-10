# Rational potentials: ħ²ψ″ = Q(z)ψ with Q = N(z)/∏_j (z − p_j)^{m_j}.
#
# Why this file exists. A polynomial Q draws its Stokes graph on a *disk* with d+2
# marked points, so the cluster types reachable through the bridge are exactly
# A_{d−1} and no more (see the M4 scope boundary in PLANNING/roadmap.md). Poles are
# the way out: in the Gaiotto/GMN dictionary the quadratic differential Q dz² on ℙ¹
# has, at each singularity, a boundary circle carrying (pole order − 2) marked
# points, so k+1 irregular singularities give a sphere with k+1 holes, and the dual
# object is an ideal triangulation of *that* surface. Two holes = the annulus =
# Ã(p,q).
#
# ── Convention ledger ──────────────────────────────────────────────────────────────
# Sources: [GMN-TBA] §6 (the quadratic-differential ↔ marked-surface dictionary: a
# pole of order p ≥ 3 is an irregular singularity carrying p − 2 marked points on a
# boundary circle, a pole of order 2 a puncture), [IN16] (poles in the exact-WKB/
# cluster dictionary), [FST08] (tagged arcs, which the refused self-folded
# configurations would need); keys resolve in PLANNING/references.md.
#
#  • POLES ARE EXACT DATA, NOT ROOTS. A pole is given by its location and order, not
#    by a denominator polynomial to be root-found. Pole orders set the marked-point
#    count, which is a *topological* datum - deriving it from a numerical root
#    multiplicity would put the surface's homeomorphism type at the mercy of a
#    clustering tolerance. Pinned by: the whole triangulation layer asserts
#    n = M + 2b − 4 exactly (`ideal_triangulation`).
#
#  • A FINITE POLE HAS ORDER ≥ 2; ORDER 1 IS REFUSED. Order 1 is an orbifold point
#    ([IN16]) and needs the surface side to grow orbifold points, which
#    ClusterSurfaces.jl does not have. Order 2 is a *regular puncture* and IS
#    supported: it is not a boundary circle with zero marked points but a
#    **vertex** of the dual triangulation, which is the whole content of the
#    puncture layer (see the ledger in `src/triangulation.jl`).
#
#  • INFINITY STAYS IRREGULAR: `Q ~ a z^e` makes the pole order of Q dz² there equal
#    e + 4, so `e ≥ −1` is required and `e = −2` (a puncture at infinity) is refused.
#    This is a scope boundary, not a mathematical one - the tracer's escape circle is
#    the one place that would have to learn to spiral outward - and a Möbius map moves
#    any such puncture to a finite point.
#
#  • A DOUBLE POLE CARRIES A WALL OF ITS OWN. Near it `√Q ~ √c/(z−p)`, so
#    `∮√Q dz = 2πi√c =: Z_p` and the trajectories are logarithmic spirals into the
#    puncture *except* at `θ ≡ arg Z_p (mod π)`, where they close up into a ring
#    domain. That is a wall in exactly the sense a saddle connection is, and it is
#    what [`ring_domain_walls`](@ref) reports. Pinned by: `stokes_graph` on the
#    once-punctured triangle changes the puncture's valence across it (3 → 2).
#
#  • THE EXPONENT AT INFINITY IS `degree`. For a polynomial problem `degree` is the
#    degree of Q, and Q ~ a z^{degree}; the same identity defines `degree` for a
#    rational problem (numerator degree minus the total pole order), so
#    `asymptotic_directions(prob, θ)` is one formula for both.
#
#  • MATHIEU LIVES HERE, NOT ON A CYLINDER. `sw_curve.jl`'s Nekrasov–Shatashvili
#    problem is V = 2Λ²cos 2x, E = u. Substituting w = e^{2ix} and ψ = w^{−1/2}φ,
#        w²ψ″ + wψ′  ⇒  ħ²φ″ = −(Q(x) + ħ²)/(4w²) · φ ,
#    so with Q(x) = Λ²(w + 1/w) − u,
#        Q̃(w) = −(Λ²w² − u w + Λ²)/(4w³)  −  ħ²/(4w²) .
#    The classical part is `mathieu_problem`: an order-3 pole at w = 0, an order-3
#    pole at w = ∞ (e = −1), and exactly TWO simple turning points
#    w = (u ± √(u²−4Λ⁴))/(2Λ²), which collide precisely at the Seiberg–Witten
#    singularities u = ±2Λ². Two order-3 poles ⇒ one marked point on each of two
#    boundary circles ⇒ the annulus Ã(1,1) ⇒ the Kronecker quiver. Pinned by:
#    `mathieu_problem` turning points vs the closed form, and
#    `triangulation_quiver == [0 2; −2 0]`.
#    The ħ² term is a pure shift u → u + ħ² of the double-pole coefficient. It is
#    NOT carried: this layer supplies classical Stokes geometry and classical
#    periods, where it does not enter. It is what an all-orders quantum period at a
#    pole would have to carry (research-questions §S10).
#
#  • THE RICCATI RECURSION STAYS POLYNOMIAL. `wkb_expansion` and everything built on
#    it (`quantum_period`, `quantization`, `double_well`, `weber`) keep their
#    `::SchrodingerProblem` signatures, so a `RationalProblem` reaching them is a
#    MethodError rather than a silent wrong answer. What a rational problem *does*
#    support is the classical layer: turning points, `period_integral`,
#    `stokes_graph`, `saddles`, `ideal_triangulation`, `triangulation_quiver`.

"""
    RationalProblem{T}

A one-dimensional Schrödinger problem ``ħ²ψ″ = Q(z)ψ`` with a **rational** `Q`,

    Q(z) = N(z) / ∏_j (z − p_j)^{m_j} ,

given by the numerator coefficients `n_coeffs` in **ascending** order and the poles
`poles` with their orders `pole_orders`. Poles are exact data, not roots of a stored
denominator (see the ledger at the top of `src/rational_potentials.jl`).

A finite pole must have order `m_j ≥ 2`, and `degree(prob) ≥ −1` at infinity (the
pole order of `Q dz²` there is `degree + 4`, so infinity must stay irregular). A
**simple** pole is an orbifold point and throws [`InvalidPotential`](@ref).

The surface this problem draws its Stokes graph on is the sphere with one singularity
per pole (including infinity). An **irregular** pole (`m_j ≥ 3`) contributes a
boundary circle carrying `m_j − 2` marked points; a **double** pole contributes a
*puncture* - a single interior vertex of the dual triangulation, carrying no marked
points, into which the Stokes lines spiral. See [`asymptotic_directions`](@ref),
[`ideal_triangulation`](@ref) and [`ring_domain_walls`](@ref).

# Examples
```julia
# Mathieu / pure SU(2) at Λ = 1, u = 3: two order-3 poles ⇒ the annulus
prob = mathieu_problem(1.0, 3.0)
turning_points(prob)                       # two simple turning points
RationalProblem([1.0, 0.0, 1.0], [0.0], [3])   # Q = (1 + z²)/z³

# a double pole ⇒ a puncture: the once-punctured triangle, cluster type D₃
RationalProblem([1.0, -0.3, 0.0, 1.0], [0.0], [2])   # Q = (1 − 0.3z + z³)/z²
```
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

# -- punctures: the residue charge and its wall --------------------------------------

"""
    RingDomainWall{F}

The wall carried by a puncture (a double pole of `Q`). `pole` is its index into
[`poles`](@ref), `central_charge` is `Z_p = ∮√Q dz = 2πi√c` for the local
`Q ~ c/(z − p)²`, and `theta` is the critical phase `arg Z_p (mod π)`.

At `theta` the trajectories around the puncture close up into a **ring domain** - a
one-parameter family of closed trajectories - instead of spiralling in, so the Stokes
graph is not generic there and [`ideal_triangulation`](@ref) has nothing to build. Off
the wall they spiral, and the puncture's valence in the dual triangulation is constant
within each chamber. This is a saddle connection's counterpart for a puncture, which
is why it is a separate type: it belongs to no turning-point pair.
"""
struct RingDomainWall{F}
    pole::Int
    central_charge::Complex{F}
    theta::F
end

pole_index(w::RingDomainWall) = w.pole
central_charge(w::RingDomainWall) = w.central_charge
theta(w::RingDomainWall) = w.theta
mass(w::RingDomainWall) = abs(w.central_charge)

"""
    ring_domain_walls(prob) -> Vector{RingDomainWall}

One [`RingDomainWall`](@ref) per puncture of `prob` (empty when there is none): the
residue central charge `Z_p = 2πi√c` of the loop around the double pole, and the phase
`arg Z_p (mod π)` at which its trajectories close into a ring domain.

`stokes_graph` refuses to trace near these phases - the spiral winds
`log(R/r)·|tan(θ − arg√c)|/2π` times before reaching the pole circle, which diverges at
the wall - so this is the list of phases a chamber scan must step over, exactly as
[`saddle_candidates`](@ref) is for turning-point pairs.
"""
function ring_domain_walls(prob::AbstractSchrodingerProblem)
    F = real(typeof(float(one(eltype(poles(prob))) + 0.0)))
    out = RingDomainWall{F}[]
    for j in puncture_indices(prob)
        Z = 2 * F(π) * im * sqrt(Complex{F}(_pole_coefficient(prob, j)))
        push!(out, RingDomainWall{F}(j, Z, mod(angle(Z), F(π))))
    end
    out
end

ring_domain_walls(::SchrodingerProblem{T}) where {T} =
    RingDomainWall{real(typeof(float(one(T))))}[]

# -- named problems ----------------------------------------------------------------

"""
    mathieu_problem(Λ, u) -> RationalProblem

The Nekrasov–Shatashvili problem of pure SU(2) - `sw_curve.jl`'s Mathieu equation
`V = 2Λ²cos 2x`, `E = u` - written on the `w = e^{2ix}` plane, where it is rational:

    Q̃(w) = −(Λ²w² − u·w + Λ²) / (4w³) .

Two order-3 poles (at `w = 0` and at `w = ∞`, since `Q̃ ~ −Λ²/(4w)`), so the surface
is the annulus with one marked point on each boundary circle: cluster type `Ã(1,1)`,
the Kronecker quiver. The two simple turning points are
`w = (u ± √(u²−4Λ⁴))/(2Λ²)`; they collide exactly at the Seiberg–Witten singularities
`u = ±2Λ²`.

The exact substitution carries a further `−ħ²/(4w²)`, i.e. a shift `u → u + ħ²` of the
double-pole coefficient (derived in the ledger at the top of
`src/rational_potentials.jl`). It is absent here because this layer is classical: it
would matter only to an all-orders quantum period at a pole.
"""
function mathieu_problem(Λ, u)
    T = promote_type(typeof(float(Λ)), typeof(float(u)))
    Λ² = T(Λ)^2
    RationalProblem(T[-Λ² / 4, T(u) / 4, -Λ² / 4], T[zero(T)], [3])
end
