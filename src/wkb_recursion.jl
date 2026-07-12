# The all-orders WKB (Riccati) recursion. With ψ = exp(∫S dz) the Schrödinger equation
# ħ²ψ″ = Qψ becomes the Riccati equation  S² + S′ = ħ⁻² Q  for S = Σ_{m≥-1} ħ^m S_m:
#
#   S₋₁ = √Q,   S₀ = −Q′/(4Q),   S_{m+1} = −(Σ_{i+j=m, i,j≥0} S_i S_j + S_m′) / (2√Q).
#
# Every S_m has the exact shape  S_m = p_m(z) · u^{ε} / Q^{k},  u = √Q,  ε = m mod 2.
# Storing (p_m, k, ε) as a RiccatiTerm keeps all algebra inside a polynomial ring -
# no rational-function field, hence no gcd, so the BigFloat ring stays numerically
# safe. u² = Q is reduced on every multiply; ε only ever moves between 0 and 1.

import AbstractAlgebra as AA

# exactness trait (local; mirrors Resurgence's _is_exact_type without depending on it)
_exact_coeff(::Type{<:Union{Integer,Rational}}) = true
_exact_coeff(::Type{Complex{S}}) where {S} = _exact_coeff(S)
_exact_coeff(::Type{<:AbstractFloat}) = false
_exact_coeff(::Type) = true

"""
    RiccatiTerm{P}

Internal representation of a WKB coefficient `S_m = num(z)·u^{sqrt_pow} / Q^{q_pow}`
with `u = √Q` and `sqrt_pow ∈ {0,1}`; `num` is a polynomial in an AbstractAlgebra
ring. Never exposed in the public API.
"""
struct RiccatiTerm{P}
    num::P
    q_pow::Int
    sqrt_pow::Int
end

# S_i S_j : parities add; u² = Q lowers the denominator power by ⌊(ε_i+ε_j)/2⌋.
function _mul_term(a::RiccatiTerm, b::RiccatiTerm)
    e = a.sqrt_pow + b.sqrt_pow
    RiccatiTerm(a.num * b.num, a.q_pow + b.q_pow - (e ÷ 2), e % 2)
end

# scalar (rational) multiple, lifted into the base ring
_scale_term(a::RiccatiTerm, r, br) = RiccatiTerm(a.num * br(r), a.q_pow, a.sqrt_pow)

# d/dz [num·u^ε/Q^k] = [num′·Q + (ε/2 − k)·num·Q′] · u^ε / Q^{k+1}
function _deriv_term(a::RiccatiTerm, Q, Qp, br)
    ε, k = a.sqrt_pow, a.q_pow
    new_num = AA.derivative(a.num) * Q + a.num * Qp * br(Rational(ε, 2) - k)
    RiccatiTerm(new_num, k + 1, ε)
end

# add two terms of equal parity over the common denominator Q^{max q_pow}
function _add_term(a::RiccatiTerm, b::RiccatiTerm, Q)
    K = max(a.q_pow, b.q_pow)
    num = a.num * Q^(K - a.q_pow) + b.num * Q^(K - b.q_pow)
    RiccatiTerm(num, K, a.sqrt_pow)
end

# divide by 2u:  1/u = u/Q, so ε=0 → ε=1 with q_pow+1; ε=1 → ε=0.
function _div_2u_term(a::RiccatiTerm, br)
    if a.sqrt_pow == 1
        RiccatiTerm(a.num * br(1 // 2), a.q_pow, 0)
    else
        RiccatiTerm(a.num * br(1 // 2), a.q_pow + 1, 1)
    end
end

"""
    WKBExpansion{P}

The WKB/Riccati expansion of a [`SchrodingerProblem`](@ref) to a given `order`:
holds `S₋₁ … S_order` as [`RiccatiTerm`](@ref)s over an AbstractAlgebra ring (exact
`QQ` or `BigFloat` `RealField`, per `arithmetic`). Build with [`wkb_expansion`](@ref);
extract with [`evaluate_s_odd`](@ref) / [`s_odd_terms`](@ref).
"""
struct WKBExpansion{P}
    prob::SchrodingerProblem
    order::Int
    arithmetic::Symbol
    Q::P
    Qp::P
    br::Any
    s::Vector{RiccatiTerm{P}}   # s[m + 2] == S_m, for m in -1:order
end

function _wkb_ring(prob::SchrodingerProblem, arithmetic::Symbol)
    exact = _exact_coeff(eltype(q_coefficients(prob)))
    mode = arithmetic === :auto ? (exact ? :exact : :bigfloat) : arithmetic
    if mode === :exact
        exact || throw(Resurgence.InvalidArgument(
            "arithmetic=:exact needs exact (rational/integer) potential coefficients"))
        R, _ = AA.polynomial_ring(AA.QQ, string(variable(prob)))
        return R, AA.QQ, mode
    elseif mode === :bigfloat
        R, _ = AA.polynomial_ring(AA.RealField, string(variable(prob)))
        return R, AA.RealField, mode
    else
        throw(Resurgence.InvalidArgument(
            "arithmetic must be :auto, :exact, or :bigfloat, got :$arithmetic"))
    end
end

"""
    wkb_expansion(prob::SchrodingerProblem; order = 10, arithmetic = :auto)
        -> WKBExpansion

Compute the Riccati coefficients `S₋₁ … S_order`. `arithmetic` selects the numerator
ring: `:exact` (`QQ`, requires rational/integer `Q`), `:bigfloat` (`RealField`), or
`:auto` (exact when the potential coefficients are exact, else BigFloat). Exact-ring
blow-up in the numerator degree is the known cost driver at high order.
"""
function wkb_expansion(prob::SchrodingerProblem; order::Integer = 10,
                       arithmetic::Symbol = :auto)
    order ≥ 0 || throw(Resurgence.InvalidArgument("order must be ≥ 0, got $order"))
    R, br, mode = _wkb_ring(prob, arithmetic)
    q = q_coefficients(prob)
    Q = R([br(c) for c in q])
    Qp = AA.derivative(Q)
    P = typeof(Q)

    s = Vector{RiccatiTerm{P}}(undef, order + 2)
    s[1] = RiccatiTerm(one(R), 0, 1)                       # S₋₁ = √Q
    s[2] = RiccatiTerm(Qp * br(-1 // 4), 1, 0)             # S₀ = −Q′/(4Q)
    Sm(m) = s[m + 2]
    for m in 0:(order - 1)
        # Σ_{i+j=m, i,j≥0} S_i S_j  (all of parity m mod 2)
        acc = _mul_term(Sm(0), Sm(m))
        for i in 1:m
            acc = _add_term(acc, _mul_term(Sm(i), Sm(m - i)), Q)
        end
        acc = _add_term(acc, _deriv_term(Sm(m), Q, Qp, br), Q)   # + S_m′
        s[m + 3] = _div_2u_term(_scale_term(acc, -1, br), br)    # S_{m+1}
    end
    WKBExpansion{P}(prob, order, mode, Q, Qp, br, s)
end

# -- internal accessors (used by periods.jl / voros.jl) ----------------------------

_s_term(w::WKBExpansion, m::Integer) = w.s[m + 2]

# numerator coefficients of S_m as plain Julia numbers (ascending), AA-free
function _numerator_coeffs(term::RiccatiTerm)
    d = AA.degree(term.num)
    d < 0 ? [] : [AA.coeff(term.num, i) for i in 0:d]
end

# Horner evaluation of an AA polynomial at a Julia complex z (→ Complex{F})
function _poly_eval(num, z::Complex{F}) where {F}
    d = AA.degree(num)
    acc = zero(Complex{F})
    for i in d:-1:0
        acc = acc * z + Complex{F}(AA.coeff(num, i))
    end
    acc
end

# value of a RiccatiTerm at z, given Q(z) and a chosen branch u = √Q(z)
function _eval_term(term::RiccatiTerm, z::Complex{F}, Qz, u) where {F}
    v = _poly_eval(term.num, z)
    term.sqrt_pow == 1 ? v * u / Qz^term.q_pow : v / Qz^term.q_pow
end

# -- public extraction -------------------------------------------------------------

"""
    order(w::WKBExpansion) -> Int

The truncation order (the highest computed `S_m`).
"""
order(w::WKBExpansion) = w.order

"""
    s_odd_terms(w::WKBExpansion) -> Vector{<:NamedTuple}

The odd WKB coefficients `S_{-1}, S_1, S_3, …` as AbstractAlgebra-free descriptions:
each a `(m, numerator, q_power)` with `numerator` the ascending coefficient vector of
`p_m(z)` (plain numbers) and `S_m = p_m(z)·√Q / Q^{q_power}` (`sqrt_pow` is always 1
for odd `m`).
"""
function s_odd_terms(w::WKBExpansion)
    [(m = m, numerator = _numerator_coeffs(_s_term(w, m)), q_power = _s_term(w, m).q_pow)
     for m in -1:2:w.order]
end

"""
    evaluate_s_odd(w::WKBExpansion, z; branch = 1, sqrt_q = nothing) -> FormalSeries

The odd part `S_odd(z, ħ) = Σ_{m odd} S_m(z) ħ^m` evaluated at `z`, returned as a
`Resurgence.FormalSeries` in `:ħ` with `power_offset = -1` and zero even coefficients -
exactly the integrand-at-a-point of the Voros/quantum-period machinery. `branch =
±1` picks the sign of the principal `√Q(z)`; pass `sqrt_q` to supply a branch value
tracked from elsewhere (used by [`period_integral`](@ref)).
"""
function evaluate_s_odd(w::WKBExpansion, z::Number; branch::Integer = 1,
                        sqrt_q = nothing)
    F = typeof(float(real(z)))
    zc = Complex{F}(z)
    Qz = w.prob(zc)
    u = sqrt_q === nothing ? branch * sqrt(Qz) : sqrt_q
    coeffs = zeros(Complex{F}, w.order + 2)                # m = -1 … order
    for m in -1:2:w.order
        coeffs[m + 2] = _eval_term(_s_term(w, m), zc, Qz, u)
    end
    Resurgence.FormalSeries(coeffs, :ħ; power_offset = -1 // 1)
end

"""
    even_odd_residual(w::WKBExpansion; z0 = nothing) -> Real

Numerical check of the exact reduction `S_even = −½ ∂_z log S_odd`, in the
division-free form `S_even·S_odd + ½ S_odd′ = 0`. Returns the largest absolute
coefficient of the left-hand side as an ħ-series evaluated at a generic point `z0`
(default a fixed off-axis complex point); ≈ 0 to the computed order.
"""
function even_odd_residual(w::WKBExpansion; z0 = nothing)
    z = z0 === nothing ? Complex{BigFloat}(3 // 10, 7 // 10) : Complex{BigFloat}(z0)
    Qz = w.prob(z)
    u = sqrt(Qz)
    N = w.order + 2
    even = zeros(Complex{BigFloat}, N)      # offset 0: S_0, S_1(skip)…; store by m
    odd = zeros(Complex{BigFloat}, N)       # offset -1
    oddp = zeros(Complex{BigFloat}, N)      # offset -1, derivative
    for m in -1:w.order
        val = _eval_term(_s_term(w, m), z, Qz, u)
        if iseven(m)
            m ≥ 0 && (even[m + 1] = val)    # S_even series, offset 0
        else
            odd[m + 2] = val                # S_odd series, offset -1
            dterm = _deriv_term(_s_term(w, m), w.Q, w.Qp, w.br)
            oddp[m + 2] = _eval_term(dterm, z, Qz, u)
        end
    end
    Se = Resurgence.FormalSeries(even, :ħ; power_offset = 0 // 1)
    So = Resurgence.FormalSeries(odd, :ħ; power_offset = -1 // 1)
    Sop = Resurgence.FormalSeries(oddp, :ħ; power_offset = -1 // 1)
    resid = Se * So + (1 // 2) * Sop
    maximum(abs, Resurgence.coefficients(resid))
end
