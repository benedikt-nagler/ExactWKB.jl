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
_wkb_float(prob::SchrodingerProblem) = _wkb_float(eltype(q_coefficients(prob)))

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

# Newton polish of a root of the exact Q, starting from `z0`; bails out at a critical
# point (multiple root, `Q′ ≈ 0`) where Newton no longer converges quadratically.
function _polish_root(prob::SchrodingerProblem, z0::Complex{F}; iters = 50) where {F}
    z = z0
    for _ in 1:iters
        qz = prob(z)
        dq = q_derivative_at(prob, z)
        (iszero(dq) || abs(dq) ≤ eps(F) * (1 + abs(qz))) && break
        Δ = qz / dq
        z -= Δ
        abs(Δ) ≤ eps(F) * (1 + abs(z)) && break
    end
    z
end

# Greedy single-linkage clustering: roots within `tol` of a cluster seed are merged;
# the representative is the cluster mean and the order is the member count.
function _cluster_roots(rts::Vector{Complex{F}}, tol) where {F}
    used = falses(length(rts))
    groups = Tuple{Complex{F},Int}[]
    for i in eachindex(rts)
        used[i] && continue
        members = Complex{F}[rts[i]]
        used[i] = true
        for j in (i + 1):length(rts)
            used[j] && continue
            if abs(rts[j] - rts[i]) ≤ tol
                push!(members, rts[j])
                used[j] = true
            end
        end
        push!(groups, (sum(members) / length(members), length(members)))
    end
    groups
end

"""
    turning_points(prob::SchrodingerProblem; cluster_tol = nothing) -> Vector{TurningPoint}

The turning points of `prob`: the zeros of ``Q = V − E``, with multiplicities. Roots
are located with `PolynomialRoots.roots` in `Complex{F}` (`F` = the working float type
of the coefficients - `BigFloat` for exact potentials, at the caller's
`setprecision`), polished by Newton against the exact `Q`, then clustered to collapse
multiple roots. Sorted deterministically by `(real, imag)`.

`cluster_tol` overrides the merge tolerance (default `8·eps(F)^{1/3}`, which separates
`O(1)`-spaced roots while merging a multiple root's `eps^{1/2}`-split copies).
"""
function turning_points(prob::SchrodingerProblem; cluster_tol = nothing)
    F = _wkb_float(prob)
    q = q_coefficients(prob)
    c = Complex{F}[Complex{F}(x) for x in q]
    rts = length(c) == 2 ? Complex{F}[-c[1] / c[2]] : PolynomialRoots.roots(c)
    polished = Complex{F}[_polish_root(prob, z) for z in rts]
    tol = cluster_tol === nothing ? 8 * eps(F)^(1 // 3) : F(cluster_tol)
    groups = _cluster_roots(polished, tol)
    tps = [TurningPoint{F}(z, m) for (z, m) in groups]
    sort!(tps; by = t -> (real(t.z), imag(t.z)))
    tps
end

"""
    turning_points(prob; ...) restricted with `simple = true`

The simple turning points only (order 1).
"""
simple_turning_points(prob::SchrodingerProblem; kwargs...) =
    filter(is_simple, turning_points(prob; kwargs...))
