# Saddles (BPS states) of a Schrödinger problem.
#
# A *saddle* is a Stokes line connecting two turning points `z_i`, `z_j` — a finite
# trajectory. It appears only at the critical phase `θ_c = arg Z_ij (mod π)`, where
#
#     Z_ij = 2 ∫_{z_i}^{z_j} √Q dz
#
# is the central charge. Under the Iwaki–Nakanishi dictionary a finite Stokes line is a
# BPS state and `|Z_ij|` is its mass — the datum the M4 cluster bridge turns into a
# spectrum. `saddle_candidates` finds the `Z_ij` of every turning-point pair from the
# periods alone; `saddles` keeps only those whose finite edge is confirmed by actually
# tracing the graph at `θ_c`.

# -- type --------------------------------------------------------------------------

"""
    Saddle{F}

A candidate BPS saddle between the simple turning points indexed `pair = (i, j)`
(`i < j`, into `simple_turning_points(prob)`). `central_charge` is `Z_ij = 2∫√Q`, and
`theta` is the critical phase `arg Z_ij (mod π)` at which the finite Stokes line
appears. Its mass is `|Z_ij|` (see [`mass`](@ref)).
"""
struct Saddle{F}
    pair::Tuple{Int,Int}
    central_charge::Complex{F}
    theta::F
end

pair(s::Saddle) = s.pair
central_charge(s::Saddle) = s.central_charge
theta(s::Saddle) = s.theta

"""
    mass(s::Saddle) -> Real

The BPS mass `|Z_ij|` of the saddle.
"""
mass(s::Saddle) = abs(s.central_charge)

# -- central charge of a turning-point pair ----------------------------------------

# Straight (or once-detoured) path from z_i to z_j for the open-path period. Endpoints
# are nudged inward off the turning points (period_integral refuses a vertex on a TP);
# the dropped end-caps are O(ε^{3/2}) and shrink out of the 1e-6 mass target. If a third
# turning point grazes the straight segment, route via a perpendicular midpoint.
function _saddle_path(prob, zi::Complex{F}, zj::Complex{F}, others;
                      nudge = 1e-6, graze = nothing) where {F}
    L = abs(zj - zi)
    gr = graze === nothing ? F(0.05) * L : F(graze)
    dir = (zj - zi) / L
    ε = F(nudge) * L
    a = zi + ε * dir
    b = zj - ε * dir

    # grazing test: distance of each other turning point to the segment [a, b]
    grazed = false
    for z in others
        t = real((z - a) * conj(b - a)) / abs2(b - a)
        if zero(F) ≤ t ≤ one(F)
            foot = a + t * (b - a)
            abs(z - foot) < gr && (grazed = true; break)
        end
    end
    if !grazed
        return Complex{F}[a, b]
    end
    # perpendicular detour: bow the midpoint out by ~half the segment length
    mid = (a + b) / 2 + im * dir * (F(0.5) * L)
    Complex{F}[a, mid, b]
end

"""
    saddle_candidates(prob::SchrodingerProblem) -> Vector{Saddle}

For every pair of simple turning points, the central charge `Z_ij = 2∫√Q` (computed
from the open-path period, with a perpendicular detour when a third turning point
grazes the straight path) and the critical phase `θ_c = arg Z_ij (mod π)`. These are
*candidates* — [`saddles`](@ref) confirms which actually form a finite Stokes line.
"""
function saddle_candidates(prob::SchrodingerProblem)
    tps = simple_turning_points(prob)
    F = isempty(tps) ? Float64 : real(typeof(float(location(first(tps)))))
    zs = [location(t) for t in tps]
    out = Saddle{F}[]
    for i in eachindex(zs), j in (i + 1):length(zs)
        others = [zs[k] for k in eachindex(zs) if k != i && k != j]
        path = _saddle_path(prob, complex(zs[i]), complex(zs[j]), others)
        half = period_integral(prob, path; closed = false)
        Z = 2 * half
        θc = mod(angle(Z), F(π))
        push!(out, Saddle{F}((i, j), Z, θc))
    end
    out
end

# -- graph queries -----------------------------------------------------------------

"""
    is_saddle(g::StokesGraph, pair::Tuple{Int,Int}) -> Bool

Whether the traced graph `g` contains a finite Stokes line between the turning points
`pair = (i, j)` (order-insensitive).
"""
function is_saddle(g::StokesGraph, pr::Tuple{Int,Int})
    (min(pr...), max(pr...)) in edges(g)
end

"""
    saddles(prob::SchrodingerProblem; verify = true, tol = 1e-8, kwargs...) -> Vector{Saddle}

The confirmed saddles of `prob`. Each [`saddle_candidates`](@ref) entry is re-checked by
tracing the Stokes graph at its critical phase `θ_c` (nudged by `tol` to land on the
knife-edge cleanly) and keeping it only if the finite edge actually appears. With
`verify = false` the raw candidates are returned. Extra `kwargs` pass through to
[`stokes_graph`](@ref).
"""
function saddles(prob::SchrodingerProblem; verify::Bool = true, tol = 1e-8, kwargs...)
    cands = saddle_candidates(prob)
    verify || return cands
    kept = eltype(cands)[]
    for s in cands
        g = stokes_graph(prob; theta = s.theta, allow_incomplete = true, kwargs...)
        is_saddle(g, s.pair) && push!(kept, s)
    end
    kept
end

# -- θ-families --------------------------------------------------------------------

"""
    stokes_graph_family(prob::SchrodingerProblem, thetas; kwargs...) -> Vector{StokesGraph}

The Stokes graph traced at each phase in `thetas` — a slice through the wall-crossing
as `θ` sweeps. Topology jumps (an edge appearing or vanishing) mark critical phases;
`allow_incomplete = true` is set so a marginal ray never aborts the sweep.
"""
function stokes_graph_family(prob::SchrodingerProblem, thetas; kwargs...)
    [stokes_graph(prob; theta = θ, allow_incomplete = true, kwargs...) for θ in thetas]
end
