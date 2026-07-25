# Stokes graphs of a Schrödinger problem.
#
# A Stokes line emanating from a simple turning point z₀ is the locus
# `Im[e^{−iθ} ∫_{z₀}^z √Q dz] = 0`. We trace it with an *augmented state* `(z, w)`,
# `w = √Q(z)`, evolving under
#
#     dz/dt = e^{iθ} / w ,      dw/dt = Q′(z) e^{iθ} / (2 w²) .
#
# Then `d(∫√Q dz)/dt = e^{iθ}`, so `e^{−iθ}∫√Q dz = t + const` stays on the real axis
# and `t` is *exactly* the mass parameter - no discrete branch choices are ever made,
# which is the whole point (branch cuts are the top correctness risk). A per-step
# re-projection `w → ±√Q(z)` snaps the tiny integration drift back onto the shell.
#
# The tracer runs in `Float64` by default: a Stokes graph's content is its *topology*,
# a discrete datum, so double precision is plenty; `bigfloat = true` opts in to the
# slow high-precision trace when a marginal saddle needs it.

import OrdinaryDiffEqTsit5 as ODE

# -- types -------------------------------------------------------------------------

"""
    StokesLine{F}

One traced Stokes line. `source` is the index (into the graph's `turning_points`) of
the simple turning point it emanates from, `direction ∈ 0:2` selects which of the
three rays, and `points` is the traced polyline in the `z`-plane. `endpoint` is one of
`:turning_point` (a saddle - `target` is the index of the turning point it runs into),
`:infinity` (an infinite ray - `target === nothing`), or `:incomplete` (tracing gave
up at `max_mass`). `mass` is the accumulated mass parameter `t` at the endpoint.
"""
struct StokesLine{F}
    source::Int
    direction::Int
    points::Vector{Complex{F}}
    endpoint::Symbol
    target::Union{Int,Nothing}
    mass::F
end

source(l::StokesLine) = l.source
direction(l::StokesLine) = l.direction
points(l::StokesLine) = l.points
"""
    endpoint(l::StokesLine) -> Symbol

How the traced line ended: `:turning_point` (a saddle connection - it ran into another
turning point), `:infinity` (it escaped) or `:incomplete` (tracing stopped early).
"""
endpoint(l::StokesLine) = l.endpoint
target(l::StokesLine) = l.target
mass(l::StokesLine) = l.mass

"""
    is_finite_line(l::StokesLine) -> Bool

`true` for a saddle trajectory (a Stokes line running between two turning points).
"""
is_finite_line(l::StokesLine) = l.endpoint === :turning_point

"""
    StokesGraph{F}

The Stokes graph of a Schrödinger problem at phase `theta`: its simple
`turning_points` and the traced `lines`. Query its topology with [`finite_lines`](@ref),
[`edges`](@ref), [`n_infinite_lines`](@ref) and [`topology_signature`](@ref).
"""
struct StokesGraph{F}
    problem::SchrodingerProblem
    theta::F
    turning_points::Vector{TurningPoint{F}}
    lines::Vector{StokesLine{F}}
end

theta(g::StokesGraph) = g.theta
turning_points(g::StokesGraph) = g.turning_points
lines(g::StokesGraph) = g.lines

# -- the tracer --------------------------------------------------------------------

# Evaluate Q and Q′ at a complex point from the (real-coefficient) problem.
@inline _Q(prob, z) = prob(z)
@inline _dQ(prob, z) = q_derivative_at(prob, z)

# Trace a single ray. Returns (points, endpoint, target_idx, mass).
function _trace_ray(prob, tps::Vector{Complex{F}}, src::Int, k::Int, θ::F;
                    seed_radius, hit_radius, escape_radius, max_mass,
                    reltol, abstol) where {F}
    z0 = tps[src]
    c = _dQ(prob, z0)                      # Q′(z₀); c ≠ 0 at a simple turning point
    dir = (2θ - angle(c) + 2 * F(π) * k) / 3
    zseed = z0 + F(seed_radius) * cis(dir)
    w = sqrt(_Q(prob, zseed))
    # Orient the seed so the ray leaves z₀ outward: dz/dt = e^{iθ}/w must point along dir.
    if real((cis(θ) / w) * cis(-dir)) < 0
        w = -w
    end

    u0 = F[real(zseed), imag(zseed), real(w), imag(w)]
    eiθ = cis(θ)

    function rhs!(du, u, p, t)
        z = complex(u[1], u[2]); ww = complex(u[3], u[4])
        dz = eiθ / ww
        dw = _dQ(prob, z) * eiθ / (2 * ww^2)
        du[1] = real(dz); du[2] = imag(dz); du[3] = real(dw); du[4] = imag(dw)
        return nothing
    end

    # Re-projection: snap w back onto √Q(z), keeping the current sheet.
    function reproject!(integ)
        u = integ.u
        z = complex(u[1], u[2]); ww = complex(u[3], u[4])
        wnew = sqrt(_Q(prob, z))
        if real(wnew * conj(ww)) < 0
            wnew = -wnew
        end
        u[3] = real(wnew); u[4] = imag(wnew)
        return nothing
    end
    reproj = ODE.DiscreteCallback((u, t, integ) -> true, reproject!;
                                  save_positions = (false, false))

    esc = ODE.ContinuousCallback(
        (u, t, integ) -> F(escape_radius)^2 - (u[1]^2 + u[2]^2), ODE.terminate!)

    others = [j for j in eachindex(tps) if j != src]
    hit = ODE.ContinuousCallback(
        function (u, t, integ)
            z = complex(u[1], u[2])
            isempty(others) && return one(F)
            minimum(abs2(z - tps[j]) for j in others) - F(hit_radius)^2
        end, ODE.terminate!)

    cbs = ODE.CallbackSet(reproj, esc, hit)
    prob_ode = ODE.ODEProblem(rhs!, u0, (zero(F), F(max_mass)))
    sol = ODE.solve(prob_ode, ODE.Tsit5(); callback = cbs,
                    reltol = F(reltol), abstol = F(abstol),
                    save_everystep = true, dense = false)

    pts = [complex(u[1], u[2]) for u in sol.u]
    zend = pts[end]
    tend = F(sol.t[end])

    if abs(zend) >= F(escape_radius) * (1 - 1e-6)
        return pts, :infinity, nothing, tend
    end
    if !isempty(others)
        d, j = findmin([abs(zend - tps[o]) for o in others])
        if d < F(hit_radius) * (1 + 1e-6)
            return pts, :turning_point, others[j], tend
        end
    end
    return pts, :incomplete, nothing, tend
end

# -- public entry ------------------------------------------------------------------

"""
    stokes_graph(prob::SchrodingerProblem; theta = 0.0, bigfloat = false, kwargs...)

Trace the Stokes graph of `prob` at phase `theta`. Three rays are seeded at each simple
turning point `z₀` in the directions `arg(z − z₀) = (2θ − arg Q′(z₀) + 2πk)/3`,
`k = 0,1,2`, and integrated until they escape to infinity, run into another turning
point (a saddle), or hit the `max_mass` safety.

Turning points of order ≥ 2 are classified fine but cannot be traced through - they
throw [`UnsupportedTurningPoint`](@ref) (Weber local models are deferred).

Keyword arguments (all with geometry-scaled defaults): `seed_radius`, `hit_radius`,
`escape_radius`, `max_mass`, `reltol`, `abstol`, and `allow_incomplete` (default
`false` - an `:incomplete` ray otherwise throws [`TracingFailed`](@ref)).
"""
function stokes_graph(prob::SchrodingerProblem; theta = 0.0, bigfloat::Bool = false,
                      seed_radius = nothing, hit_radius = nothing,
                      escape_radius = nothing, max_mass = nothing,
                      reltol = nothing, abstol = nothing,
                      allow_incomplete::Bool = false)
    F = bigfloat ? BigFloat : Float64
    all_tps = turning_points(prob)
    simple = filter(is_simple, all_tps)
    bad = filter(t -> !is_simple(t), all_tps)
    if !isempty(bad)
        b = first(bad)
        throw(UnsupportedTurningPoint(location(b), order(b)))
    end
    isempty(simple) && throw(InvalidPotential("no simple turning points to trace"))

    tps_F = [TurningPoint{F}(Complex{F}(location(t)), order(t)) for t in simple]
    tps_z = [t.z for t in tps_F]
    θ = F(theta)

    # Geometry scale: pairwise turning-point spread, floored at 1.
    scale = one(F)
    for i in eachindex(tps_z), j in (i + 1):length(tps_z)
        scale = max(scale, abs(tps_z[i] - tps_z[j]))
    end
    scale = max(scale, maximum(abs, tps_z; init = one(F)))

    sr = F(something(seed_radius, 1e-2 * scale))
    hr = F(something(hit_radius, sr / 20))
    er = F(something(escape_radius, 5 * scale))
    # Mass to reach `er` along an asymptotic ray grows like ``√|a_d| R^{d/2+1}/(d/2+1)``
    # (``Q ~ a_d z^d`` at infinity, mass = ∫|√Q|d|z|). `max_mass` must comfortably
    # exceed it or an escaping ray would be cut off mid-flight and misread as
    # `:incomplete`; the curved path buys the ×4 headroom.
    qc = q_coefficients(prob)
    d = length(qc) - 1
    p = F(d) / 2 + 1
    escape_mass = sqrt(abs(F(last(qc)))) * er^p / p
    mm = F(something(max_mass, 4 * escape_mass))
    rt = F(something(reltol, bigfloat ? 1e-20 : 1e-9))
    at = F(something(abstol, bigfloat ? 1e-24 : 1e-12))

    ls = StokesLine{F}[]
    for s in eachindex(tps_z), k in 0:2
        pts, ep, tgt, m = _trace_ray(prob, tps_z, s, k, θ;
                                     seed_radius = sr, hit_radius = hr,
                                     escape_radius = er, max_mass = mm,
                                     reltol = rt, abstol = at)
        if ep === :incomplete && !allow_incomplete
            throw(TracingFailed("ray $k from turning point $s reached max_mass = " *
                                "$mm without escaping or hitting a turning point"))
        end
        push!(ls, StokesLine{F}(s, k, pts, ep, tgt, m))
    end
    StokesGraph{F}(prob, θ, tps_F, ls)
end

# -- topology API ------------------------------------------------------------------

"""
    finite_lines(g::StokesGraph) -> Vector{StokesLine}

The saddle trajectories (Stokes lines connecting two turning points). Each saddle is
traced twice - once from each endpoint - so this list double-counts; use [`edges`](@ref)
for the deduplicated turning-point pairs.
"""
finite_lines(g::StokesGraph) = filter(is_finite_line, g.lines)

"""
    edges(g::StokesGraph) -> Vector{Tuple{Int,Int}}

The deduplicated saddle edges as sorted turning-point index pairs `(i, j)`, `i < j`.
"""
function edges(g::StokesGraph)
    es = Set{Tuple{Int,Int}}()
    for l in finite_lines(g)
        i, j = l.source, l.target::Int
        push!(es, (min(i, j), max(i, j)))
    end
    sort!(collect(es))
end

"""
    n_infinite_lines(g::StokesGraph) -> Int

Number of Stokes lines escaping to infinity.
"""
n_infinite_lines(g::StokesGraph) = count(l -> l.endpoint === :infinity, g.lines)

"""
    topology_signature(g::StokesGraph) -> Tuple

A canonical, tolerance-robust invariant of the graph's topology: the number of simple
turning points, the number of infinite rays, and the sorted saddle edges (turning
points ordered canonically by position). Two graphs with the same signature have the
same Stokes topology - this is the datum the cluster bridge will read.
"""
function topology_signature(g::StokesGraph)
    (length(g.turning_points), n_infinite_lines(g), edges(g))
end
