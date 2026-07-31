# Stokes graphs of a Schrödinger problem.
#
# A Stokes line emanating from a turning point z₀ is the locus
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
# Degenerate turning points. At a zero of order `m`, `Q ≈ c_m (z−z₀)^m` with
# `c_m = Q^{(m)}(z₀)/m!`, so `∫√Q ≈ √c_m (z−z₀)^{m/2+1}/(m/2+1)` and the Stokes
# condition puts `m+2` rays at `arg(z−z₀) = (2θ − arg c_m + 2πk)/(m+2)`, `k = 0…m+1`
# (`m = 1` is the familiar trivalent case). For **even** `m` the point is not a branch
# point of `w² = Q` - `√Q` is single-valued there and the spectral curve has a node -
# and opposite rays sit exactly `π` apart, so a double turning point is a *crossing of
# two smooth Stokes curves* rather than a vertex. The augmented-state tracer needs no
# other change: it never makes a discrete branch choice, which is exactly why it
# generalizes. See [BNR82] for the new Stokes curves such points emit.
#
# The tracer runs in `Float64` by default: a Stokes graph's content is its *topology*,
# a discrete datum, so double precision is plenty; `bigfloat = true` opts in to the
# slow high-precision trace when a marginal saddle needs it.

import OrdinaryDiffEqTsit5 as ODE

# -- types -------------------------------------------------------------------------

"""
    StokesLine{F}

One traced Stokes line. `source` is the index (into the graph's `turning_points`) of
the turning point it emanates from, `direction ∈ 0:(m+1)` selects which of the `m+2`
rays of an order-`m` point, and `points` is the traced polyline in the `z`-plane
(`m = 1` gives the usual three rays). `endpoint` is one of
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

# The leading Taylor coefficient c_m = Q^{(m)}(z₀)/m! of an order-m turning point, in
# the tracer's float type. `q_taylor_at` is exact when the coefficients are.
function _leading_coeff(prob, z0::Complex{F}, m::Int) where {F}
    Complex{F}(q_taylor_at(prob, z0)[m + 1])
end

# Trace a single ray. Returns (points, endpoint, target_idx, mass). `nrays = m + 2` and
# `c` is the order-m leading coefficient; `m = 1` recovers `nrays = 3`, `c = Q′(z₀)`.
function _trace_ray(prob, tps::Vector{Complex{F}}, src::Int, k::Int, θ::F;
                    nrays::Int, c::Complex{F},
                    seed_radius, hit_radius, escape_radius, max_mass,
                    reltol, abstol) where {F}
    z0 = tps[src]
    dir = (2θ - angle(c) + 2 * F(π) * k) / nrays
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

Trace the Stokes graph of `prob` at phase `theta`. At a turning point `z₀` of order `m`
(so `Q ≈ c_m (z−z₀)^m`), `m + 2` rays are seeded in the directions
`arg(z − z₀) = (2θ − arg c_m + 2πk)/(m+2)`, `k = 0 … m+1`, and integrated until they
escape to infinity, run into another turning point (a saddle), or hit the `max_mass`
safety. A simple turning point gives the familiar three rays with `c_1 = Q′(z₀)`.

**Degenerate turning points** (order ≥ 2) are traced, not refused: see
[`is_degenerate`](@ref). An order-2 point emits four rays and is a *crossing* of two
smooth Stokes curves rather than a vertex, because `√Q` is single-valued there. The
downstream triangulation layer still requires an all-simple graph.

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
    isempty(all_tps) && throw(InvalidPotential("no turning points to trace"))

    tps_F = [TurningPoint{F}(Complex{F}(location(t)), order(t)) for t in all_tps]
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
    for s in eachindex(tps_z)
        ms = order(tps_F[s])
        nrays = ms + 2
        c = _leading_coeff(prob, tps_z[s], ms)
        for k in 0:(nrays - 1)
            pts, ep, tgt, m = _trace_ray(prob, tps_z, s, k, θ; nrays, c,
                                         seed_radius = sr, hit_radius = hr,
                                         escape_radius = er, max_mass = mm,
                                         reltol = rt, abstol = at)
            if ep === :incomplete && !allow_incomplete
                throw(TracingFailed("ray $k from turning point $s reached max_mass = " *
                                    "$mm without escaping or hitting a turning point"))
            end
            push!(ls, StokesLine{F}(s, k, pts, ep, tgt, m))
        end
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
    is_degenerate(g::StokesGraph) -> Bool

`true` when `g` has a turning point of order ≥ 2. Such a graph is a *wall*: the
degenerate point splits into simple ones under any generic perturbation, and the
`(m+2)`-gon it is dual to admits several triangulations. The triangulation and bridge
layers require `!is_degenerate(g)`.
"""
is_degenerate(g::StokesGraph) = any(!is_simple, g.turning_points)

"""
    turning_point_orders(g::StokesGraph) -> Vector{Int}

The multiplicities of `g`'s turning points, sorted ascending. `all(== 1)` for a
generic graph; see [`is_degenerate`](@ref).
"""
turning_point_orders(g::StokesGraph) = sort!([order(t) for t in g.turning_points])

"""
    topology_signature(g::StokesGraph) -> Tuple

A canonical, tolerance-robust invariant of the graph's topology: the sorted
[`turning_point_orders`](@ref), the number of infinite rays, and the sorted saddle
edges (turning points ordered canonically by position). Two graphs with the same
signature have the same Stokes topology - this is the datum the cluster bridge reads.
The orders are carried so a degenerate graph is never silently equated with a generic
one of the same turning-point count.
"""
function topology_signature(g::StokesGraph)
    (turning_point_orders(g), n_infinite_lines(g), edges(g))
end
