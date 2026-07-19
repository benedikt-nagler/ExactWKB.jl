# The signed frame — the Iwaki–Nakanishi orientation layer on top of the DDP decay rule.
#
# `_oriented_contour` (src/charge_lattice.jl) fixes a *representative* of ±γ_e by the
# uniform decay rule, not the physical charge. The physical cycle is
#
#     γ_e^phys = ε_e · γ_e^dec,      ε_e ∈ {±1},
#
# and everything downstream follows with no ledger change: `Z_j^phys = ε_j Z_j^dec`,
# `⟨γ_i^phys, γ_j^phys⟩ = ε_i ε_j P[i,j]`, and ledger item 5 (`B = −P`, src/ddp.jl)
# applies verbatim to the *signed* pairing. Contours are never reoriented — ledger
# item 4 and the M5 orientation ledger (src/quantum_periods.jl) stay literal; the sign
# is bookkeeping layered on an unchanged numerical frame.
#
# ── how ε is determined ────────────────────────────────────────────────────────────
# By inverting the keystone identity, which the M4 layer had weakened to `|P| = |B|`:
#
#     P[i,j] == −ε_i ε_j · B[i,j],      B = triangulation_quiver(t).B
#
# Each quiver edge forces the product ε_i ε_j, so ε is the solution of a linear system
# over 𝔽₂ on the quiver graph — obtained by spanning-tree propagation, unique up to one
# global sign per connected component. Solvability is a genuine theorem, not an
# accounting identity: it requires the product of edge signs around every cycle of the
# quiver graph to be +1. Tree chambers have no cycles, so it is automatic there; a
# *cyclic* chamber (a triangle of the triangulation with all three edges internal,
# generic for degree ≥ 5) imposes one real cocycle condition, and it is exactly this
# condition that the M4 layer could not meet with the uniform decay frame. It holds in
# every chamber of the cubic/quartic/quintic fixtures — see test_signed_frame.jl.
#
# ── the gauge, and why the walk lives in the tests ─────────────────────────────────
# The residual global sign per component is *cosmetic*: it leaves the signed pairing,
# the seed, every mass `|Z|` and every wall phase `arg Z mod π` invariant, so the BPS
# spectrum does not depend on it — only the displayed sign of a central charge does.
# `signed_frame` therefore fixes it by convention (`ε = +1` on the first diagonal of
# each component) and does no extra numerics.
#
# ε is *chamber-local* data. There is no transport rule that carries it across a wall:
# a flip changes the charge lattice basis itself (γ'(2,4) = γ(2,3) − γ(3,4) at the
# quintic (2,3) wall), so the ε-ratio across a flip is not a function of `b_kj`; and
# even a diagonal whose turning-point pair *persists* across the wall can change its ε
# (cubic θ ≈ 0.785), because the seed is a different object on the two sides even
# though the cycle is literally the same. What connects the chambers is therefore the
# combinatorial walk itself, `verify_signed_frame`: from the reference chamber it
# flips its way to the target and lands on exactly the triangulation the tracer
# produces there. That is a test-suite oracle, not a hot path.
#
# ── the reference chamber (a pin, not a derivation) ────────────────────────────────
# The reference is the *top chamber* (θ_max, π), θ_max the largest confirmed wall
# phase below π, where **ε ≡ +1**. This is pinned, not derived: it is the gauge in
# which DDP ledger items 4/5/6 were originally fixed by the cubic oracle, and as
# θ₀ → π⁻ the chamber-relative MGS order degenerates to the absolute θ-decreasing
# order of ledger item 6. It is self-checking — were it wrong, the transported frame
# would fail the tightened keystone at the very first chamber.

"""
    chamber_walls(prob::SchrodingerProblem; kwargs...) -> Vector{Tuple{F,Vector{Tuple{Int,Int}}}}

The walls of the θ-circle: one entry `(θ_c, pairs)` per *confirmed* saddle phase, in
increasing order, listing the turning-point pairs whose saddles sit at that phase.
Several pairs at one phase is a *degenerate* wall — crossing it flips several
diagonals at once (they must commute in the quiver).

Phases are reduced **mod π** (the wall circle has circumference π, as in
[`phase`](@ref)): a saddle at `θ_c = π` is the same wall as one at `θ_c = 0`, and the
two are merged into one degenerate wall at `0`.

Only confirmed [`saddles`](@ref) are walls; an unconfirmed
[`saddle_candidates`](@ref) phase leaves the triangulation unchanged.
"""
function chamber_walls(prob::SchrodingerProblem; kwargs...)
    F = _wkb_float(prob)
    sads = saddles(prob; kwargs...)
    isempty(sads) && return Tuple{F,Vector{Tuple{Int,Int}}}[]
    tol = sqrt(eps(F))
    # reduce onto the wall circle [0, π), snapping a phase within tol of π back to 0
    reduce_phase(θ) = (r = mod(F(θ), F(π)); abs(r - F(π)) ≤ tol * F(π) ? zero(F) : r)
    entries = [(reduce_phase(theta(s)), pair(s)) for s in sads]
    sort!(entries; by = first)
    walls = Tuple{F,Vector{Tuple{Int,Int}}}[]
    for (θ, p) in entries
        if !isempty(walls) && abs(walls[end][1] - θ) ≤ tol * (1 + abs(θ))
            push!(walls[end][2], p)
        else
            push!(walls, (θ, [p]))
        end
    end
    walls
end

"""
    reference_theta(prob::SchrodingerProblem; kwargs...) -> Real

The phase of the *reference chamber* — the midpoint of `(θ_max, π)`, where `θ_max` is
the largest wall phase below `π`. This is the chamber in which the signed frame is
pinned to `ε ≡ +1`. Throws [`ChamberError`](@ref) if the top chamber is too narrow to
trace generically.
"""
function reference_theta(prob::SchrodingerProblem; kwargs...)
    F = _wkb_float(prob)
    walls = chamber_walls(prob; kwargs...)
    θmax = isempty(walls) ? zero(F) : maximum(w[1] for w in walls)
    gap = F(π) - θmax
    gap > sqrt(eps(F)) || throw(ChamberError(
        "the top chamber (θ_max, π) has width $(Float64(gap)) — too narrow to trace " *
        "generically, so the signed frame has no reference chamber to pin ε ≡ +1"))
    θmax + gap / 2
end

# Solve `P[i,j] == −ε_i ε_j B[i,j]` for ε by spanning-tree propagation over the quiver
# graph. Returns `nothing` if an edge magnitude disagrees or a cycle fails to close
# (the cocycle condition); otherwise ε with `+1` on the first vertex of each component.
function _solve_signs(P::AbstractMatrix{<:Integer}, B::AbstractMatrix{<:Integer})
    m = size(B, 1)
    ε = zeros(Int, m)
    for root in 1:m
        ε[root] == 0 || continue
        ε[root] = 1
        stack = [root]
        while !isempty(stack)
            i = pop!(stack)
            for j in 1:m
                B[i, j] == 0 && continue
                abs(P[i, j]) == abs(B[i, j]) || return nothing   # magnitude mismatch
                σ = -P[i, j] ÷ B[i, j]                           # = ε_i ε_j
                abs(σ) == 1 || return nothing
                if ε[j] == 0
                    ε[j] = σ * ε[i]
                    push!(stack, j)
                elseif ε[j] != σ * ε[i]
                    return nothing                               # cocycle failure
                end
            end
        end
    end
    ε
end

"""
    signed_frame(P::AbstractMatrix, t::IdealTriangulation) -> Vector{Int}
    signed_frame(prob::SchrodingerProblem, θ; kwargs...) -> Vector{Int}

The orientation signs `ε` of the charge-lattice basis: `γ_e^phys = ε_e · γ_e^dec`.
Obtained by inverting the keystone identity `P[i,j] == −ε_i ε_j · B[i,j]` against the
combinatorial exchange matrix of the triangulation, with the residual (physically
inert) gauge fixed to `ε = +1` on the first diagonal of each connected component.

Throws [`ContourError`](@ref) when the identity has no solution — either a magnitude
mismatch (the numerical contours do not represent the combinatorial classes) or a
failure of the cocycle condition around a cycle of the quiver graph.
"""
function signed_frame(P::AbstractMatrix{<:Integer}, t::IdealTriangulation)
    B = triangulation_quiver(t).B
    ε = _solve_signs(P, B)
    ε === nothing && throw(ContourError(
        "no orientation ε solves the keystone identity P = −εBε with P = $P and " *
        "B = $B: either |P| ≠ |B| (the numerical contours do not represent the " *
        "combinatorial classes) or the sign cocycle fails to close around a cycle " *
        "of the quiver graph"))
    ε
end

function signed_frame(prob::SchrodingerProblem, θ::Real; margin = nothing,
                      n::Integer = 32, rtol = nothing, kwargs...)
    t = ideal_triangulation(stokes_graph(prob; theta = θ, kwargs...))
    signs(charge_basis(prob, t; margin, n, rtol))
end

"""
    verify_signed_frame(prob::SchrodingerProblem, θ; kwargs...) -> NamedTuple

Prove that the signed frame at `θ` is connected to the reference chamber: walk from
[`reference_theta`](@ref) down to `θ` across the intervening walls, flipping
combinatorially ([`flip`](@ref)), and at each wall transport the frame's gauge by
continuity of the central charge.

Returns `(; n_walls, triangulation_matches, walked)`, where `triangulation_matches`
says the walk lands on exactly the triangulation the tracer produces at `θ` — the
substantive claim, and what connects every chamber to the verified reference one.

**Two things are deliberately not asserted, because neither is a theorem.**

  * *"the transported ε equals the local solve up to component sign"* is vacuous: the
    local solve is only defined up to that sign.
  * *"ε is continuous across a wall on diagonals whose turning-point pair persists"*
    is false, and measurably so (cubic θ ≈ 0.785; quintic below the top chamber).
    The cycle γ_p does persist — turning points do not move with θ — but the *seed*
    is a different object in each chamber, so nothing forces the sign attaching γ_p
    to its vertex to carry over. **ε is chamber-local data, fixed by the keystone.**

So the gauge is fixed by convention, not transported; it is inert anyway (see the
header). The content of ε is verified by the keystone identity in every chamber and,
physically, by chamber-independence of the BPS spectrum in `test_bps.jl`.
"""
function verify_signed_frame(prob::SchrodingerProblem, θ::Real; margin = nothing,
                             n::Integer = 32, rtol = nothing, kwargs...)
    θref = reference_theta(prob; kwargs...)
    walls = chamber_walls(prob; kwargs...)
    crossed = sort!([w for w in walls if θ < w[1] < θref]; by = first, rev = true)

    t = ideal_triangulation(stokes_graph(prob; theta = θref, kwargs...))
    ε = ones(Int, n_diagonals(t))          # the reference-chamber pin

    for (θw, pairs) in crossed
        ks = [findfirst(==(p), t.diagonal_tp_pair) for p in pairs]
        any(isnothing, ks) && throw(ChamberError(
            "the wall at θ = $(Float64(θw)) carries saddle pairs $pairs, but the " *
            "triangulation being walked has diagonals $(t.diagonal_tp_pair)"))
        Bq = triangulation_quiver(t).B
        all(Bq[a, b] == 0 for a in ks, b in ks if a != b) || throw(ChamberError(
            "degenerate wall at θ = $(Float64(θw)): the diagonals $ks flipped " *
            "together do not commute in the quiver, so the crossing is not a " *
            "well-defined simultaneous flip (marginal stability)"))

        tnew = t
        for k in ks
            tnew = flip(tnew, k; direction = -1)
        end
        tnew, _ = canonical_reorder(tnew)

        ε, t = signs(charge_basis(prob, tnew; margin, n, rtol)), tnew
    end

    ttrace = ideal_triangulation(stokes_graph(prob; theta = θ, kwargs...))
    (; n_walls = length(crossed),
       triangulation_matches = t.diagonal_tp_pair == ttrace.diagonal_tp_pair &&
                               t.edge_endpoints == ttrace.edge_endpoints,
       walked = ε)
end
