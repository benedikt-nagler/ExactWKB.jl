# BPS spectra from maximal green sequences - the payoff of the cluster bridge.
#
# In physics language this file computes the BPS particle spectrum of a 4d N = 2
# quantum field theory: a degree-d polynomial potential realizes the Argyres–Douglas
# theory of type (A₁, A_{d−1}), whose BPS states are the saddle connections of the
# Stokes graph. A `BPSState` is a BPS hypermultiplet with charge γ in the chamber's
# charge lattice, N = 2 central charge Z_γ, mass |Z_γ| and BPS phase arg Z_γ; its
# `omega` is the Donaldson–Thomas invariant Ω(γ) (= 1 throughout finite type), and
# the phase-ordered spectrum is the Kontsevich–Soibelman wall-crossing product read
# off a maximal green sequence.
#
# The spectrum of a polynomial problem is enumerated cluster-algebraically: the seed
# of a chamber is `bridge_seed(charge_basis(...))`, and the physical maximal green
# sequence is constructed *greedily in phase order* (ledger item 6 of src/ddp.jl in
# its chamber-relative form: walls are crossed as θ decreases from the chamber's θ₀,
# so the MGS charge order is by increasing `mod(θ₀ − θ_c, π)`): among the currently
# green vertices, mutate the one whose c-vector charge `c` has the wall phase of
# `Z(c) = Σ c_j Z_j^phys` nearest below θ₀. This produces THE physical
# sequence in `#states` mutations - no enumeration of all maximal green sequences
# (A₄ already has hundreds). Two walls tied in phase are legal only when they
# commute (the current exchange-matrix entry vanishes - e.g. the two disjoint wells
# of a symmetric double well); a non-commuting tie is genuine marginal stability and
# throws `ChamberError`.
#
# Everything here is integer cluster machinery plus the numerical central charges;
# Voros symbols never enter. The definitive physical validation is built in
# (`verify = true`): the state multiset (mass, wall phase) must coincide with the
# confirmed `saddles(prob)` - the finite Stokes lines actually traced. Further
# oracles live in test_bps.jl: Ω ≡ 1 in finite type, Zamolodchikov periodicity
# h + 2, DT closure, chamber independence.
#
# ── the signed frame closed both original-bridge findings (2026-07-19) ─────────────────────
# The original bridge layer used the uniform all-decay frame and the *absolute* θ-decreasing order.
# Both are top-chamber specializations, and both failed elsewhere: cyclic chambers had
# no seed at all, and the quintic's real-axis fan chamber (θ ∈ (0, 0.103)) swept a
# cluster-consistent but unphysical 6-state sequence against 7 saddles. The signed
# frame (`src/signed_frame.jl`, ε from the keystone P = −εBε) and the chamber-relative
# order fix both, independently of each other.
#
# So the saddle cross-check `_matches_saddles` changed meaning rather than code: it
# was the *gate* that chose a chamber, and is now the *verification* of a frame that
# should be correct by construction. A failure is a bug, not a chamber property - and
# the default θ is simply the widest wall-free gap, no longer a search.

"""
    BPSState{F}

One BPS state of a [`BPSSpectrum`](@ref): its integer `charge` (a c-vector in the
chamber's [`ChargeBasis`](@ref)), its `central_charge` ``Z_γ = Σ c_j Z_j``, and its
DT invariant `omega` (``Ω(γ) = 1`` throughout finite type). Its mass is
[`mass`](@ref)`(s) = |Z_γ|` and its wall phase [`phase`](@ref)`(s) ∈ [0, π)`.
"""
struct BPSState{F}
    charge::Vector{Int}
    central_charge::Complex{F}
    omega::Int
end

charge(s::BPSState) = s.charge
central_charge(s::BPSState) = s.central_charge
omega(s::BPSState) = s.omega

"""
    mass(s::BPSState) -> Real

The BPS mass `|Z_γ|`.
"""
mass(s::BPSState) = abs(s.central_charge)

# wall phase in [0, π), robust at the real axis: a numerically-real Z is a wall at
# θ_c = 0, never at π − eps
function _wall_phase(Z::Complex{F}) where {F}
    tol = sqrt(eps(F)) * (1 + abs(Z))
    abs(imag(Z)) ≤ tol && return zero(F)
    mod(angle(Z), F(π))
end

"""
    phase(s::BPSState) -> Real

The wall phase ``θ_c = \\arg Z_γ \\pmod π ∈ [0, π)`` at which the state's finite
Stokes line appears.
"""
phase(s::BPSState) = _wall_phase(s.central_charge)

"""
    BPSSpectrum{F}

The BPS spectrum of a Schrödinger problem, enumerated from the chamber traced at
`theta0`: `states` in θ-decreasing wall-phase order (= the maximal green sequence
order, ledger item 6), the physical green `sequence` itself, and the chamber's
[`ChargeBasis`](@ref) `basis` in which the charges are expressed. Build with
[`bps_spectrum`](@ref).
"""
struct BPSSpectrum{F}
    states::Vector{BPSState{F}}
    sequence::Vector{Int}
    basis::ChargeBasis{F}
    theta0::F
end

"""
    n_states(sp::BPSSpectrum) -> Int

Number of BPS states in the spectrum.
"""
n_states(sp::BPSSpectrum) = length(sp.states)

"""
    charges(sp::BPSSpectrum) -> Vector{Vector{Int}}

The state charges in spectrum (θ-decreasing) order.
"""
charges(sp::BPSSpectrum) = [s.charge for s in sp.states]

central_charges(sp::BPSSpectrum) = [s.central_charge for s in sp.states]

# The greedy phase-ordered sweep on a chamber's charge basis. Returns the swept
# states and the green sequence; with `verify` the sweep is checked to be a genuine
# maximal green sequence with θ-decreasing phases and the Ω(γ) are attached.
#
# Ledger item 6 in its chamber-relative form: walls are crossed as θ DECREASES from
# the chamber's θ₀, so the sweep order is by increasing `mod(θ₀ − θ_c, π)`. As
# θ₀ → π⁻ this degenerates to the absolute θ-decreasing order in which item 6 was
# originally pinned - the top chamber is that special case, which is why the original
# absolute rule worked there and nowhere else.
function _sweep_chamber(cb::ChargeBasis{F}; verify::Bool) where {F}
    m = n_charges(cb)
    m == 0 && return BPSState{F}[], Int[]
    seed0 = bridge_seed(cb)
    Zb = physical_charges(cb)
    θ0 = F(cb.triangulation.theta)
    # distance travelled downwards in θ from θ₀ to a wall at phase p (period π)
    _rel(p) = mod(θ0 - p, F(π))
    tie_tol = 100 * sqrt(eps(F))
    seq = Int[]
    charges_swept = Vector{Vector{Int}}()
    zs = Complex{F}[]
    current = seed0
    while true
        greens = [k for k in 1:m if ClusterAlgebras.is_green(current, k)]
        isempty(greens) && break
        length(seq) ≤ 10_000 || throw(ChamberError(
            "the greedy phase sweep did not close after 10000 mutations"))
        cs = [ClusterAlgebras.c_vector(current, k) for k in greens]
        Zs = [sum(c[j] * Zb[j] for j in 1:m) for c in cs]
        ps = [_rel(_wall_phase(Z)) for Z in Zs]     # θ travelled down from θ₀
        pmin = minimum(ps)
        tied = [i for i in eachindex(greens) if ps[i] < pmin + tie_tol]
        if length(tied) > 1
            B = current.quiver.B
            all(B[greens[a], greens[b]] == 0 for a in tied, b in tied if a != b) ||
                throw(ChamberError(
                    "marginal stability: two non-commuting walls share the phase " *
                    "θ_c ≈ $(Float64(mod(θ0 - pmin, F(π)))) - perturb the potential " *
                    "or the energy"))
        end
        i = tied[argmin([greens[j] for j in tied])]    # deterministic among commuting ties
        push!(seq, greens[i])
        push!(charges_swept, cs[i])
        push!(zs, Zs[i])
        current = ClusterAlgebras.mutate(current, greens[i])
    end

    omegas = fill(1, length(seq))
    if verify
        ps = [_rel(_wall_phase(Z)) for Z in zs]
        all(ps[i + 1] ≥ ps[i] - tie_tol for i in 1:(length(ps) - 1)) ||
            throw(ChamberError(
                "the swept wall phases do not recede monotonically from θ₀ = " *
                "$(Float64(θ0)) - the greedy sequence is not the physical chamber " *
                "order (ledger item 6, chamber-relative form)"))
        cvs = try
            ClusterAlgebras.ordered_c_vectors(seed0, seq)
        catch err
            throw(ChamberError(
                "the greedy sweep is not a maximal green sequence: " *
                sprint(showerror, err)))
        end
        cvs == charges_swept || throw(ChamberError(
            "c-vector mismatch between the greedy sweep and its replay"))
        om = ClusterAlgebras.omega(ClusterAlgebras.quantum_dilog_word(seed0, seq))
        omegas = [get(om, c, 0) for c in charges_swept]
        all(>(0), omegas) || throw(ChamberError(
            "a swept charge is missing from the quantum dilogarithm word"))
    end
    [BPSState{F}(charges_swept[i], zs[i], omegas[i]) for i in eachindex(seq)], seq
end

# The definitive physical validation: the swept states must reproduce the confirmed
# saddles as a multiset of (mass, wall phase) - matched bijectively.
function _matches_saddles(states, sads; mrtol = 1e-5, ptol = 1e-5)
    length(states) == length(sads) || return false
    used = falses(length(sads))
    for s in states
        i = findfirst(eachindex(sads)) do k
            used[k] && return false
            Zs = central_charge(sads[k])
            abs(abs(Zs) - mass(s)) ≤ mrtol * abs(Zs) || return false
            d = abs(_wall_phase(complex(Zs)) - phase(s))
            min(d, π - d) ≤ ptol
        end
        i === nothing && return false
        used[i] = true
    end
    true
end

"""
    bps_spectrum(prob::SchrodingerProblem; theta = nothing, margin = nothing, n = 32,
                 rtol = nothing, verify = true, kwargs...) -> BPSSpectrum

The BPS spectrum of `prob`, computed cluster-algebraically through the cluster bridge:
trace the Stokes graph at `theta` (default: the midpoint of the largest wall-free
gap), build the [`ideal_triangulation`](@ref) → [`charge_basis`](@ref) →
[`bridge_seed`](@ref), then construct the physical maximal green sequence greedily
in chamber-relative wall-phase order. Singularity positions of the Borel-plane walls
are the returned central charges ``Z_γ``; the Stokes constants are the integer
``Ω(γ)·⟨γ,γ'⟩`` (verify numerically with [`verify_ddp`](@ref)). In physics terms
this is the BPS spectrum of the ``(A_1, A_{d-1})`` Argyres–Douglas theory at the
Coulomb-branch point selected by `Q`.

**Every chamber works.** With the signed frame (`src/signed_frame.jl`) the sweep is
correct in cyclic chambers and in the tree chambers that the original bridge layer
got wrong, so the default `theta` is simply the most numerically comfortable chamber
rather than a
search for one that happens to be right, and any explicit `theta` off a wall is
equally valid - the spectrum is chamber-independent.

`margin`/`n`/`rtol` pass to [`charge_basis`](@ref), remaining `kwargs` to
[`stokes_graph`](@ref)/[`saddles`](@ref). With `verify = true` (default) the wall
phases are checked to recede monotonically from `θ₀`, the greedy sequence is
re-validated as a maximal green sequence (`ClusterAlgebras.ordered_c_vectors`), the
DT invariants are attached from `ClusterAlgebras.omega`, and - the definitive
physical gate - the state multiset (mass, wall phase) must reproduce the confirmed
[`saddles`](@ref). A failure of that gate is now a bug rather than a property of the
chamber, and throws [`ChamberError`](@ref), as do a non-commuting phase tie (marginal
stability) and a sweep that fails to close.
"""
function bps_spectrum(prob::SchrodingerProblem; theta = nothing, margin = nothing,
                      n::Integer = 32, rtol = nothing, verify::Bool = true, kwargs...)
    F = _wkb_float(prob)
    sads = verify || theta === nothing ? saddles(prob; kwargs...) : Saddle{F}[]

    if theta !== nothing
        θ0 = F(theta)
        t = ideal_triangulation(stokes_graph(prob; theta = θ0, kwargs...))
        cb = charge_basis(prob, t; margin, n, rtol, verify)
        states, seq = _sweep_chamber(cb; verify)
        verify && !_matches_saddles(states, sads) && throw(ChamberError(
            "the chamber at θ = $(Float64(θ0)) swept $(length(states)) states that " *
            "do not reproduce the $(length(sads)) confirmed saddles in (mass, " *
            "phase). With the signed frame every chamber should be physical, so " *
            "this is an internal inconsistency, not a reason to pick another θ"))
        return BPSSpectrum{F}(states, seq, cb, θ0)
    end

    if isempty(sads) && isempty(saddle_candidates(prob))
        θ0 = F(π) / 2                            # e.g. Airy: no cycles, empty spectrum
        t = ideal_triangulation(stokes_graph(prob; theta = θ0, kwargs...))
        cb = charge_basis(prob, t; margin, n, rtol, verify)
        return BPSSpectrum{F}(BPSState{F}[], Int[], cb, θ0)
    end

    # Any chamber is physical, so take the widest wall-free gap - the most comfortable
    # θ to trace, not a search. (The original bridge layer looped over gaps looking for one whose
    # frame happened to be right; the signed frame makes that unnecessary.)
    phases = sort(unique([s.theta for s in saddle_candidates(prob)]))
    gaps = [(i == length(phases) ? phases[1] + F(π) - phases[end] :
             phases[i + 1] - phases[i], i) for i in eachindex(phases)]
    gap, i = maximum(gaps)
    gap > sqrt(eps(F)) || throw(ChamberError(
        "the walls of this potential leave no θ-gap wider than $(sqrt(eps(F))): its " *
        "saddles are numerically marginal, so no chamber can be traced generically"))
    θ0 = mod(phases[i] + gap / 2, F(π))
    t = ideal_triangulation(stokes_graph(prob; theta = θ0, kwargs...))
    cb = charge_basis(prob, t; margin, n, rtol, verify)
    states, seq = _sweep_chamber(cb; verify)
    verify && !_matches_saddles(states, sads) && throw(ChamberError(
        "the chamber at θ = $(Float64(θ0)) swept $(length(states)) states that do " *
        "not reproduce the $(length(sads)) confirmed saddles in (mass, phase) - " *
        "with the signed frame this is an internal inconsistency, not a chamber to " *
        "skip over"))
    BPSSpectrum{F}(states, seq, cb, θ0)
end
