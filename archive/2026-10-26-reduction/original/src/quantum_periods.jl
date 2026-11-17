# Energy-parametrized quantum periods: the bridge from a fixed problem to
# its spectral family. At an energy E the real simple turning points of a real
# polynomial potential cut the real axis into classically allowed (Q < 0, "well") and
# forbidden (Q > 0, "barrier") intervals; each well, and each barrier strictly between
# two wells, carries a period cycle. The cycle's quantum period is the Voros symbol of
# its encircling contour - the full Stokes-graph pipeline re-run at E (`with_energy` re-solve,
# no new series type).
#
# ── Spectral orientation ledger (pinned by the harmonic oracle in test_quantum_periods) ──
# Cycles are built as `encircling_contour(left, right)` with the turning points in
# ascending real order. With that orientation:
#   * a WELL cycle has classical period v₋₁ = −i·J(E), J(E) = 2∫√(E−V) dz > 0 the
#     classical action (harmonic: v₋₁ = −iπE, so |v₋₁|/(2πħ) = n + ½ quantizes);
#   * a BARRIER cycle has real negative v₋₁ = −2∫√(V−E) dz (minus twice the
#     tunneling/instanton action - the decaying orientation of the DDP layer).
# Contours are never reoriented downstream (DDP ledger item 4).

"""
    SpectralCycle{F}

A period cycle of a [`SchrodingerProblem`](@ref) family at a fixed energy: `kind`
is `:well` (classically allowed interval, ``Q < 0``) or `:barrier` (forbidden
interval between two wells, ``Q > 0``), `left`/`right` are the bounding real simple
turning points (ascending), and `contour` is the closed encircling contour the
quantum period integrates along. Build with [`spectral_cycles`](@ref); feed to
[`quantum_period`](@ref).
"""
struct SpectralCycle{F}
    kind::Symbol
    left::Complex{F}
    right::Complex{F}
    contour::Vector{Complex{F}}
end

"""
    kind(c::SpectralCycle) -> Symbol

`:well` or `:barrier`.
"""
kind(c::SpectralCycle) = c.kind

"""
    endpoints(c::SpectralCycle) -> Tuple{Complex,Complex}

The bounding turning points `(left, right)` in ascending real order.
"""
endpoints(c::SpectralCycle) = (c.left, c.right)

Base.:(==)(a::SpectralCycle, b::SpectralCycle) =
    a.kind == b.kind && a.left == b.left && a.right == b.right && a.contour == b.contour

# real polynomial potential guard (scope: real potentials, real spectra)
function _require_real_potential(prob::SchrodingerProblem, what)
    all(x -> isreal(x) || iszero(imag(x)), prob.v_coeffs) && isreal(prob.energy) ||
        throw(Resurgence.InvalidArgument(
            "$what needs a real polynomial potential and real energy (real-spectrum scope)"))
end

# distance from point p to the closed segment [a, b] (all complex)
function _segment_distance(p, a, b)
    d = b - a
    L2 = abs2(d)
    iszero(L2) && return abs(p - a)
    t = clamp(real(conj(d) * (p - a)) / L2, 0, 1)
    abs(p - (a + t * d))
end

"""
    spectral_cycles(prob::SchrodingerProblem, E; margin = nothing, n = 32,
                    real_tol = nothing, coalescence_tol = nothing)
        -> Vector{SpectralCycle}

The period cycles of the family `Q = V − E` at energy `E`, for a **real polynomial
potential** (the real-spectrum scope): the real simple turning points are sorted; classically
allowed intervals (`Q < 0`) give `:well` cycles and forbidden intervals strictly
between two wells give `:barrier` cycles, in ascending position order. Each cycle
carries an [`encircling_contour`](@ref) sized to enclose exactly its two turning
points (`margin` overrides the automatic padding, `n` the polygon resolution).

Throws [`CoalescentTurningPoints`](@ref) when `E` is too close to a critical value
of `V` - a multiple turning point, or two real turning points closer than
`coalescence_tol` (default `eps^{1/4}` of the working float type, scaled).
`real_tol` sets the reality tolerance for turning points (default `√eps`, scaled).
Returns an empty vector when no interval between real turning points is allowed
(e.g. `E` below the potential minimum).
"""
function spectral_cycles(prob::SchrodingerProblem, E; margin = nothing,
                         n::Integer = 32, real_tol = nothing, coalescence_tol = nothing)
    _require_real_potential(prob, "spectral_cycles")
    probE = with_energy(prob, E)
    tps = turning_points(probE)
    F = _wkb_float(probE)
    scale = 1 + maximum(abs ∘ location, tps; init = zero(F))
    for tp in tps
        is_simple(tp) ||
            throw(CoalescentTurningPoints(E, location(tp), zero(F)))
    end
    rtol = real_tol === nothing ? sqrt(eps(F)) * scale : F(real_tol)
    real_tps = [tp for tp in tps if abs(imag(location(tp))) ≤ rtol]
    sort!(real_tps; by = real ∘ location)
    ctol = coalescence_tol === nothing ? eps(F)^(1 // 4) * scale : F(coalescence_tol)
    for k in 1:(length(real_tps) - 1)
        a, b = location(real_tps[k]), location(real_tps[k + 1])
        sep = abs(b - a)
        sep ≤ ctol && throw(CoalescentTurningPoints(E, (a + b) / 2, sep))
    end
    length(real_tps) < 2 && return SpectralCycle{F}[]

    # classify the consecutive intervals by the sign of Q at the midpoint
    kinds = Symbol[]
    for k in 1:(length(real_tps) - 1)
        a, b = location(real_tps[k]), location(real_tps[k + 1])
        qm = real(probE((a + b) / 2))
        push!(kinds, qm < 0 ? :well : :barrier)
    end

    others(k) = [location(tp) for tp in tps
                 if tp !== real_tps[k] && tp !== real_tps[k + 1]]
    cycles = SpectralCycle{F}[]
    for k in eachindex(kinds)
        kinds[k] === :barrier &&
            !(k > 1 && kinds[k - 1] === :well && k < length(kinds) &&
              kinds[k + 1] === :well) && continue
        # Strip the O(eps)-imaginary noise of the (certified real) turning points:
        # the contour's first vertex then lies exactly on the real axis, where the
        # principal √Q branch is deterministic. With a residual ±1e-150im jitter,
        # Q(z₀) < 0 lands on either side of the √ branch cut and the whole cycle's
        # period flips sign nondeterministically across an energy stencil.
        a = Complex{F}(real(location(real_tps[k])))
        b = Complex{F}(real(location(real_tps[k + 1])))
        pad = margin
        if pad === nothing
            dmin = minimum((_segment_distance(z, a, b) for z in others(k));
                           init = F(Inf))
            pad = min(abs(b - a) / 2, isfinite(dmin) ? F(9 // 20) * dmin : F(Inf))
        end
        ct = encircling_contour(a, b; margin = pad, n)
        # Orientation pins (see the ledger above): with the deterministic principal
        # branch at the exact-real start vertex, the counterclockwise contour gives
        # a well cycle v₋₁ = −iJ (the harmonic oracle) but a barrier cycle +2S_I;
        # barrier cycles are therefore traversed clockwise (same start vertex), so
        # their classical period is −2S_I < 0 - the decaying orientation.
        kinds[k] === :barrier && (ct = vcat(ct[1:1], reverse(ct[2:end])))
        push!(cycles, SpectralCycle{F}(kinds[k], a, b, ct))
    end
    cycles
end

# ── the same two cycles, continued through the barrier top ──────────────────────
#
# `spectral_cycles` reads the classical regions off the real axis, so it stops existing
# where they do: at the barrier top the two inner turning points merge and above it they
# leave the real axis, taking the barrier cycle's classical *meaning* with them. The
# cycles themselves are unbothered - they are cycles on `w² = Q`, and `E` moves their
# branch points around without destroying them. `uniform_cycles` builds that
# continuation, which is what the uniform quantization condition needs.
#
# ── Continuation ledger (pinned in test_quantization.jl) ────────────────────────
#
# 1. WHICH INNER TURNING POINT the well cycle runs to. Below the top it is the one
#    adjacent to the left outer turning point (`b` in `a < b < c < d`). Above the top
#    the pair has become `∓iδ`, and the continuation is the one with POSITIVE imaginary
#    part: only that choice gives `φ = Φ₀ − iπε/2` with `Φ₀ = (1/ħ)∫_a^0 √(E−V) > 0`,
#    the real half-well phase. The other choice differs by the sign of `Re v₋₁` and
#    would put the barrier term in with the wrong sign. PINNED by
#    `−Im v₋₁/(2ħ) = ∫_a^0 √(E−V)/ħ` continuing smoothly across `E_top`.
# 2. THE BARRIER MARGIN IS NOT CAPPED BY THE SEPARATION. `spectral_cycles` pads a cycle
#    by `min(|b−a|/2, 0.45·d_other)`; here the first term is dropped for the inner pair,
#    because it is exactly the term that shrinks to zero as the pair merges. The
#    contour stays a fixed size while the cycle it carries vanishes - which is the whole
#    point, since `v₋₁ → 0` continuously and only the contour would have degenerated.
# 3. THE SIGN OF ε IS GEOMETRIC, NOT ANALYTIC. `ε = ∓|Re log V_A|/2π`, positive when the
#    inner pair is real (a genuine barrier, `V_A = e^{−2S_I/ħ} < 1`) and negative when it
#    is a conjugate pair (over-barrier). The principal `√Q` branch at a contour's first
#    vertex flips when the pair rotates from real to imaginary, so reading the sign off
#    the integral instead would flip it exactly at the top - the one place continuity
#    matters. PINNED by `ε → 0` continuously from both sides.

"""
    uniform_cycles(prob::SchrodingerProblem, E; margin = nothing, n = 32,
                   real_tol = nothing) -> NamedTuple

The well and vanishing cycles of a symmetric double well at energy `E`, **continued
through the barrier top**: a named tuple `(well, barrier, below)` of two
[`SpectralCycle`](@ref)s and the flag `below = E < E_top`.

Below the barrier top these are the left well and the barrier of
[`spectral_cycles`](@ref). Above it the inner turning points have left the real axis
and neither classical region exists any more, but both cycles do: `well` runs from the
left outer turning point to the inner one with positive imaginary part, and `barrier`
encircles the (now conjugate) inner pair - the cycle that *vanishes* at the top. This
is the input of the uniform quantization condition
([`quantization_condition`](@ref) with `uniform = true`).

Requires a four-turning-point symmetric double-well family, all turning points simple
(so not the barrier top itself, nor the well bottom, where a pair has merged and the
well cycle would enclose one branch point).
"""
function uniform_cycles(prob::SchrodingerProblem, E; margin = nothing, n::Integer = 32,
                        real_tol = nothing)
    _require_real_potential(prob, "uniform_cycles")
    probE = with_energy(prob, E)
    tps = turning_points(probE)
    F = _wkb_float(probE)
    # merged pairs first: at the barrier top (and at the well bottom) the count is
    # short precisely BECAUSE a pair has merged, and that is the informative error
    for tp in tps
        is_simple(tp) || throw(CoalescentTurningPoints(E, location(tp), zero(F)))
    end
    length(tps) == 4 || throw(QuantizationError(
        "uniform_cycles needs a four-turning-point double-well family, found " *
        "$(length(tps)) turning points at E = $E"))
    zs = [location(tp) for tp in tps]
    scale = 1 + maximum(abs, zs)
    rtol = real_tol === nothing ? sqrt(eps(F)) * scale : F(real_tol)
    reals = [z for z in zs if abs(imag(z)) ≤ rtol]
    length(reals) in (2, 4) || throw(QuantizationError(
        "uniform_cycles found $(length(reals)) real turning points at E = $E; a " *
        "symmetric double well has four below the barrier top and two above it"))
    below = length(reals) == 4
    sort!(reals; by = real)
    a = Complex{F}(real(first(reals)))          # left outer turning point
    d = Complex{F}(real(last(reals)))
    inner = if below
        [Complex{F}(real(z)) for z in reals[2:3]]
    else
        cs = [z for z in zs if abs(imag(z)) > rtol]
        length(cs) == 2 || throw(QuantizationError(
            "uniform_cycles: $(length(cs)) non-real turning points at E = $E, expected " *
            "a conjugate pair"))
        sort(cs; by = imag)
    end
    # ledger item 1: the inner turning point the well cycle runs to
    b = below ? inner[1] : inner[2]
    c = below ? inner[2] : inner[1]
    # `a` and `b` are rebuilt (the real ones with their O(eps) imaginary noise stripped,
    # as in `spectral_cycles`), so the other turning points are selected by distance,
    # never by identity - the latter silently keeps `a` and `b` themselves and collapses
    # the padding to zero.
    ztol = sqrt(eps(F)) * scale
    others = [z for z in zs if abs(z - a) > ztol && abs(z - b) > ztol]

    padW = margin
    if padW === nothing
        dmin = minimum((_segment_distance(z, a, b) for z in others); init = F(Inf))
        padW = min(abs(b - a) / 2, isfinite(dmin) ? F(9 // 20) * dmin : F(Inf))
    end
    ctW = encircling_contour(a, b; margin = padW, n)

    # ledger item 2: no `|c − b|/2` cap - that is the term that would vanish with the
    # cycle. The remaining cap is the distance from the OUTER turning points to the
    # pair's segment, which is what keeps the ellipse (semi-major `|c−b|/2 + pad`) from
    # swallowing them; measuring to the pair's midpoint instead is not enough, and
    # silently encloses `a` and `d` at low energy, where the pair is wide.
    padB = margin
    if padB === nothing
        padB = F(9 // 20) * min(_segment_distance(a, b, c), _segment_distance(d, b, c))
    end
    ctB = encircling_contour(inner[1], inner[2]; margin = padB, n)

    (well = SpectralCycle{F}(:well, a, b, ctW),
     barrier = SpectralCycle{F}(:barrier, Complex{F}(inner[1]), Complex{F}(inner[2]),
                                ctB),
     below = below)
end

"""
    quantum_period(prob::SchrodingerProblem, cycle::SpectralCycle, E;
                   order = 12, arithmetic = :auto) -> VorosSymbol
    quantum_period(prob::SchrodingerProblem, kind::Symbol, E; kwargs...) -> VorosSymbol

The quantum period (Voros symbol) of `cycle` at energy `E`: the full WKB pipeline
re-run on `with_energy(prob, E)` - [`wkb_expansion`](@ref) to `order`, then
[`voros_symbol`](@ref) along the cycle's contour. The `kind` form (`:well` or
`:barrier`) resolves the cycle via [`spectral_cycles`](@ref) and requires it to be
unique of its kind (a symmetric double well has two `:well` cycles - pass the
`SpectralCycle` explicitly there).
"""
function quantum_period(prob::SchrodingerProblem, cycle::SpectralCycle, E;
                        order::Integer = 12, arithmetic::Symbol = :auto, rtol = nothing)
    w = wkb_expansion(with_energy(prob, E); order, arithmetic)
    voros_symbol(w, cycle.contour; rtol)
end

function quantum_period(prob::SchrodingerProblem, cyclekind::Symbol, E;
                        order::Integer = 12, arithmetic::Symbol = :auto,
                        rtol = nothing, kwargs...)
    cycles = spectral_cycles(prob, E; kwargs...)
    matches = [c for c in cycles if c.kind === cyclekind]
    length(matches) == 1 || throw(Resurgence.InvalidArgument(
        "quantum_period(prob, :$cyclekind, E) needs exactly one :$cyclekind cycle " *
        "at E = $E, found $(length(matches)); pass the SpectralCycle explicitly"))
    quantum_period(prob, matches[1], E; order, arithmetic, rtol)
end
