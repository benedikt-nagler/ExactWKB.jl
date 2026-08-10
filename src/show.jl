# Compact and `text/plain` display for the public types. Kept terse - a Stokes graph
# never dumps its traced polylines, a Voros symbol defers to its FormalSeries.

# round a real/complex to a few significant digits for readable display
_disp(x::Real) = round(float(x); sigdigits = 5)
_disp(z::Complex) = complex(round(float(real(z)); sigdigits = 5),
                            round(float(imag(z)); sigdigits = 5))

# -- SchrodingerProblem ------------------------------------------------------------

function Base.show(io::IO, prob::SchrodingerProblem)
    print(io, "SchrodingerProblem(Q deg ", degree(prob), ", E = ", _disp(energy(prob)), ")")
end

function Base.show(io::IO, ::MIME"text/plain", prob::SchrodingerProblem)
    q = q_coefficients(prob)
    v = variable(prob)
    terms = String[]
    for (k, c) in pairs(q)
        iszero(c) && continue
        pw = k == 1 ? "" : k == 2 ? string(v) : "$v^$(k-1)"
        push!(terms, isempty(pw) ? string(_disp(c)) : "$(_disp(c))*$pw")
    end
    body = isempty(terms) ? "0" : replace(join(terms, " + "), "+ -" => "- ")
    print(io, "SchrodingerProblem (ħ²ψ″ = Q ψ, Q = V − E)\n",
          "  Q(", v, ") = ", body, "\n  E = ", _disp(energy(prob)))
end

# -- RationalProblem ---------------------------------------------------------------

function Base.show(io::IO, prob::RationalProblem)
    print(io, "RationalProblem(Q ~ z^", degree(prob), " at ∞, ",
          n_finite_poles(prob), " finite pole",
          n_finite_poles(prob) == 1 ? "" : "s", ")")
end

function Base.show(io::IO, ::MIME"text/plain", prob::RationalProblem)
    v = variable(prob)
    terms = String[]
    for (k, c) in pairs(q_numerator(prob))
        iszero(c) && continue
        pw = k == 1 ? "" : k == 2 ? string(v) : "$v^$(k-1)"
        push!(terms, isempty(pw) ? string(_disp(c)) : "$(_disp(c))*$pw")
    end
    num = isempty(terms) ? "0" : replace(join(terms, " + "), "+ -" => "- ")
    den = join(["($v − $(_disp(p)))^$m"
                for (p, m) in zip(poles(prob), pole_orders(prob))], " ")
    print(io, "RationalProblem (ħ²ψ″ = Q ψ)\n",
          "  Q(", v, ") = [", num, "] / ", isempty(den) ? "1" : den, "\n",
          "  poles     : ")
    if n_finite_poles(prob) == 0
        print(io, "∞ only (order ", degree(prob) + 4, ")")
    else
        print(io, join(["$(_disp(p)) (order $m)"
                        for (p, m) in zip(poles(prob), pole_orders(prob))], ", "),
              ", ∞ (order ", degree(prob) + 4, ")")
    end
    print(io, "\n  surface   : ", n_finite_poles(prob) + 1, " boundary circle",
          n_finite_poles(prob) == 0 ? "" : "s", ", ",
          sum(m - 2 for m in pole_orders(prob); init = 0) + degree(prob) + 2,
          " marked points")
end

# -- TurningPoint ------------------------------------------------------------------

function Base.show(io::IO, tp::TurningPoint)
    kind = is_simple(tp) ? "simple" : "order $(order(tp))"
    print(io, "TurningPoint(", _disp(location(tp)), ", ", kind, ")")
end

# -- WKBExpansion ------------------------------------------------------------------

function Base.show(io::IO, w::WKBExpansion)
    print(io, "WKBExpansion(order ", w.order, ", ", w.arithmetic, ")")
end

function Base.show(io::IO, ::MIME"text/plain", w::WKBExpansion)
    print(io, "WKBExpansion - all-orders WKB (Riccati) solution\n",
          "  order      : ", w.order, "  (S₋₁ … S", _sub(w.order), ")\n",
          "  arithmetic : ", w.arithmetic, "\n  potential  : ", w.prob)
end

# -- VorosSymbol -------------------------------------------------------------------

function Base.show(io::IO, vs::VorosSymbol)
    print(io, "VorosSymbol(v₋₁ = ", _disp(classical_period(vs)), ", ", quantum_series(vs), ")")
end

function Base.show(io::IO, ::MIME"text/plain", vs::VorosSymbol)
    print(io, "VorosSymbol - quantum period ∮ S dz\n",
          "  classical v₋₁ = ", _disp(classical_period(vs)), "\n",
          "  quantum series: ", quantum_series(vs))
end

# -- SpectralCycle -----------------------------------------------------------------

function Base.show(io::IO, c::SpectralCycle)
    print(io, "SpectralCycle(:", c.kind, ", ", _disp(real(c.left)), " ↔ ",
          _disp(real(c.right)), ")")
end

function Base.show(io::IO, ::MIME"text/plain", c::SpectralCycle)
    print(io, "SpectralCycle - ", c.kind === :well ? "classically allowed (well)" :
          "forbidden (barrier/tunneling)", " period cycle\n",
          "  turning points : ", _disp(c.left), " ↔ ", _disp(c.right), "\n",
          "  contour        : ", length(c.contour), "-gon")
end

# -- StokesLine --------------------------------------------------------------------

function Base.show(io::IO, l::StokesLine)
    dest = l.endpoint === :turning_point ? "→ TP $(l.target)" :
           l.endpoint === :infinity ? "→ ∞" :
           l.endpoint === :pole ? "→ pole $(l.target)" : "→ incomplete"
    print(io, "StokesLine(TP ", l.source, ", ray ", l.direction, " ", dest,
          ", mass ", _disp(l.mass), ")")
end

# -- StokesGraph -------------------------------------------------------------------

function Base.show(io::IO, g::StokesGraph)
    print(io, "StokesGraph(θ = ", _disp(g.theta), ", ", length(g.turning_points),
          " TPs, ", length(edges(g)), " saddles, ", n_infinite_lines(g), " rays→∞)")
end

function Base.show(io::IO, ::MIME"text/plain", g::StokesGraph)
    print(io, "StokesGraph at θ = ", _disp(g.theta), "\n",
          "  turning points : ", length(g.turning_points), "\n",
          "  Stokes lines   : ", length(g.lines),
          " (", n_infinite_lines(g), " to ∞, ", count(is_finite_line, g.lines), " saddle-traced)\n",
          "  saddle edges   : ", isempty(edges(g)) ? "none" : join(edges(g), ", "), "\n",
          "  signature      : ", topology_signature(g))
end

# -- Saddle ------------------------------------------------------------------------

function Base.show(io::IO, s::Saddle)
    print(io, "Saddle(", s.pair[1], "–", s.pair[2], ", |Z| = ", _disp(mass(s)),
          ", θ_c = ", _disp(s.theta), ")")
end

function Base.show(io::IO, ::MIME"text/plain", s::Saddle)
    print(io, "Saddle (BPS state) between turning points ", s.pair[1], " and ", s.pair[2], "\n",
          "  central charge Z = ", _disp(s.central_charge), "\n",
          "  mass |Z|         = ", _disp(mass(s)), "\n",
          "  critical phase θ = ", _disp(s.theta), " (mod π)")
end

# -- RingDomainWall ------------------------------------------------------------------

function Base.show(io::IO, w::RingDomainWall)
    print(io, "RingDomainWall(pole ", w.pole, ", |Z| = ", _disp(mass(w)),
          ", θ_c = ", _disp(w.theta), ")")
end

function Base.show(io::IO, ::MIME"text/plain", w::RingDomainWall)
    print(io, "RingDomainWall at the puncture (double pole ", w.pole, ")\n",
          "  residue charge Z = ", _disp(w.central_charge), "  (= ∮√Q)\n",
          "  mass |Z|         = ", _disp(mass(w)), "\n",
          "  critical phase θ = ", _disp(w.theta), " (mod π): there the trajectories\n",
          "                     close into a ring domain instead of spiralling in")
end

# -- IdealTriangulation --------------------------------------------------------------

# The surface, named the way a reader would name it. A disk is "the m-gon" (the only
# case before poles existed), and the general one has to say what it is: `n_boundaries`
# counts singularities, of which `n_punctures` are punctures rather than circles.
function _surface_desc(nm, nsing, npunc)
    circles, m = nsing - npunc, nm - npunc
    npunc == 0 && circles == 1 && return "$(m)-gon"
    npunc == 0 && circles == 2 && return "annulus ($m marked points)"
    npunc == 1 && circles == 1 && return "once-punctured $(m)-gon"
    "sphere with $circles boundary circle$(circles == 1 ? "" : "s") " *
        "($m marked point$(m == 1 ? "" : "s")) and $npunc puncture" *
        (npunc == 1 ? "" : "s")
end
_surface_desc(t) = _surface_desc(t.n_marked, n_boundaries(t), n_punctures(t))

function Base.show(io::IO, t::IdealTriangulation)
    print(io, "IdealTriangulation(", _surface_desc(t), ", ", length(t.triangles),
          " triangles, ", n_diagonals(t), " diagonals)")
end

function Base.show(io::IO, ::MIME"text/plain", t::IdealTriangulation)
    print(io, "IdealTriangulation of the ", _surface_desc(t), " (dual Stokes regions)\n",
          "  triangles : ", length(t.triangles), " (one per turning point)\n",
          "  diagonals : ")
    if n_diagonals(t) == 0
        print(io, "none")
    else
        parts = ["γ$(t.diagonal_tp_pair[e]) ↔ $(t.edge_endpoints[e])"
                 for e in 1:n_diagonals(t)]
        print(io, join(parts, ", "))
    end
end

# -- PolygonDecomposition ------------------------------------------------------------

function Base.show(io::IO, d::PolygonDecomposition)
    print(io, "PolygonDecomposition(", _surface_desc(d), ", ", n_cells(d), " cells ",
          Tuple(cell_sizes(d)), ", ", n_diagonals(d), " diagonals)")
end

function Base.show(io::IO, ::MIME"text/plain", d::PolygonDecomposition)
    print(io, "PolygonDecomposition of the ", _surface_desc(d),
          " (dual Stokes regions)\n",
          "  cells     : ", join(string.(cell_sizes(d)), "-gon, "), "-gon",
          " (one per turning point)\n",
          "  diagonals : ", n_diagonals(d), "\n",
          "  refines to: ", n_refinements(d), " ideal triangulation",
          n_refinements(d) == 1 ? "" : "s")
end

# -- ChargeBasis ---------------------------------------------------------------------

function Base.show(io::IO, cb::ChargeBasis)
    print(io, "ChargeBasis(", n_charges(cb), " cycles)")
end

function Base.show(io::IO, ::MIME"text/plain", cb::ChargeBasis)
    print(io, "ChargeBasis - signed frame on the triangulation diagonals")
    Zp = physical_charges(cb)
    for e in 1:n_charges(cb)
        print(io, "\n  γ", cb.triangulation.diagonal_tp_pair[e],
              " : ε = ", cb.signs[e] > 0 ? "+1" : "−1",
              ", Z = ", _disp(Zp[e]))
        cb.signs[e] < 0 && print(io, "  (decay rep. ", _disp(cb.central_charges[e]), ")")
    end
    n_charges(cb) > 0 && print(io, "\n  signed pairing = ", signed_pairing(cb))
end

# -- BPSState / BPSSpectrum ----------------------------------------------------------

function Base.show(io::IO, s::BPSState)
    print(io, "BPSState(", s.charge, ", Z = ", _disp(s.central_charge),
          ", Ω = ", s.omega, ")")
end

function Base.show(io::IO, sp::BPSSpectrum)
    print(io, "BPSSpectrum(", n_states(sp), " states, θ₀ = ", _disp(sp.theta0), ")")
end

function Base.show(io::IO, ::MIME"text/plain", sp::BPSSpectrum)
    print(io, "BPSSpectrum - ", n_states(sp), " BPS states (chamber θ₀ = ",
          _disp(sp.theta0), ", θ-decreasing)")
    for s in sp.states
        print(io, "\n  ", rpad(string(s.charge), 4 + 3 * length(s.charge)),
              " |Z| = ", rpad(string(_disp(mass(s))), 10),
              " θ_c = ", rpad(string(_disp(phase(s))), 10), " Ω = ", s.omega)
    end
end

# -- TBASystem / TBASolution ---------------------------------------------------------

function Base.show(io::IO, sys::TBASystem)
    print(io, "TBASystem(", n_states(sys), " states)")
end

function Base.show(io::IO, sol::TBASolution)
    print(io, "TBASolution(", n_states(sol.system), " states, ",
          length(sol.grid), "-point grid, ", sol.iterations,
          " sweeps, residual ", _disp(sol.residual), ")")
end

# -- SeibergWittenSU2 --------------------------------------------------------------

function Base.show(io::IO, sw::SeibergWittenSU2)
    print(io, "SeibergWittenSU2(Λ = ", _disp(dynamical_scale(sw)), ")")
end

function Base.show(io::IO, ::MIME"text/plain", sw::SeibergWittenSU2)
    s = sw_singularities(sw)
    print(io, "SeibergWittenSU2 - pure SU(2) N=2 Seiberg–Witten geometry\n",
          "  Λ            : ", _disp(dynamical_scale(sw)), "\n",
          "  singularities: monopole u = ", _disp(s.monopole),
          ", dyon u = ", _disp(s.dyon))
end

# subscript helper for small non-negative integers
function _sub(n::Integer)
    n < 0 && return "₋" * _sub(-n)
    digs = collect("₀₁₂₃₄₅₆₇₈₉")
    join(digs[d + 1] for d in reverse(digits(n)))
end
