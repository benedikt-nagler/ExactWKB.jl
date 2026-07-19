# Compact and `text/plain` display for the public types. Kept terse — a Stokes graph
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
    print(io, "WKBExpansion — all-orders WKB (Riccati) solution\n",
          "  order      : ", w.order, "  (S₋₁ … S", _sub(w.order), ")\n",
          "  arithmetic : ", w.arithmetic, "\n  potential  : ", w.prob)
end

# -- VorosSymbol -------------------------------------------------------------------

function Base.show(io::IO, vs::VorosSymbol)
    print(io, "VorosSymbol(v₋₁ = ", _disp(classical_period(vs)), ", ", quantum_series(vs), ")")
end

function Base.show(io::IO, ::MIME"text/plain", vs::VorosSymbol)
    print(io, "VorosSymbol — quantum period ∮ S dz\n",
          "  classical v₋₁ = ", _disp(classical_period(vs)), "\n",
          "  quantum series: ", quantum_series(vs))
end

# -- SpectralCycle -----------------------------------------------------------------

function Base.show(io::IO, c::SpectralCycle)
    print(io, "SpectralCycle(:", c.kind, ", ", _disp(real(c.left)), " ↔ ",
          _disp(real(c.right)), ")")
end

function Base.show(io::IO, ::MIME"text/plain", c::SpectralCycle)
    print(io, "SpectralCycle — ", c.kind === :well ? "classically allowed (well)" :
          "forbidden (barrier/tunneling)", " period cycle\n",
          "  turning points : ", _disp(c.left), " ↔ ", _disp(c.right), "\n",
          "  contour        : ", length(c.contour), "-gon")
end

# -- StokesLine --------------------------------------------------------------------

function Base.show(io::IO, l::StokesLine)
    dest = l.endpoint === :turning_point ? "→ TP $(l.target)" :
           l.endpoint === :infinity ? "→ ∞" : "→ incomplete"
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

# -- IdealTriangulation --------------------------------------------------------------

function Base.show(io::IO, t::IdealTriangulation)
    print(io, "IdealTriangulation(", t.n_marked, "-gon, ", length(t.triangles),
          " triangles, ", n_diagonals(t), " diagonals)")
end

function Base.show(io::IO, ::MIME"text/plain", t::IdealTriangulation)
    print(io, "IdealTriangulation of the ", t.n_marked, "-gon (dual Stokes regions)\n",
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

# -- ChargeBasis ---------------------------------------------------------------------

function Base.show(io::IO, cb::ChargeBasis)
    print(io, "ChargeBasis(", n_charges(cb), " cycles)")
end

function Base.show(io::IO, ::MIME"text/plain", cb::ChargeBasis)
    print(io, "ChargeBasis — signed frame on the triangulation diagonals")
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
    print(io, "BPSSpectrum — ", n_states(sp), " BPS states (chamber θ₀ = ",
          _disp(sp.theta0), ", θ-decreasing)")
    for s in sp.states
        print(io, "\n  ", rpad(string(s.charge), 4 + 3 * length(s.charge)),
              " |Z| = ", rpad(string(_disp(mass(s))), 10),
              " θ_c = ", rpad(string(_disp(phase(s))), 10), " Ω = ", s.omega)
    end
end

# subscript helper for small non-negative integers
function _sub(n::Integer)
    n < 0 && return "₋" * _sub(-n)
    digs = collect("₀₁₂₃₄₅₆₇₈₉")
    join(digs[d + 1] for d in reverse(digits(n)))
end
