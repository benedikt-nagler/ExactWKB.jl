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

# subscript helper for small non-negative integers
function _sub(n::Integer)
    n < 0 && return "₋" * _sub(-n)
    digs = collect("₀₁₂₃₄₅₆₇₈₉")
    join(digs[d + 1] for d in reverse(digits(n)))
end
