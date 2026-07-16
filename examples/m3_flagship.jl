# M3 flagship chain: a double-well Schrödinger problem run through the whole exact-WKB
# pipeline — turning points → all-orders WKB → Voros symbol on the inner (instanton)
# cycle → Borel/Padé from Resurgence on the quantum series → Stokes graph at θ = 0 with
# the finite saddle edge → a rendered PNG.
#
# Run:  julia --project=test examples/m3_flagship.jl

using ExactWKB
import Resurgence
using CairoMakie

# (z²−1)² − 1/4 : symmetric double well, four simple turning points.
prob = SchrodingerProblem([3//4, 0, -2, 0, 1])
println(prob)

tps = turning_points(prob)
println("\nTurning points:")
foreach(t -> println("  ", t), tps)

# All-orders WKB, then the Voros symbol on the inner cycle (turning points 2 & 3).
w = wkb_expansion(prob; order = 6)
inner = encircling_contour(tps[2], tps[3])
vs = voros_symbol(w, inner)
println("\n", sprint((io, v) -> show(io, MIME"text/plain"(), v), vs))

# Hand the quantum series to Resurgence: Borel transform, then Padé (reduce = true to
# dodge the alternating-zero degeneracy) and its poles = Borel singularities.
Φ = quantum_series(vs)
B = Resurgence.borel(Φ)
P = Resurgence.pade(B; reduce = true)          # reduce = true dodges the ħ-parity degeneracy
println("\nResurgence handoff:")
println("  Borel series : ", B)
println("  Padé poles   : ", round.(Resurgence.poles(P); sigdigits = 4))

# Stokes graph at θ = 0 — the instanton saddle connects the inner turning points.
g = stokes_graph(prob; theta = 0.0)
println("\n", sprint((io, v) -> show(io, MIME"text/plain"(), v), g))

sad = saddles(prob)
println("\nBPS saddles (mass |Z|, critical phase θ_c):")
foreach(s -> println("  ", s), sad)

# Render the Stokes graph.
fig = plot_stokes_graph(g)
out = joinpath(@__DIR__, "m3_double_well_stokes.png")
CairoMakie.save(out, fig)
println("\nWrote ", out)
