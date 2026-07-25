# ExactWKB.jl

[![CI](https://github.com/benedikt-nagler/ExactWKB.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benedikt-nagler/ExactWKB.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Exact WKB analysis of one-dimensional Schrödinger equations with polynomial potentials —
equivalently, of second-order linear ODEs with a small parameter: turning points, the
all-orders WKB expansion, quantum periods and Voros symbols, Stokes graphs, exact
quantization conditions, BPS spectra and wall-crossing, and the dictionary to cluster
algebras.

The object is the Schrödinger equation in its complex-analytic form

$$\hbar^2 \psi''(z) = Q(z)\,\psi(z), \qquad Q(z) = V(z) - E,$$

with $Q$ a polynomial and $z$ complex — which is also the general second-order linear ODE
with a large parameter, once it is put in normal form. The WKB ansatz $\psi = \exp\left(\hbar^{-1}\int S\right)$
turns it into a Riccati equation solved by a power series in $\hbar$ — a series that diverges
for every $\hbar \neq 0$. **Exact WKB** is the theory that makes it usable anyway: Borel
summed, it produces genuine solutions rather than an asymptotic approximation, and the
combinatorial object that records which solution you are on — the *Stokes graph*, the level
set $\mathrm{Im}(e^{-i\theta}\int\sqrt{Q}) = 0$ — changes discontinuously as the summation
direction $\theta$ rotates.

What that buys you, concretely:

- **spectra beyond leading order.** `wkb_eigenvalue` solves an all-orders Bohr–Sommerfeld
  condition rather than the leading-order one, and the double-well machinery reaches
  exponentially small quantities — the level splitting $\sim e^{-A/\hbar}$ — that no order of
  ordinary perturbation theory can see.
- **the Stokes structure of the equation.** `stokes_graph` gives the Stokes regions and how
  they reconnect as $\theta$ varies; `saddles` finds the special directions where the
  topology jumps. This is the classical connection-problem use of WKB for any second-order
  linear ODE with a large parameter, made exact.
- **periods as functions, not just numbers.** `voros_symbol` keeps a period as an
  $\hbar$-series, `voros_value` Borel sums it on either side of a Stokes ray, and
  `ddp_transform` applies the Delabaere–Dillinger–Pham jump between the two.
- **a physics layer on top:** BPS spectra and their central charges, the
  Gaiotto–Moore–Neitzke TBA, and the pure $SU(2)$ Seiberg–Witten geometry.

The Stokes graph is also where cluster algebras enter. Iwaki and Nakanishi observed
that a generic Stokes graph is dual to an ideal triangulation of a polygon, that the jump of
the Voros symbols across a wall is exactly a cluster $y$-mutation, and that a full circuit of
the $\theta$-plane is a maximal green sequence. The package computes both sides and checks
that they agree — the combinatorics is exact, so it certifies the numerics.

**Scope.** $Q$ must be polynomial. That covers the harmonic, quartic, double-well and cubic
oscillators and everything else usually used to test the theory; the Seiberg–Witten layer
adds the Mathieu problem through its own period module. Exact `Rational` input flows through
to `BigFloat` at whatever precision the caller has set, so the same problem runs cheaply in
`Float64` or carefully in high precision.

`ExactWKB.jl` is the bridge package of a small ecosystem: divergent-series machinery comes
from [Resurgence.jl](https://github.com/benedikt-nagler/Resurgence.jl), cluster combinatorics
from [ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl).

## Installation

None of the three packages is registered in General yet. From the Julia REPL:

```julia
pkg> add https://github.com/benedikt-nagler/ClusterAlgebras.jl
pkg> add https://github.com/benedikt-nagler/Resurgence.jl
pkg> add https://github.com/benedikt-nagler/ExactWKB.jl
```

Requires Julia 1.10 or later.

## Quick start

```julia
using ExactWKB

# Double well, Q = (z² − 1)² − 1/4, coefficients in ascending order
prob = SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0])

tps = turning_points(prob)           # four simple turning points
w   = wkb_expansion(prob; order = 6) # Riccati recursion for S_odd
vs  = voros_symbol(w, encircling_contour(tps[2], tps[3]))

classical_period(vs)                 # the ħ⁻¹ term, ∮√Q ≈ 1.8399
g = stokes_graph(prob; theta = 0.3)  # traced Stokes lines
saddles(prob)                        # 3 saddle connections, with central charges
```

Exact quantization gives energy levels rather than just series:

```julia
harm = SchrodingerProblem([0.0, 0.0, 1.0])   # V = z², so Eₙ = 2ħ(n + ½)
wkb_eigenvalue(harm, 0, 0.1; order = 6)      # 0.1
wkb_eigenvalue(harm, 3, 0.1; order = 6)      # 0.7
```

For the harmonic oscillator the WKB series truncates, so this is exact rather than
asymptotic — which makes it a good convention check.

And the cluster bridge, on the cubic:

```julia
cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])  # Q = z³ − z
g     = stokes_graph(cubic; theta = 0.3)
t     = ideal_triangulation(g)       # dual triangulation of the (d+2)-gon
q     = triangulation_quiver(t)      # a ClusterAlgebras.Quiver — here type A₂
cb    = charge_basis(cubic, t; margin = 0.4)   # charges, central charges, pairing
signed_pairing(cb) == -q.B           # the keystone identity

sp = bps_spectrum(cubic; theta = 0.3, margin = 0.4)
charges(sp)                          # [[1, 0], [0, 1]] — BPS charges in phase order
sp.sequence                          # [1, 2] — the same data as a green sequence
```

The BPS spectrum is enough to reconstruct the periods a second time, with no series and no
Padé anywhere in the pipeline:

```julia
sol = solve_tba(sp)                  # conformal-limit GMN integral equations
```

## Examples

[`examples/ddp_demo.jl`](examples/ddp_demo.jl) runs the Voros-jump layer end to end: Borel
summing a Voros symbol on both sides of a Stokes ray, applying the
Delabaere–Dillinger–Pham jump, and checking it against a cluster $y$-mutation.

## Functionality

**Schrödinger problems.** `SchrodingerProblem` holds $V$ by its coefficients and an energy;
`with_energy` re-tunes it.

**Turning points.** `turning_points` finds the roots of $Q$ with polished multiplicities;
`is_simple`, `order` and `location` inspect them. Higher-order turning points are detected
and rejected where the theory needs simple ones.

**WKB expansion.** `wkb_expansion` runs the Riccati recursion $S^2 + S' = \hbar^{-2}Q$ to any
order, keeping terms in the exact form $p(z)\,u^{\varepsilon}/Q^k$ rather than in a fraction
field. `even_odd_residual` is the built-in oracle that the recursion actually solves the
equation.

**Periods and Voros symbols.** `period_integral` integrates along piecewise-linear contours
with continuous branch tracking of $\sqrt{Q}$; `encircling_contour` builds a cycle around a
pair of turning points. `voros_symbol` assembles the resulting series, with the classical
period stored separately from the quantum corrections.

**Stokes graphs.** `stokes_graph` traces the lines $\mathrm{Im}(e^{-i\theta}\int\sqrt{Q}) = 0$
using an augmented state $(z, \sqrt{Q})$, which avoids discrete branch choices altogether.
`finite_lines`, `edges` and `topology_signature` describe the result combinatorially, and
`stokes_graph_family` sweeps $\theta$. `saddles` finds the saddle connections — the special
values of $\theta$ where the topology jumps — with their central charges $Z = 2\int\sqrt{Q}$.

**Spectra and exact quantization.** `spectral_cycles` picks out the well and barrier cycles
at a given energy and `quantum_period` re-runs the whole pipeline there.
`quantization_condition` implements all-orders Bohr–Sommerfeld for a single well and the
parity-factorized condition for a double well; `wkb_eigenvalue` solves it by Newton
iteration, and `spectral_determinant` gives the same information as a function whose zeros
are the spectrum.

**Double wells.** `perturbative_b` and `instanton_a` compute the two $\hbar$-series of the
Zinn-Justin–Jentschura framework, `energy_splitting` gives the even/odd gap, and `verify_zjj`
checks the P/NP relation between the two series order by order.

**Delabaere–Dillinger–Pham layer.** `voros_value` Borel sums a Voros symbol on either side of
a Stokes ray; `ddp_transform` applies the jump; `verify_ddp` confirms it numerically.
`ddp_seed` turns the intersection pairing into a `ClusterAlgebras.Seed`, and
`verify_ddp_mutation` checks that the jump agrees with cluster $y$-mutation.

**The cluster bridge.** `ideal_triangulation` converts a generic Stokes graph into a
triangulation of the $(d+2)$-gon by walking the faces of the planar map;
`triangulation_quiver` reads off the quiver; `flip` performs the (handed) combinatorial flip.
`charge_basis` builds the charge lattice with its central charges and intersection pairing,
and `bps_spectrum` produces the BPS states in phase order — which is a maximal green sequence
of the quiver, verified by replaying the c-vectors and cross-checked against the confirmed
saddles.

**Signed frames.** The physical charge lattice differs from the decay-oriented one by a sign
per basis cycle. `signed_frame` solves for those signs by inverting the keystone identity
$P = -\varepsilon_i\varepsilon_j B$ over the quiver graph, `verify_signed_frame` checks the
result, and `chamber_walls` / `reference_theta` locate the chamber the frame belongs to. The
signs are chamber-local: solvability in a chamber whose triangulation has an internal
triangle is a theorem, not bookkeeping.

**TBA.** `tba_system` and `solve_tba` implement the conformal-limit Gaiotto–Moore–Neitzke
integral equations, which determine the Borel-summed Voros symbols from the BPS data
$(Z_\gamma, \Omega(\gamma), \langle\gamma,\gamma'\rangle)$ alone. This is a genuinely
independent route to the same quantum periods — no WKB series, no Padé — and therefore the
sharpest available check on the series pipeline. `n_iterations` and `residual` report the
fixed-point solve.

**Seiberg–Witten $SU(2)$.** `SeibergWittenSU2` carries the pure $\mathcal{N}=2$ geometry by
its dynamical scale. `sw_periods` gives the electric and magnetic periods anywhere on the
$u$-plane in closed form, `quantum_sw_periods` adds the Nekrasov–Shatashvili $\hbar$
corrections, `continue_periods` transports them along an arbitrary path by the Picard–Fuchs
equation (loops around the singularities reproduce the exact integer `sw_monodromy`
matrices), and `ms_wall` / `sw_chamber` trace and classify by the wall of marginal stability.
`su2_bps_states` gives the BPS content in either chamber — two states strongly coupled, the
two dyon towers plus the $\Omega = -2$ W-boson weakly — and `verify_su2_wall_crossing` checks
the Kontsevich–Soibelman identity that relates them. This layer stands on its own: it uses
period quadratures directly, not the polynomial WKB engine.

**Numerical hygiene.** Structs are immutable and every operation returns a new object;
invalid input throws a typed error (`NonGenericGraph`, `CoalescentTurningPoints`, `TBAError`,
…) rather than a plausible-looking wrong number; precision comes from the argument type
rather than a global setting.

## Extensions

Plotting is loaded on demand.

```julia
using CairoMakie                     # or GLMakie
plot_stokes_graph(g)                 # one graph, or a vector for a θ-family
plot_triangulation(g, t)             # the dual triangulation on top of the graph
```

## Related packages

This is the bridge of a family of Julia packages for **exact and asymptotic methods** — where
a discrete, exactly computable structure controls a continuous, only-asymptotically-defined
one. [Resurgence.jl](https://github.com/benedikt-nagler/Resurgence.jl) is the continuous
foundation, providing the Borel–Padé summation, transseries and alien calculus used here;
[ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl) is the discrete
one, providing quivers, seed mutation, green sequences and DT invariants. Neither depends on
the other, or on this package.

## License

MIT
