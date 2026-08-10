# ExactWKB.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://benedikt-nagler.github.io/ExactWKB.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://benedikt-nagler.github.io/ExactWKB.jl/dev)
[![CI](https://github.com/benedikt-nagler/ExactWKB.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benedikt-nagler/ExactWKB.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Exact WKB analysis of one-dimensional Schrödinger equations with polynomial or rational
potentials - equivalently, of second-order linear ODEs with a small parameter: turning points,
the all-orders WKB expansion, quantum periods and Voros symbols, Stokes graphs, exact
quantization conditions, BPS spectra and wall-crossing, hyperkähler metrics, and the
dictionary to cluster algebras.

The object is the Schrödinger equation in its complex-analytic form

$$\hbar^2 \psi''(z) = Q(z)\,\psi(z), \qquad Q(z) = V(z) - E,$$

with $z$ complex - equivalently the general second-order linear ODE with a large parameter, in
normal form. The WKB ansatz turns it into a Riccati equation solved by a power series in
$\hbar$ that diverges for every $\hbar \neq 0$. **Exact WKB** makes it usable anyway: Borel
summed, it gives genuine solutions rather than an approximation, and the *Stokes graph* - the
level set $\mathrm{Im}(e^{-i\theta}\int\sqrt{Q}) = 0$ - records which solution you are on,
jumping discontinuously as the summation direction $\theta$ rotates.

That last jump is where cluster algebras enter. By the Iwaki–Nakanishi dictionary a generic
Stokes graph is dual to an ideal triangulation, the jump of the Voros symbols across a wall is
a cluster $y$-mutation, and a full circuit of the $\theta$-plane is a maximal green sequence.
The package computes both sides and checks they agree - the combinatorics is exact, so it
certifies the numerics.

**Scope.** $Q$ may be a polynomial or a rational function. The polynomial case is the harmonic,
quartic, double-well and cubic oscillators and everything else usually used to test the
theory. The rational case adds the poles, and with them the surfaces: an **irregular** pole
(order $\geq 3$) is a boundary circle of the surface the triangulation lives on, a **regular**
one (a double pole) is a puncture, so a rational $Q$ reaches cluster types a polynomial one
cannot - the annulus gives the Kronecker quiver and the once-punctured $m$-gon gives $D_m$.
The Mathieu equation is one of these, which is how the WKB side and the Seiberg–Witten curve
side meet on the same problem. Exact `Rational` input flows through to `BigFloat` at whatever
precision the caller has set.

`ExactWKB.jl` is the bridge package of a small ecosystem: divergent-series machinery comes
from [Resurgence.jl](https://github.com/benedikt-nagler/Resurgence.jl), cluster combinatorics
from [ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl).

## Installation

Not yet registered in General, but both of its foundations are, so Pkg resolves them by
name. From the Julia REPL:

```julia
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
asymptotic - which makes it a good convention check.

And the cluster bridge, on the cubic:

```julia
cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])  # Q = z³ − z
g     = stokes_graph(cubic; theta = 0.3)
t     = ideal_triangulation(g)       # dual triangulation of the (d+2)-gon
q     = triangulation_quiver(t)      # a ClusterAlgebras.Quiver - here type A₂
cb    = charge_basis(cubic, t; margin = 0.4)   # charges, central charges, pairing
signed_pairing(cb) == -q.B           # the keystone identity

sp = bps_spectrum(cubic; theta = 0.3, margin = 0.4)
charges(sp)                          # [[1, 0], [0, 1]] - BPS charges in phase order
sp.sequence                          # [1, 2] - the same data as a green sequence
```

The BPS spectrum is enough to reconstruct the periods a second time, with no series and no
Padé anywhere in the pipeline:

```julia
sol = solve_tba(sp)                  # conformal-limit GMN integral equations
```

A rational $Q$ puts the same machinery on a surface with holes and punctures. A double pole is
a puncture, and the once-punctured square is cluster type $D_4$ - a type no polynomial
potential reaches:

```julia
using ClusterAlgebras

p = RationalProblem([1.0, 0.0, 0.0, 0.0, 1.0], [0.0], [2])  # Q = (1 + z⁴)/z², pole order 2
t = ideal_triangulation(stokes_graph(p; theta = 0.3))
n_punctures(t)                       # 1
mutation_type(triangulation_quiver(t))     # D4
```

The Mathieu equation is the rational problem behind Seiberg–Witten $SU(2)$: on the
$w = e^{2ix}$ plane it has an order-3 pole at each end, the origin and infinity, so the
surface is an annulus and the quiver is Kronecker.

```julia
m = mathieu_problem(1.0, 3.0)        # Λ = 1, u = 3, on the w = e^{2ix} plane
triangulation_quiver(ideal_triangulation(stokes_graph(m; theta = 0.3))).B   # [0 -2; 2 0]
```

## Documentation

The [documentation](https://benedikt-nagler.github.io/ExactWKB.jl/dev) has a manual page per
layer and three worked tutorials that run the machinery end to end:

- **The cubic, end to end** - one potential from turning points to Borel plane, Stokes graph,
  triangulation, $A_2$ quiver, BPS spectrum, wall-crossing and TBA.
- **The double well** - spectra beyond all orders: exact eigenvalues, the $10^{-6}$ level
  splitting, and the Dunne–Ünsal relation between the perturbative and instanton series.
- **Seiberg–Witten $SU(2)$** - periods on the $u$-plane, monodromy measured by transport, the
  wall of marginal stability, both BPS chambers, and Kontsevich–Soibelman wall-crossing.

Build it locally with `julia --project=docs docs/make.jl`.

## Examples

[`examples/`](examples/) holds runnable scripts, one per layer:
[`ddp_demo.jl`](examples/ddp_demo.jl) (Borel summing a Voros symbol on both sides of a Stokes
ray, applying the Delabaere–Dillinger–Pham jump, and checking it against a cluster
$y$-mutation), [`cluster_bridge_demo.jl`](examples/cluster_bridge_demo.jl),
[`spectra_demo.jl`](examples/spectra_demo.jl),
[`double_well_demo.jl`](examples/double_well_demo.jl), [`tba_demo.jl`](examples/tba_demo.jl)
and the Seiberg–Witten trio.

[`examples/su2_four_package_loop.jl`](examples/su2_four_package_loop.jl) is the one that
closes the circle: pure $SU(2)$ computed twice by disjoint routes - Mathieu through the Stokes
graph, the annulus, the Kronecker quiver and back to a surface, against the same periods read
off the Seiberg–Witten curve by elliptic integrals. It is a script rather than a test because
it is the only program in which all four packages meet.

## Functionality

**Schrödinger problems.** `SchrodingerProblem` holds $V$ by its coefficients and an energy;
`with_energy` re-tunes it. Every later layer dispatches on the abstract
`AbstractSchrodingerProblem`, so a potential class only has to supply its value, its Taylor
jets and its turning points.

**Rational problems.** `RationalProblem(numerator, poles, orders)` is the second such class:
$Q = N(z)/\prod (z - p_j)^{m_j}$, with the poles carried as **exact** data rather than
root-found, because a pole order is a topological datum of the surface and not a numerical
one. `poles` / `pole_orders` / `n_finite_poles` describe them, `puncture_indices` and
`n_punctures` single out the double poles, and `ring_domain_walls` locates the directions
$\theta \equiv \arg(2\pi i\sqrt{c})$ where a puncture's Stokes spirals close into a ring
domain - a wall of its own, belonging to no turning-point pair. `mathieu_problem(Λ, u)` is the
Mathieu equation as one of these.

**Turning points.** `turning_points` finds the roots of $Q$ with polished multiplicities.
Higher-order turning points are traced rather than refused - a turning point of order $m$
emits $m + 2$ rays - and are rejected only by the layers whose theorems need simple ones.

**WKB expansion.** `wkb_expansion` runs the Riccati recursion $S^2 + S' = \hbar^{-2}Q$ to any
order, keeping terms in the exact form $p(z)\,u^{\varepsilon}/Q^k$ rather than in a fraction
field. `even_odd_residual` is the built-in oracle that it solves the equation.

**Periods and Voros symbols.** `period_integral` integrates along piecewise-linear contours
with continuous branch tracking of $\sqrt{Q}$, `encircling_contour` builds a cycle around a
pair of turning points, and `voros_symbol` assembles the series with the classical period kept
separate from the quantum corrections.

**Stokes graphs.** `stokes_graph` traces the lines using an augmented state $(z, \sqrt{Q})$,
avoiding discrete branch choices altogether; `finite_lines`, `edges` and `topology_signature`
describe the result combinatorially and `stokes_graph_family` sweeps $\theta$. `saddles` finds
the saddle connections - where the topology jumps - with their central charges. A pole of $Q$
is a second escape boundary alongside infinity, so a ray may end there instead, and
`is_degenerate` / `turning_point_orders` report a graph whose turning points have collided.

**Degenerate graphs.** `polygon_decomposition` is the dual of a Stokes graph that is *not*
generic: where a triangulation has triangles, it has an $(m+2)$-gon per turning point of order
$m$, and `is_triangulation` tells the two apart. `refinements` enumerates every ideal
triangulation refining it (`n_refinements` counts them, a product of Catalan numbers). This is
where the statement **a double turning point is a flip wall** lives: a square cell has exactly
two refinements, they differ by one `flip`, and they are the triangulations traced on either
side of the critical energy - so the quivers differ by one mutation.

**The Weber model.** `weber_model` is the parabolic-cylinder local model of a *merging* pair
of turning points, with `weber_index` $\nu = Z/(2\pi i\hbar)$, the exact Stirling-type
`weber_voros_series` and its Borel sum `weber_voros_sum` (valid on the imaginary axis too,
where the singularities sit on the Laplace ray), `weber_connection`, and
`weber_barrier_amplitude` - the connection constant a merging pair contributes where a simple
turning point contributes the Airy $\tfrac{1}{2}$.

**Spectra and exact quantization.** `quantization_condition` implements all-orders
Bohr–Sommerfeld for a single well and the parity-factorized condition for a double well;
`wkb_eigenvalue` solves it by Newton iteration and `spectral_determinant` gives the same
information as a function whose zeros are the spectrum. `spectral_cycles` and `quantum_period`
supply the cycles at a given energy. Passing `uniform = true` swaps the Airy connection
constant for the Weber one and runs on `uniform_cycles`, which continue both cycles *through*
the barrier top - so levels above the barrier, where the classical turning-point layout stops
existing, become computable rather than out of range.

**Double wells.** `perturbative_b` and `instanton_a` compute the two $\hbar$-series of the
Zinn-Justin–Jentschura framework, `energy_splitting` gives the even/odd gap, and `verify_zjj`
checks the P/NP relation order by order.

**Delabaere–Dillinger–Pham layer.** `voros_value` Borel sums a Voros symbol on either side of
a Stokes ray, `ddp_transform` applies the jump, `verify_ddp` confirms it numerically, and
`ddp_seed` / `verify_ddp_mutation` check the jump against cluster $y$-mutation.

**The cluster bridge.** `ideal_triangulation` converts a generic Stokes graph into a
triangulation of the surface it lives on, `triangulation_quiver` reads off the quiver, and
`flip` performs the (handed) combinatorial flip. For a polynomial $Q$ that surface is the
$(d+2)$-gon; in general it is a sphere carrying one boundary circle per irregular singularity
of $Q\,dz^2$ (with `pole order − 2` marked points each) and one interior vertex per regular
one, which `n_punctures` and `puncture_valence` report. `charge_basis` builds the charge
lattice with its central charges and pairing - on a punctured surface from each diagonal's own
strip region, since two diagonals can share a turning-point pair and differ by a loop around a
hole. `bps_spectrum` produces the BPS states in phase order - a maximal green sequence of the
quiver, verified by replaying the c-vectors and cross-checked against the confirmed saddles.

**Signed frames.** The physical charge lattice differs from the decay-oriented one by a sign
per basis cycle. `signed_frame` solves for those signs by inverting the keystone identity
$P = -\varepsilon_i\varepsilon_j B$, `verify_signed_frame` checks the result, and
`chamber_walls` / `reference_theta` locate the chamber - the signs are chamber-local.

**TBA.** `tba_system` and `solve_tba` implement the conformal-limit Gaiotto–Moore–Neitzke
integral equations, which fix the Borel-summed Voros symbols from the BPS data
$(Z_\gamma, \Omega(\gamma), \langle\gamma,\gamma'\rangle)$ alone - an independent route to the
same quantum periods, with no WKB series and no Padé, and therefore the sharpest available
check on the series pipeline.

**Hyperkähler metrics.** The finite-radius companion of the TBA layer, following
Gaiotto–Moore–Neitzke's twistor construction. `gmn_torus` packages lattice-level BPS data,
`solve_gmn` solves the integral equations for the twistor coordinates $\mathcal{X}_\gamma$,
and `metric_point` / `symplectic_expansion` / `hk_metric` assemble the metric from the
holomorphic symplectic form - the Laurent-degree fit of $\varpi(\zeta)$ *is* the twistor
theorem run as a numerical test. `semiflat_metric` is the uncorrected comparison, and the
Ooguri–Vafa layer (`ooguri_vafa_xi`, `ooguri_vafa_instantons`) solves the one-state case a
second, independent way in closed form. `su2_torus` feeds the Seiberg–Witten spectrum in
directly, and the resulting metric agrees across the wall of marginal stability - which is
Kontsevich–Soibelman wall-crossing at the level of a metric.

**Seiberg–Witten $SU(2)$.** `sw_periods` gives the electric and magnetic periods anywhere on
the $u$-plane in closed form, `quantum_sw_periods` adds the Nekrasov–Shatashvili $\hbar$
corrections to **any** order (one universal cycle-independent operator per order, with no
Voros regularization anywhere), and `continue_periods` transports them along any path by the
Picard–Fuchs equation (loops reproduce the exact integer `sw_monodromy` matrices). `ms_wall` /
`sw_chamber` handle the wall of marginal stability - including the NS-deformed wall, which
moves outward as $\hbar^2$ - `su2_bps_states` gives the BPS content in either chamber, and
`verify_su2_wall_crossing` checks the Kontsevich–Soibelman identity relating them. This layer
uses period quadratures directly, not the polynomial WKB engine; `mathieu_problem` is the same
physics coming the other way, through the Stokes graph.

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

This is the bridge of a family of Julia packages for **exact and asymptotic methods** - where
a discrete, exactly computable structure controls a continuous, only-asymptotically-defined
one. [Resurgence.jl](https://github.com/benedikt-nagler/Resurgence.jl) is the continuous
foundation, providing the Borel–Padé summation, transseries and alien calculus used here;
[ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl) is the discrete
one, providing quivers, seed mutation, green sequences and DT invariants. Neither depends on
the other, or on this package.
[ClusterSurfaces.jl](https://github.com/benedikt-nagler/ClusterSurfaces.jl) sits alongside as
the geometric front end of the same core, and turns the quiver of a Stokes graph back into the
marked surface the graph lives on - not a dependency here, but the return leg of the
four-package loop above.

## License

MIT
