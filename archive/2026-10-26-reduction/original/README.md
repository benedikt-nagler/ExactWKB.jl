# ExactWKB.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://benedikt-nagler.github.io/ExactWKB.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://benedikt-nagler.github.io/ExactWKB.jl/dev)
[![CI](https://github.com/benedikt-nagler/ExactWKB.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benedikt-nagler/ExactWKB.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Exact WKB analysis of one-dimensional Schrödinger equations in Julia: turning points, the
all-orders WKB expansion, quantum periods, Stokes graphs, exact quantization conditions, BPS
spectra, hyperkähler metrics, and the dictionary to cluster algebras.

The object is

$$\hbar^2 \psi''(z) = Q(z)\,\psi(z), \qquad Q(z) = V(z) - E,$$

with $z$ complex. The WKB ansatz turns this into a Riccati equation whose power series in
$\hbar$ diverges for every $\hbar \neq 0$. Borel summed, the series gives genuine solutions.
The Stokes graph, the level set $\mathrm{Im}(e^{-i\theta}\int\sqrt{Q}) = 0$, records which
solution you are on, and jumps as $\theta$ rotates.

Cluster algebras enter at that jump. By the Iwaki–Nakanishi dictionary a generic Stokes graph
is dual to an ideal triangulation, a Voros-symbol jump is a cluster $y$-mutation, and a full
circuit of the $\theta$-plane is a maximal green sequence. The package computes both sides and
checks that they agree.

$Q$ may be a polynomial or a rational function. Poles of a rational $Q$ put the triangulation
on a surface: an irregular pole (order $\geq 3$) is a boundary circle, a double pole is a
puncture. This reaches cluster types no polynomial does, such as Kronecker from the annulus
and $D_m$ from the once-punctured $m$-gon.

Divergent-series machinery comes from
[Resurgence.jl](https://github.com/benedikt-nagler/Resurgence.jl), cluster combinatorics from
[ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl).

## Installation

Not yet registered, but both foundations are, so Pkg resolves them by name.

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

Exact quantization gives energy levels:

```julia
harm = SchrodingerProblem([0.0, 0.0, 1.0])   # V = z², so Eₙ = 2ħ(n + ½)
wkb_eigenvalue(harm, 0, 0.1; order = 6)      # 0.1
wkb_eigenvalue(harm, 3, 0.1; order = 6)      # 0.7
```

The cluster bridge, on the cubic:

```julia
cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])  # Q = z³ − z
g     = stokes_graph(cubic; theta = 0.3)
t     = ideal_triangulation(g)       # dual triangulation of the (d+2)-gon
q     = triangulation_quiver(t)      # a ClusterAlgebras.Quiver, here type A₂
cb    = charge_basis(cubic, t; margin = 0.4)
signed_pairing(cb) == -q.B           # the keystone identity

sp = bps_spectrum(cubic; theta = 0.3, margin = 0.4)
charges(sp)                          # [[1, 0], [0, 1]]
sp.sequence                          # [1, 2] - the same data as a green sequence
solve_tba(sp)                        # the periods again, from the spectrum alone
```

A rational $Q$ reaches types a polynomial one cannot:

```julia
using ClusterAlgebras

p = RationalProblem([1.0, 0.0, 0.0, 0.0, 1.0], [0.0], [2])  # Q = (1 + z⁴)/z²
t = ideal_triangulation(stokes_graph(p; theta = 0.3))
n_punctures(t)                             # 1
mutation_type(triangulation_quiver(t))     # D4

m = mathieu_problem(1.0, 3.0)              # Λ = 1, u = 3, on the w = e^{2ix} plane
triangulation_quiver(ideal_triangulation(stokes_graph(m; theta = 0.3))).B   # [0 -2; 2 0]
```

## Features

- **Problems.** `SchrodingerProblem` holds $V$ by its coefficients and an energy;
  `with_energy` re-tunes it. `RationalProblem(numerator, poles, orders)` carries its poles as
  exact data, since a pole order is topological. `poles`, `pole_orders`, `puncture_indices`,
  `n_punctures`, `ring_domain_walls`, `mathieu_problem`.
- **Turning points.** `turning_points` finds the roots of $Q$ with polished multiplicities. A
  turning point of order $m$ emits $m + 2$ rays.
- **WKB expansion.** `wkb_expansion` runs $S^2 + S' = \hbar^{-2}Q$ to any order, keeping terms
  as $p(z)\,u^{\varepsilon}/Q^k$. `even_odd_residual` is the built-in oracle.
- **Periods and Voros symbols.** `period_integral` integrates along piecewise-linear contours
  with continuous branch tracking. `encircling_contour`, `voros_symbol`.
- **Stokes graphs.** `stokes_graph` traces lines in the augmented state $(z, \sqrt{Q})$, so no
  discrete branch choice is ever made. `finite_lines`, `edges`, `topology_signature`,
  `stokes_graph_family`, `saddles`, `is_degenerate`.
- **Degenerate graphs.** `polygon_decomposition` is the dual of a non-generic graph, with an
  $(m+2)$-gon per turning point of order $m$. `refinements` enumerates the triangulations
  refining it, and `n_refinements` counts them.
- **The Weber model.** `weber_model` is the parabolic-cylinder model of a merging pair of
  turning points. `weber_index`, `weber_voros_series`, `weber_voros_sum`, `weber_connection`,
  `weber_barrier_amplitude`.
- **Spectra.** `quantization_condition` implements all-orders Bohr–Sommerfeld and the
  parity-factorized double-well condition. `wkb_eigenvalue` solves it by Newton;
  `spectral_determinant`, `spectral_cycles`, `quantum_period`. Passing `uniform = true`
  reaches levels above the barrier top.
- **Parameter derivatives.** The mathematics is differentiated, not the program.
  `wkb_derivative` differentiates the Riccati recursion, `period_derivative` differentiates
  under the integral sign on a frozen closed cycle, and `quantization_derivative` returns
  residual and derivatives from one evaluation. `eigenvalue_sensitivity` gives
  $\partial E/\partial v_k$.
- **Double wells.** `perturbative_b`, `instanton_a`, `energy_splitting`, and `verify_zjj` for
  the Zinn-Justin–Jentschura relation.
- **Wall-crossing.** `voros_value` Borel sums on either side of a Stokes ray, `ddp_transform`
  applies the Delabaere–Dillinger–Pham jump, and `ddp_seed` / `verify_ddp_mutation` check it
  against cluster $y$-mutation.
- **The cluster bridge.** `ideal_triangulation`, `triangulation_quiver`, `flip`,
  `charge_basis`, `bps_spectrum`. Charges on a punctured surface come from each diagonal's own
  strip region.
- **Signed frames.** `signed_frame` solves for the signs relating the physical charge lattice
  to the decay-oriented one, with `verify_signed_frame`, `chamber_walls` and
  `reference_theta`. The signs are chamber-local.
- **TBA.** `tba_system` and `solve_tba` implement the conformal-limit Gaiotto–Moore–Neitzke
  equations. They fix the Voros symbols from BPS data alone, with no series and no Padé.
- **Hyperkähler metrics.** `gmn_torus`, `solve_gmn`, `metric_point`, `symplectic_expansion`,
  `hk_metric`, `semiflat_metric`, and the closed-form Ooguri–Vafa layer
  (`ooguri_vafa_xi`, `ooguri_vafa_instantons`).
- **Seiberg–Witten $SU(2)$.** `sw_periods` gives both periods in closed form on the whole
  $u$-plane, `quantum_sw_periods` adds Nekrasov–Shatashvili corrections to any order, and
  `continue_periods` transports them by Picard–Fuchs. `ms_wall`, `sw_chamber`,
  `su2_bps_states`, `verify_su2_wall_crossing`.
- **Typed errors.** `NonGenericGraph`, `CoalescentTurningPoints`, `TBAError` and friends,
  instead of a plausible wrong number. Precision comes from the argument type.

## Plotting

```julia
using CairoMakie                     # or GLMakie
plot_stokes_graph(g)                 # one graph, or a vector for a θ-family
plot_triangulation(g, t)             # the dual triangulation on top of the graph
```

## Documentation

The [documentation](https://benedikt-nagler.github.io/ExactWKB.jl/dev) has a manual page per
layer and three worked tutorials:

- **The cubic, end to end.** Turning points to Borel plane, Stokes graph, triangulation, $A_2$
  quiver, BPS spectrum, wall-crossing and TBA.
- **The double well.** Exact eigenvalues, the $10^{-6}$ level splitting, and the Dunne–Ünsal
  relation.
- **Seiberg–Witten $SU(2)$.** Periods on the $u$-plane, monodromy by transport, the wall of
  marginal stability, both chambers, and Kontsevich–Soibelman wall-crossing.

Build locally with `julia --project=docs docs/make.jl`.

[`examples/`](examples/) holds one runnable script per layer, including
[`ddp_demo.jl`](examples/ddp_demo.jl),
[`cluster_bridge_demo.jl`](examples/cluster_bridge_demo.jl),
[`spectra_demo.jl`](examples/spectra_demo.jl),
[`double_well_demo.jl`](examples/double_well_demo.jl) and
[`tba_demo.jl`](examples/tba_demo.jl).
[`su2_four_package_loop.jl`](examples/su2_four_package_loop.jl) computes pure $SU(2)$ twice by
disjoint routes and is the only program in which all four packages meet.

## Related packages

This is the bridge of a family of Julia packages for exact and asymptotic methods.
[Resurgence.jl](https://github.com/benedikt-nagler/Resurgence.jl) supplies Borel–Padé
summation, transseries and alien calculus.
[ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl) supplies quivers,
seed mutation, green sequences and DT invariants. Neither depends on the other, or on this
package. [ClusterSurfaces.jl](https://github.com/benedikt-nagler/ClusterSurfaces.jl) turns the
quiver of a Stokes graph back into a marked surface.

Names from those packages are not re-exported. Use `import Resurgence` /
`import ClusterAlgebras` and qualify.

## License

MIT
