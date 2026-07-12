# ExactWKB.jl

Exact WKB analysis of one-dimensional Schrödinger problems: turning points, the
all-orders WKB recursion, quantum periods and Voros symbols, Stokes graphs and their
BPS saddles.

## Introduction

For a one-dimensional Schrödinger equation in the Iwaki–Nakanishi normalization,

$$\hbar^2 \psi''(z) = Q(z)\,\psi(z), \qquad Q = V - E,$$

the WKB ansatz $\psi = \exp\!\big(\hbar^{-1}\!\int S\,dz\big)$ turns the equation into a
Riccati recursion whose all-orders solution $S = \sum_{m\ge -1} \hbar^{m} S_m$ is exact
but divergent. Its **quantum periods** $\oint S_m\,dz$ — assembled into the **Voros
symbol** — are the observables of the theory: Bohr–Sommerfeld quantization, the exact
spectrum, and (through Borel summation, provided by
[`Resurgence.jl`](../Resurgence.jl)) the non-perturbative completions.

The complex $z$-plane is organized by the **Stokes graph**: the trajectories
$\mathrm{Im}\big[e^{-i\theta}\!\int^z \sqrt{Q}\,dz\big] = 0$ emanating from the turning
points. A trajectory that connects two turning points is a **saddle** — a BPS state
with central charge $Z = 2\int\sqrt{Q}\,dz$ and mass $|Z|$. As $\theta$ sweeps, the
topology jumps at the critical phases $\arg Z \bmod \pi$: wall-crossing. This
Iwaki–Nakanishi dictionary between Stokes graphs and cluster algebras is the payoff the
next milestone (M4) builds on — which is why this package depends on both
`Resurgence.jl` and `ClusterAlgebras.jl`.

## Features (M3)

- **Potentials & turning points**: `SchrodingerProblem` for polynomial $Q = V - E$;
  `turning_points` via polished, clustered polynomial roots (exact-rational potentials
  are computed in `BigFloat` at the caller's `setprecision`).
- **All-orders WKB**: `wkb_expansion` runs the Riccati recursion in an exact ring
  (rational or `BigFloat`); the Riccati residual and hand-derived Airy $S_m$ are test
  oracles.
- **Quantum periods & Voros symbols**: `period_integral` with continuous $\sqrt{Q}$
  branch tracking; `voros_symbol` packages the series in exactly the shape
  `Resurgence.jl`'s `borel`/`pade`/`laplace_sum` consume.
- **Stokes graphs**: `stokes_graph` traces trajectories with an augmented state
  $(z,\sqrt{Q})$ — the mass parameter is the integration time and branch choices never
  arise — with `finite_lines`, `edges`, `n_infinite_lines`, `topology_signature`.
- **Saddles / BPS states**: `saddles` returns the confirmed saddle connections with
  their central charges, masses, and critical phases; `stokes_graph_family` sweeps
  $\theta$.
- **Plotting**: `plot_stokes_graph` (Makie extension) renders a graph or a
  $\theta$-family.

`examples/m3_flagship.jl` runs the full chain on a double well.

## Installation

`Resurgence.jl` and `ClusterAlgebras.jl` are unregistered — develop all three from a
local checkout:

```julia
using Pkg
Pkg.develop([PackageSpec(path = "../Resurgence.jl"),
             PackageSpec(path = "../ClusterAlgebras.jl"),
             PackageSpec(path = ".")])
```

## Quick start

```julia
using ExactWKB

# symmetric double well  Q = (z²−1)² − 1/4
prob = SchrodingerProblem([3//4, 0, -2, 0, 1])
tps  = turning_points(prob)                       # four simple turning points

w  = wkb_expansion(prob; order = 6)               # all-orders WKB
vs = voros_symbol(w, encircling_contour(tps[2], tps[3]))   # inner (instanton) cycle

g = stokes_graph(prob; theta = 0.0)               # Stokes graph at θ = 0
topology_signature(g)                             # (4, 10, [(2, 3)])
saddles(prob)                                     # 3 BPS states: two well cycles + instanton
```

## Status

M3 (Schrödinger / Stokes-graph core) is complete. M4 — the cluster bridge — is the next
milestone; see `PLANNING/roadmap.md` and the ecosystem `../PLANNING/roadmap.md`.

## License

MIT.
