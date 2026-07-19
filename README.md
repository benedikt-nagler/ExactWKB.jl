# ExactWKB.jl

[![CI](https://github.com/benedikt-nagler/ExactWKB.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/benedikt-nagler/ExactWKB.jl/actions/workflows/CI.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Exact WKB analysis for one-dimensional Schrödinger equations with polynomial potentials:
turning points, the all-orders WKB expansion, quantum periods and Voros symbols, Stokes
graphs, exact quantization conditions, and the dictionary to cluster algebras.

The starting point is

$$\hbar^2 \psi''(z) = Q(z)\,\psi(z), \qquad Q(z) = V(z) - E,$$

with $Q$ a polynomial. The usual WKB ansatz $\psi = \exp\left(\hbar^{-1}\int S\right)$ turns
this into a Riccati equation whose solution is a power series in $\hbar$. That series
diverges, and exact WKB is the statement that it can nevertheless be Borel summed. Doing so
carefully gives real answers: Voros symbols (Borel sums of period integrals) that satisfy
exact quantization conditions, and Stokes graphs whose combinatorics changes discontinuously
as one rotates the summation direction $\theta$.

The last part is where cluster algebras enter. Iwaki and Nakanishi observed that a generic
Stokes graph is dual to an ideal triangulation of a polygon, that the jump of the Voros
symbols across a wall is exactly a cluster $y$-mutation, and that a full circuit of the
$\theta$-plane is a maximal green sequence. The package computes both sides and checks that
they agree.

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

## Functionality

**Schrödinger problems.** `SchrodingerProblem` holds $V$ by its coefficients and an energy;
`with_energy` re-tunes it. Exact `Rational` input flows through to `BigFloat` at whatever
precision the caller has set, so the same problem can be solved cheaply in `Float64` or
carefully in high precision.

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

**Stokes graphs.** `stokes_graph` traces the lines $\operatorname{Im}(e^{-i\theta}\int\sqrt{Q}) = 0$
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
`signed_frame` fixes the orientation signs in each chamber, and `bps_spectrum` produces the
BPS states in phase order — which is a maximal green sequence of the quiver, verified by
replaying the c-vectors and cross-checked against the confirmed saddles.

## Extensions

Plotting is loaded on demand.

```julia
using CairoMakie                     # or GLMakie
plot_stokes_graph(g)                 # one graph, or a vector for a θ-family
plot_triangulation(g, t)             # the dual triangulation on top of the graph
```

## Examples

- [`examples/stokes_graphs_intro.ipynb`](examples/stokes_graphs_intro.ipynb) — the basics,
  from turning points to Stokes graphs
- [`examples/m3_flagship.jl`](examples/m3_flagship.jl) — the full WKB chain on the double well
- [`examples/m4_flagship.jl`](examples/m4_flagship.jl) — the cluster bridge end to end
- [`examples/m5_flagship.jl`](examples/m5_flagship.jl) — exact quantization and level splitting

## Related packages

[Resurgence.jl](https://github.com/benedikt-nagler/Resurgence.jl) provides the Borel–Padé
summation, transseries and alien calculus used here.
[ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl) provides quivers,
seed mutation, green sequences and DT invariants.

## References

- A. Voros, *The return of the quartic oscillator: the complex WKB method*, Ann. Inst. H. Poincaré A **39** (1983), 211–338.
- E. Delabaere, H. Dillinger, F. Pham, *Résurgence de Voros et périodes des courbes hyperelliptiques*, Ann. Inst. Fourier **43** (1993), 163–199.
- J. Zinn-Justin, U. D. Jentschura, *Multi-instantons and exact results I, II*, Ann. Phys. **313** (2004), 197–267, 269–325.
- T. Kawai, Y. Takei, *Algebraic Analysis of Singular Perturbation Theory*, Transl. Math. Monogr. **227**, AMS (2005).
- D. Gaiotto, G. W. Moore, A. Neitzke, *Wall-crossing, Hitchin systems, and the WKB approximation*, Adv. Math. **234** (2013), 239–403.
- G. V. Dunne, M. Ünsal, *Uniform WKB, multi-instantons, and resurgent trans-series*, Phys. Rev. D **89** (2014), 105009.
- K. Iwaki, T. Nakanishi, *Exact WKB analysis and cluster algebras*, J. Phys. A **47** (2014), 474009.
- T. Bridgeland, I. Smith, *Quadratic differentials as stability conditions*, Publ. Math. IHÉS **121** (2015), 155–278.

## License

MIT
