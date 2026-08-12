```@meta
CurrentModule = ExactWKB
```

# The cluster bridge

From a polynomial potential alone, with no hand-supplied combinatorics, the package produces
the cluster algebra that governs its wall-crossing.
[The cubic, end to end](@ref) walks through the whole chain on one example. This page is the
reference for each step.

The chain is

```
stokes_graph → ideal_triangulation → triangulation_quiver
                     ↓
              charge_basis  (cycles, central charges, pairing, signs)
                     ↓
              bps_spectrum  (a maximal green sequence of that quiver)
```

## From graph to triangulation

A generic Stokes graph of a degree-``d`` potential is dual to an ideal triangulation of a
``(d+2)``-gon: each turning point gives a triangle, each Stokes region a vertex.
[`ideal_triangulation`](@ref) performs the planar face walk and asserts every genericity
invariant on the way, throwing [`NonGenericGraph`](@ref) instead of guessing when the graph has
a saddle connection or a degenerate face. [`triangulation_quiver`](@ref) then reads the quiver
off the oriented internal triangles. For polynomial potentials with simple turning points it is
always finite type ``A_{d-1}``, and it equals ``-\,``[`signed_pairing`](@ref) exactly.

[`flip`](@ref) is the combinatorial flip of a diagonal. It is *handed*: `direction = ±1`
distinguishes the two ways a diagonal can rotate, which matters because the sign frame below is
not transported across a wall.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/triangulation.jl"]
```

## The charge lattice

[`charge_basis`](@ref) attaches numbers to the combinatorics: one cycle per diagonal, its
central charge ``Z_\gamma = \oint\sqrt{Q}``, and the intersection pairing ``P``. The keystone
identity

```math
P_{ij} = -\varepsilon_i \varepsilon_j B_{ij}
```

relates the numerically computed pairing to the combinatorially computed quiver, and
`charge_basis` enforces it strictly. A failure raises an error, since the alternative is a
plausible-looking BPS spectrum.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/charge_lattice.jl"]
```

## The signed frame

The decay rule that orients a cycle fixes a representative of ``\pm\gamma``, not the physical
charge. The physical cycle is ``\gamma^{\text{phys}}_e = \varepsilon_e \gamma_e``, and the
signs ``\varepsilon`` are solved by inverting the keystone: a spanning-tree propagation over
the quiver graph, plus a cocycle condition that must close around every cycle. That condition
is vacuous on tree chambers and genuine on cyclic ones, which is the case an unsigned frame
cannot handle. With it, every chamber is supported.

``\varepsilon`` is chamber-*local*: it is not transported, and it is not continuous across a
wall even where the cycle itself persists. The residual gauge freedom, one global sign per
connected component, is physically inert.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/signed_frame.jl"]
```

## The spectrum of stable states

The saddle connections of the Stokes graph are the Bogomolny-Prasad-Sommerfield (BPS) states of
the problem. [`bps_spectrum`](@ref) sweeps the walls in the order they are crossed as
``\theta`` decreases from the chamber's own ``\theta_0``, collecting the charges as
``c``-vectors. The result is a maximal green sequence of the triangulation quiver, and it is
verified three ways: by replaying the ``c``-vectors through `ClusterAlgebras`, by attaching the
Donaldson-Thomas invariants ``\Omega``, and by cross-checking against the independently traced
[saddle connections](stokes.md). The spectrum does not depend on which chamber you start in.
Only its ordering does.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/bps.jl"]
```

## Degenerate graphs

A Stokes graph whose turning points have collided has no triangulation, but it still has a
dual. [`polygon_decomposition`](@ref) builds it: where a triangulation has triangles, a
[`PolygonDecomposition`](@ref) has an ``(m+2)``-gon per turning point of order ``m``, and
[`is_triangulation`](@ref) tells the two apart ([`n_cells`](@ref), [`cell_sizes`](@ref)
describe it). The same face walk produces both, since it is written in terms of the valences
``m+2`` instead of the constant 3, so a degenerate graph is not a special case in the code.

[`refinements`](@ref) enumerates every ideal triangulation refining a decomposition
([`n_refinements`](@ref) counts them). This is where the statement *a double turning point is
a flip wall* becomes checkable: a square cell has exactly two refinements, they differ by a
single [`flip`](@ref), and they are the triangulations traced on either side of the critical
energy, so the quivers on the two sides differ by one mutation.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/polygon_decomposition.jl"]
```

## The refined (motivic) bridge

Read in the quantum torus instead of its classical limit, a BPS chamber becomes a product of
quantum dilogarithms and [`refined_dt`](@ref) is its refined Donaldson–Thomas invariant,
evaluated through `ClusterAlgebras.jl`'s truncated Kontsevich–Soibelman algebra with the skew
form supplied by the bridge's own `signed_pairing`. [`verify_refined_wall_crossing`](@ref)
checks that two chambers give the same product, and [`refined_ddp_jump`](@ref) is the refined
form of the Delabaere–Dillinger–Pham jump, namely conjugation by a quantum dilogarithm, whose
``v \to 1`` limit is the classical factor.

That asymmetry is the content of the refinement: the *transformation* has a classical limit,
while the *product* does not. On the WKB side a Voros symbol is a Borel-summed number and its
lateral jump is a number, so no ``v`` is measurable. The refinement is algebraic, and a
genuinely ``q``-deformed WKB side would need non-commuting Voros symbols, that is, quantum
Fock–Goncharov coordinates.
[`refined_quiver_iso`](@ref) supplies the identification used when comparing chambers reached
from *different* potentials.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/refined.jl"]
```
