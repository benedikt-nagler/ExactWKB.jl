```@meta
CurrentModule = ExactWKB
```

# The cluster bridge

This is the payoff layer: from a polynomial potential alone, with no hand-supplied
combinatorics, the package produces the cluster algebra that governs its wall-crossing. [The cubic, end to end](@ref) walks through the
whole chain on one example; this page is the reference for each step.

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
invariant on the way, throwing [`NonGenericGraph`](@ref) if the graph has a saddle connection
or a degenerate face rather than guessing. [`triangulation_quiver`](@ref) then reads the quiver
off the oriented internal triangles - for polynomial potentials with simple turning points this
is always finite type ``A_{d-1}``, and it equals ``-\,``[`signed_pairing`](@ref) exactly.

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
`charge_basis` enforces it strictly - if it fails, you get an error rather than a plausible
BPS spectrum.

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
is vacuous on tree chambers and genuine on cyclic ones - which is exactly the case an
unsigned frame cannot handle. With it, every chamber is supported.

``\varepsilon`` is chamber-*local*: it is not transported, and it is not continuous across a
wall even where the cycle itself persists. The residual gauge freedom (one global sign per
connected component) is physically inert.

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
[saddle connections](stokes.md). The spectrum
does not depend on which chamber you start in - only its ordering does.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/bps.jl"]
```
