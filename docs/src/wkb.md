```@meta
CurrentModule = ExactWKB
```

# The WKB expansion

The Wentzel-Kramers-Brillouin ansatz ``\psi = \exp\left(\hbar^{-1}\int S\,dz\right)`` turns the Schrödinger equation
into the Riccati equation

```math
S^2 + \hbar\, S' = \hbar^{-2} Q,
```

solved order by order from ``S_{-1} = \sqrt{Q}`` and ``S_0 = -Q'/(4Q)``. Every term is exactly
of the form ``p_m(z)\,u^{\varepsilon_m}/Q^{k_m}`` with ``u = \sqrt{Q}``, so
[`wkb_expansion`](@ref) stores the triple ``(p_m, k_m, \varepsilon_m)`` and never builds a
rational-function field - which is what keeps the `BigFloat` coefficient ring gcd-free and
numerically safe.

Only the odd part is needed: ``S_{\text{even}} = -\tfrac12 \partial \log S_{\text{odd}}`` is a
total derivative and integrates to zero around a cycle. That reduction is never *used* to
compute anything - it is enforced as an oracle by [`even_odd_residual`](@ref), which should
vanish to working precision on any expansion you build.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/wkb_recursion.jl"]
```

## Periods

[`period_integral`](@ref) integrates along a piecewise-linear contour with continuous branch
tracking of ``\sqrt{Q}`` - each vertex carries a reference value, so no discrete sign choice is
ever made. [`encircling_contour`](@ref) builds a cycle around a pair of turning points.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/periods.jl"]
```

## Voros symbols

A [`VorosSymbol`](@ref) is the quantum period ``\oint S\,dz`` on a cycle, with the classical
part ``\hbar^{-1}\oint\sqrt{Q}`` kept separate from the quantum corrections. The classical
period is the central charge of the corresponding BPS state; the quantum series is a divergent
`Resurgence.FormalSeries` in ``\hbar``, ready for Borel–Padé summation
([Wall-crossing](ddp.md)) - and its Borel singularities land on the central-charge lattice of
the [BPS spectrum](bridge.md).

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/voros.jl"]
```
