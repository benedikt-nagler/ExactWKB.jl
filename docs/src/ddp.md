```@meta
CurrentModule = ExactWKB
```

# Wall-crossing

Borel summing a Voros symbol requires a direction in the Borel plane. When that direction hits
a singularity, that is a Stokes ray, the two lateral sums differ, and the difference is fixed:
the Delabaere–Dillinger–Pham formula says the jump across the ray of a saddle ``\gamma'`` is

```math
s_+(V_\gamma) = s_-(V_\gamma)\,\bigl(1 + V_{\gamma'}\bigr)^{\langle\gamma,\gamma'\rangle},
```

with ``\langle\gamma,\gamma'\rangle`` the intersection pairing of the two cycles, an
*integer*. The Stokes constant is therefore exact combinatorial data, while the summed values
on either side are numerical. [`verify_ddp`](@ref) measures the exponent from the two lateral
sums and returns it next to its residual, so the integer is a result and not an input.

Written in terms of ``y_j = V_{\gamma_j}`` the same formula *is* cluster ``y``-mutation. With
the exchange matrix ``B = -P`` built from the pairing ([`ddp_seed`](@ref)), crossing the wall
at ``k`` acts on the ``y``-variables exactly as Fomin–Zelevinsky's ``\mu_k``, which
[`verify_ddp_mutation`](@ref) checks against `ClusterAlgebras`. Mutating in the wrong direction
fails the check by an order of magnitude, so this tests the orientation and is not a tautology.

The Voros symbols themselves never move. Crossing ``\theta_c`` only rotates the Laplace ray in
the Borel plane, so the jump is a property of the summation and not of the series.

Every sign in this layer is pinned by an oracle; see the convention ledger at the head of
`src/ddp.jl`. In particular the wall symbol must *decay* along the wall
(``\mathrm{Re}(Z e^{-i\theta_c}) < 0``), and a reversed one raises a typed error instead of
being silently flipped.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/ddp.jl"]
```
