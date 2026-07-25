```@meta
CurrentModule = ExactWKB
```

# Thermodynamic Bethe ansatz

The spectrum of stable (Bogomolny-Prasad-Sommerfield) states is a finite list of integers and
complex numbers: ``(Z_\gamma, \Omega(\gamma), \langle\gamma,\gamma'\rangle)``. Remarkably, that list is enough
to reconstruct the Borel-summed Voros symbols exactly - no Wentzel-Kramers-Brillouin series,
no Padé anywhere in the pipeline. The mechanism is the conformal limit of the Gaiotto–Moore–Neitzke construction: the
symbols ``\mathcal X_\gamma(\hbar)`` are the solution of a system of integral equations of
thermodynamic-Bethe-ansatz type,

```math
\log \mathcal X_\gamma(\theta) = \frac{-Z_\gamma}{\hbar}
  + \sum_{\gamma'} \frac{\Omega(\gamma')\,\langle\gamma,\gamma'\rangle}{4\pi i}
    \int_{\ell_{\gamma'}} \frac{d\theta'}{\theta - \theta'}\,
    \log\bigl(1 - \mathcal X_{\gamma'}(\theta')\bigr),
```

whose only input is the spectrum. [`tba_system`](@ref) extracts that state-level data from a
[`BPSSpectrum`](@ref) and [`solve_tba`](@ref) iterates the fixed point on shared rapidity
grids, folding ``\pm\gamma`` pairs so that the kernel divergences cancel against each other.

This makes the integral-equation route the sharpest available check on the series pipeline:
two entirely independent computations of the same quantum period, agreeing to ``10^{-8}`` on
the cubic. It also reproduces the wall-crossing: the lateral jump of [`voros_value`](@ref)
across a ray is exactly the Delabaere-Dillinger-Pham factor, with the same integer Stokes
constant, obtained here as a residue rather than as a Borel discontinuity.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/tba.jl"]
```
