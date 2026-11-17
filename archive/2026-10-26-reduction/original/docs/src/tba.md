```@meta
CurrentModule = ExactWKB
```

# Thermodynamic Bethe ansatz

The spectrum of stable (Bogomolny-Prasad-Sommerfield) states is a finite list of integers and
complex numbers: ``(Z_\gamma, \Omega(\gamma), \langle\gamma,\gamma'\rangle)``. That list is
enough to reconstruct the Borel-summed Voros symbols exactly, with no
Wentzel-Kramers-Brillouin series and no Padé anywhere in the pipeline. The mechanism is the
conformal limit of the Gaiotto–Moore–Neitzke construction: the symbols
``\mathcal X_\gamma(\hbar)`` solve a system of integral equations of
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

The integral-equation route is therefore the sharpest available check on the series pipeline:
two independent computations of the same quantum period, agreeing to ``10^{-8}`` on the cubic.
It also reproduces the wall-crossing. The lateral jump of [`voros_value`](@ref) across a ray is
exactly the Delabaere-Dillinger-Pham factor with the same integer Stokes constant, obtained
here as a residue instead of as a Borel discontinuity.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/tba.jl"]
```

## Hyperkähler metrics

The finite-radius companion of the TBA layer, following Gaiotto–Moore–Neitzke's twistor
construction. Where `solve_tba` takes the conformal limit, [`solve_gmn`](@ref) keeps the
radius ``R`` finite and solves for the twistor coordinates ``\mathcal{X}_\gamma``
([`xi_value`](@ref)), whose deviation from [`semiflat_xi`](@ref) is the instanton correction.

[`gmn_torus`](@ref) packages the input, and it takes **lattice-level BPS data only**: charges,
``\Omega``, ``Z`` and a Gram matrix. It therefore serves any source of a BPS spectrum, whether
traced from a potential or supplied in closed form by [`su2_torus`](@ref). From a solved
[`MetricPoint`](@ref), [`symplectic_expansion`](@ref) fits the holomorphic symplectic form
``\varpi(\zeta) = -\tfrac{i}{2\zeta}\omega_+ + \omega_3 - \tfrac{i}{2}\zeta\omega_-``
and [`hk_metric`](@ref) assembles ``g = \omega_1\omega_3^{-1}\omega_2``. That ``\varpi`` has
exactly this Laurent degree *is* the twistor theorem, so the residual of the fit is the theorem
run as a numerical test. [`semiflat_metric`](@ref) is the uncorrected comparison.

The Ooguri–Vafa case of a single BPS state is solved a second, independent way in closed form
([`ooguri_vafa_torus`](@ref), [`ooguri_vafa_xi`](@ref), [`ooguri_vafa_instantons`](@ref) as a
Bessel-``K`` sum), because with one state the integral equation degenerates into a plain
integral. Two solvers reaching the same metric is the strongest check in the layer.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/hyperkahler.jl"]
```

## Moduli derivatives

[`solve_gmn_derivative`](@ref) differentiates a converged [`GMNSolution`](@ref) with
respect to the moduli, by implicit differentiation of the fixed point rather than a finite
difference of [`solve_gmn`](@ref). The nine solver calls a finite difference of
[`metric_point`](@ref) needs collapse into one linear solve, so the accuracy of the metric
stops being set by a step size. The input `dZdu` is the holomorphic ``u``-derivative of the
*basis* central charges. State charges follow by linearity.

The moduli are ordered ``(\operatorname{Re} u, \operatorname{Im} u, \theta_1, \dots,
\theta_r)``, and [`n_parameters`](@ref) counts them. [`log_xi_derivative`](@ref) evaluates
``\partial_\mu \log \mathcal{X}_\gamma(\zeta)`` at fixed ``\zeta``. The BPS data is held
constant, so the derivative is undefined on a wall of marginal stability.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/gmn_derivative.jl"]
```
