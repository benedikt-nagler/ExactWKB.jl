```@meta
CurrentModule = ExactWKB
```

# Spectra

Everything so far treats ``Q = V - E`` at one fixed energy. Turning the machinery into a
spectral method means letting ``E`` vary and asking which values make the Voros symbols
consistent, which is the exact quantization condition. Its answer is a number and not an
asymptotic series, as accurate as the Borel–Padé summation behind it.

## Energy-dependent cycles

[`spectral_cycles`](@ref) finds the cycles of a real potential at energy ``E``: one `:well`
cycle per classically allowed interval (``V < E``) and one `:barrier` cycle per forbidden
interval between wells. [`quantum_period`](@ref) then re-runs the whole pipeline on
`with_energy(prob, E)` and returns the [`VorosSymbol`](@ref) of that cycle. Near an energy
where two turning points collide the construction is refused with
[`CoalescentTurningPoints`](@ref), since the alternative is a number that quietly means
nothing.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/quantum_periods.jl"]
```

## Exact quantization

[`quantization_condition`](@ref) is a real residual whose zero in ``E`` is an eigenvalue. Two
layouts are supported:

- **single well**: all-orders Bohr-Sommerfeld,
  ``\mathrm{Im}\log V_B^{\mathrm{med}} + 2\pi(n+\tfrac12) = 0``, where ``V_B`` is the
  median-summed well symbol. For the harmonic oscillator the WKB series truncates and this is
  exact at every order, which makes it the package's sharpest convention check.
- **symmetric double well**: the parity-factorized condition
  ``\cos\varphi = \tfrac{\sigma}{2}\sqrt{V_A/(1+V_A)}`` with ``\sigma = \text{parity}\cdot(-1)^n``,
  where ``V_A`` is the barrier (instanton) symbol. The two parities give the two members of a
  doublet, and their difference is the level splitting.

[`wkb_eigenvalue`](@ref) solves the condition by backtracking Newton from a bracketed
Bohr–Sommerfeld seed; [`spectral_determinant`](@ref) packages the same information as a
function whose zeros are the spectrum.

The Newton step uses an exact derivative and not a finite difference.
[`quantization_derivative`](@ref) returns the residual together with its derivatives with
respect to the energy or the potential coefficients, all from a single evaluation, so an
iteration costs one evaluation and carries no step-size parameter. Differentiating the
condition also answers a question the eigenvalue alone cannot:
[`eigenvalue_sensitivity`](@ref) gives ``\partial E_n/\partial v_k`` by the implicit function
theorem on ``F(E, v) = 0``, and returns ``\partial F/\partial E`` alongside it, which is the
conditioning of that answer.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/quantization.jl"]
```

## Double wells and the P/NP relation

The symmetric double well has two ``\hbar``-series: the perturbative one ``B(E,\hbar)``
(the well cycle) and the instanton one ``A(E,\hbar)`` (the barrier cycle). They are not
independent. The Zinn-Justin–Jentschura / Dunne–Ünsal relation

```math
\frac{1}{\partial B/\partial E} = -c\,\hbar^3\,
\left.\frac{\partial A}{\partial \hbar}\right|_{B}
```

determines the entire non-perturbative sector from perturbation theory alone.
[`verify_zjj`](@ref) checks it order by order with ``c`` fitted and returns the residuals. For
``V = (z^2-1)^2`` the fit lands on ``c = 3/2`` and the residuals vanish to working precision.
[`energy_splitting`](@ref) is the even/odd gap of a doublet, which the same data predicts to
one-instanton accuracy.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/double_well.jl"]
```

## The Weber model

A [`WeberModel`](@ref) is the parabolic-cylinder local model of a *merging* pair of turning
points, the situation an Airy-type local model cannot describe because it assumes one simple
turning point in isolation. [`weber_model`](@ref) builds it, [`weber_index`](@ref) is
``\nu = Z/(2\pi i \hbar)``, and [`weber_connection`](@ref) supplies the connection constant
that a merging pair contributes where a simple turning point contributes the Airy
``\tfrac12``.

The model comes with its own exact series: [`weber_voros_series`](@ref) is Stirling-type, and
[`weber_voros_sum`](@ref) Borel sums it, including on the imaginary axis, where the
singularities sit directly on the Laplace ray and a naive sum has nothing to integrate.
[`weber_barrier_amplitude`](@ref) computes the barrier transmission from the package's own
Borel-summed ``\Gamma``, which is what makes `uniform = true` in
[`quantization_condition`](@ref) work at and above the barrier top.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/weber.jl"]
```
