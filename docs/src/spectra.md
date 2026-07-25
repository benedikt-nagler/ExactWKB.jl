```@meta
CurrentModule = ExactWKB
```

# Spectra

Everything so far treats ``Q = V - E`` at one fixed energy. Turning the machinery into a
spectral method means letting ``E`` vary and asking which values make the Voros symbols
consistent - the exact quantization condition. The answer is not an asymptotic series but a
number, as accurate as the Borel–Padé summation behind it.

## Energy-dependent cycles

[`spectral_cycles`](@ref) finds the cycles of a real potential at energy ``E``: one `:well`
cycle per classically allowed interval (``V < E``) and one `:barrier` cycle per forbidden
interval between wells. [`quantum_period`](@ref) then re-runs the whole pipeline on
`with_energy(prob, E)` and returns the [`VorosSymbol`](@ref) of that cycle. Near an energy
where two turning points collide the construction is refused with
[`CoalescentTurningPoints`](@ref) rather than returning a number that quietly means nothing.

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
  exact at every order - which makes it the package's sharpest convention check.
- **symmetric double well**: the parity-factorized condition
  ``\cos\varphi = \tfrac{\sigma}{2}\sqrt{V_A/(1+V_A)}`` with ``\sigma = \text{parity}\cdot(-1)^n``,
  where ``V_A`` is the barrier (instanton) symbol. The two parities give the two members of a
  doublet, and their difference is the level splitting.

[`wkb_eigenvalue`](@ref) solves the condition by backtracking Newton from a bracketed
Bohr–Sommerfeld seed; [`spectral_determinant`](@ref) packages the same information as a
function whose zeros are the spectrum.

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
[`verify_zjj`](@ref) checks it order by order with ``c`` fitted, and returns residuals - for
``V = (z^2-1)^2`` the fit lands on ``c = 3/2`` and the residuals vanish to working precision.
[`energy_splitting`](@ref) is the even/odd gap of a doublet, which the same data predicts to
one-instanton accuracy.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/double_well.jl"]
```
