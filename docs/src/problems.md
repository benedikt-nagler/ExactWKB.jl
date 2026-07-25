```@meta
CurrentModule = ExactWKB
```

# Problems and turning points

Everything in the package starts from a [`SchrodingerProblem`](@ref): a polynomial potential
``V`` given by its coefficients in ascending order, together with an energy ``E``. The object
of the theory is

```math
\hbar^2 \psi''(z) = Q(z)\,\psi(z), \qquad Q(z) = V(z) - E,
```

the Iwaki–Nakanishi normalization, and it is fixed throughout - every turning point, period,
Voros symbol and Stokes line below is defined against this ``Q``. [`with_energy`](@ref)
re-tunes ``E`` without rebuilding the potential, which is what the energy-dependent
[Spectra](spectra.md) layer does internally.

Coefficients may be `Int`, `Rational` or any float type. Exact input flows through the
Wentzel-Kramers-Brillouin recursion in exact arithmetic and lands in `BigFloat` at whatever
precision the caller has set; `Float64` input stays in `Float64`. There is no global precision
switch.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/potentials.jl"]
```

## Turning points

The zeros of ``Q``. [`turning_points`](@ref) finds them with `PolynomialRoots`, polishes them
by Newton, and assigns multiplicities by clustering, so a degenerate configuration is reported
as such rather than silently returned as a pair of nearby simple roots. Most of the theory
here needs *simple* turning points; the layers that do (Stokes graphs, triangulations) throw
[`UnsupportedTurningPoint`](@ref) rather than proceed.

A double turning point is not an error in itself - it is what a symmetric double well has at
``E = 0``, and it becomes two simple points as soon as the energy moves off the well bottom.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/turning_points.jl"]
```
