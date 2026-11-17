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

the Iwaki–Nakanishi normalization. It is fixed throughout: every turning point, period,
Voros symbol and Stokes line below is defined against this ``Q``. [`with_energy`](@ref)
re-tunes ``E`` without rebuilding the potential, which is what the energy-dependent
[Spectra](spectra.md) layer does internally.

Coefficients may be `Int`, `Rational` or any float type. Exact input flows through the
Wentzel-Kramers-Brillouin recursion in exact arithmetic and lands in `BigFloat` at whatever
precision the caller has set. `Float64` input stays in `Float64`, and there is no global
precision switch.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/potentials.jl"]
```

## Turning points

The zeros of ``Q``. [`turning_points`](@ref) finds them with `PolynomialRoots`, polishes them
by Newton, and assigns multiplicities by clustering, so a degenerate configuration is reported
as such and not silently returned as a pair of nearby simple roots. Most of the theory here
needs *simple* turning points, and the layers that do (Stokes graphs, triangulations) throw
[`UnsupportedTurningPoint`](@ref) instead of proceeding.

A double turning point is not an error in itself. It is what a symmetric double well has at
``E = 0``, and it becomes two simple points as soon as the energy moves off the well bottom.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/turning_points.jl"]
```

## Rational potentials

A [`RationalProblem`](@ref) carries ``Q = N(z)/\prod_j (z-p_j)^{m_j}``, the second problem
class the package understands. Its poles are held as **exact** data instead of being found
numerically, since a pole order is a topological datum of the surface the problem lives on and
not a numerical one. [`poles`](@ref), [`pole_orders`](@ref) and [`n_finite_poles`](@ref)
report them, and [`puncture_indices`](@ref) / [`n_punctures`](@ref) single out the *double*
poles, which behave differently from the rest.

An irregular pole (order ``\ge 3``) is a second escape boundary for a Stokes line, alongside
infinity. A double pole is not: ``\sqrt{Q} \sim \sqrt{c}/(z-p)`` makes ``\int\sqrt{Q}``
logarithmic, so trajectories spiral in and every ray that reaches it ends at the same place.
That is a *regular puncture*, a single point of the surface and not a boundary circle. Such a
puncture carries a wall of its own, belonging to no turning-point pair: at
``\theta \equiv \arg(2\pi i \sqrt{c})`` the spirals close into a ring domain, which
[`ring_domain_walls`](@ref) locates and [`RingDomainWall`](@ref) represents.

[`mathieu_problem`](@ref) is the Mathieu equation in this form. It lives here and not in a
periodic-potential layer of its own because substituting ``w = e^{2ix}`` turns Mathieu into a
*rational* ``Q`` with two order-3 poles. The periodic case and the rational case are one case,
and the same Stokes-graph machinery covers both.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/rational_potentials.jl"]
```
