# ExactWKB.jl

Compute finite WKB expansions, numerical periods, Voros sums and Stokes trajectories
for selected one-dimensional Schrödinger problems. The package is unregistered.

## Installation

```julia
using Pkg
Pkg.develop(path="/path/to/ExactWKB.jl")
```

```julia
using ExactWKB
p = SchrodingerProblem([-1.0, 0.0, 1.0])
w = wkb_expansion(p; order=4)
v = voros_symbol(w, encircling_contour(turning_points(p)...))
classical_period(v) # approximately -pi*im
```

- **Formal expansions:** real polynomial potentials, finite Riccati coefficients and their parameter derivatives.
- **Periods:** continued square roots along supplied piecewise linear contours.
- **Summation:** finite Borel-Padé approximations of Voros tails on explicit sector rays.
- **Stokes graphs:** numerical trajectories with event tolerances and finite tracing budgets.
- **Quantization:** exact harmonic results and approximate single-well quadratic/quartic calculations.
- **DDP comparison:** supplied cycles and pairings; no automatic BPS inference.
- **Weber reference:** the quadratic normal form and its gamma asymptotics.

Rational potentials support classical periods and tracing, not the polynomial Riccati pipeline.
Anharmonic root residuals are not spectral error bounds. Multiple zeros require algebraic evidence.

See the [manual](https://benedikt-nagler.github.io/ExactWKB.jl/) for supported domains.
Historical applications are [archived](archive/2026-10-26-reduction/STATUS.md) and unvalidated.
Related packages: Resurgence.jl and ClusterAlgebras.jl. See [LICENSE](LICENSE).
