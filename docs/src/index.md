# ExactWKB.jl

Compute finite Riccati expansions, numerical periods and selected Stokes graphs for
`hbar^2 * psi'' = (V-E) * psi`. Numerical output is not a certified Borel sum.

Real-coefficient polynomial potentials support the formal WKB recursion.
Complex evaluation points are supported; nonreal coefficients are refused. Rational potentials support
classical periods and numerical trajectory tracing only. Pole orders are supplied data.
A root multiplicity belongs to the exact polynomial represented by the input numbers.
The root locations remain numerical.

Quantization supports confining quadratics and symmetric quartics
`V0 + a*z^2 + b*z^4`, with `a >= 0` and `b > 0`. Harmonic eigenvalues are exact.
Anharmonic phase conditions are finite-order approximations and may miss
nonperturbative contributions. A small root residual does not bound their spectral error.

Voros summation samples explicit rays `theta +/- tilt`; stability in order and tilt
must be checked separately. A finite approximant does not certify physical laterals.

The DDP comparison takes supplied cycles and intersection data. The jump formula
requires simple zeros and one isolated regular saddle (Iwaki-Nakanishi 2014, Theorem 3.4).
The routines do not establish these hypotheses or infer a complete BPS spectrum.

The Weber routines refer only to the quadratic normal form and its gamma asymptotics.
Stokes trajectory endpoints use numerical event tolerances and finite integration budgets.
Incomplete tracing raises `TracingFailed` unless partial output is explicitly requested.

Removed applications are preserved in the repository's `archive/2026-10-26-reduction/`.
They are experimental and unvalidated, and are excluded from loading, tests and this manual.
