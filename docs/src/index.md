```@meta
CurrentModule = ExactWKB
```

# ExactWKB.jl

Exact Wentzel-Kramers-Brillouin analysis of one-dimensional Schrödinger equations with
polynomial potentials - equivalently, of second-order linear differential equations with a
small parameter: turning points, the all-orders expansion, quantum periods and Voros symbols,
Stokes graphs, exact quantization conditions, spectra of stable states and their wall-crossing,
and the dictionary to cluster algebras.

The object is the Schrödinger equation in its complex-analytic form

```math
\hbar^2 \psi''(z) = Q(z)\,\psi(z), \qquad Q(z) = V(z) - E,
```

with ``Q`` a polynomial and ``z`` complex. The
[Wentzel-Kramers-Brillouin ansatz](https://en.wikipedia.org/wiki/WKB_approximation) turns it
into a Riccati equation solved by a power series in ``\hbar`` that diverges for every
``\hbar \neq 0``. **Exact** Wentzel-Kramers-Brillouin analysis makes it usable anyway: Borel
summed, the series gives genuine solutions rather than an approximation, and the *Stokes graph*
- the level set ``\mathrm{Im}\left(e^{-i\theta}\int\sqrt{Q}\right) = 0`` - records which
solution you are on, jumping discontinuously as the summation direction ``\theta`` rotates.

That last jump is where cluster algebras enter. By the Iwaki-Nakanishi dictionary a generic
Stokes graph is dual to an ideal triangulation, the jump of the Voros symbols across a wall is
a cluster ``y``-mutation, and a full circuit of the ``\theta``-plane is a maximal green
sequence. The package computes both sides and checks that they agree - the combinatorics is
exact, so it certifies the numerics.

**Scope.** ``Q`` must be polynomial: the harmonic, quartic, double-well and cubic oscillators
and everything else usually used to test the theory, with the Seiberg-Witten layer adding the
Mathieu problem through its own period module. Exact `Rational` input flows through to
`BigFloat` at whatever precision the caller has set.

## Installation

Not yet registered in General, but both of its foundations are, so Pkg resolves them by
name. From the Julia REPL:

```julia
pkg> add https://github.com/benedikt-nagler/ExactWKB.jl
```

Requires Julia 1.10 or later.

## Quick start

```julia
using ExactWKB

# Double well, Q = (z² − 1)² − 1/4, coefficients in ascending order
prob = SchrodingerProblem([0.75, 0.0, -2.0, 0.0, 1.0])

tps = turning_points(prob)           # four simple turning points
w   = wkb_expansion(prob; order = 6) # Riccati recursion for S_odd
vs  = voros_symbol(w, encircling_contour(tps[2], tps[3]))

classical_period(vs)                 # the ħ⁻¹ term, ∮√Q ≈ 1.8399
g = stokes_graph(prob; theta = 0.3)  # traced Stokes lines
saddles(prob)                        # 3 saddle connections, with central charges
```

Exact quantization gives energy levels rather than just series:

```julia
harm = SchrodingerProblem([0.0, 0.0, 1.0])   # V = z², so Eₙ = 2ħ(n + ½)
wkb_eigenvalue(harm, 0, 0.1; order = 6)      # 0.1
wkb_eigenvalue(harm, 3, 0.1; order = 6)      # 0.7
```

And the cluster bridge, on the cubic:

```julia
cubic = SchrodingerProblem([0.0, -1.0, 0.0, 1.0])  # Q = z³ − z
t = ideal_triangulation(stokes_graph(cubic; theta = 0.3))
q = triangulation_quiver(t)          # a ClusterAlgebras.Quiver - here type A₂

sp = bps_spectrum(cubic; theta = 0.3, margin = 0.4)
charges(sp)                          # [[1, 0], [0, 1]] - charges in phase order
solve_tba(sp)                        # the same periods again, from the spectrum alone
```

## Tutorials

Three worked problems, each stated as a question and answered end to end in code. They are the
best way in: together they exercise every layer of the package and show what each check is
worth.

**[The cubic, end to end](@ref)** - *what is a divergent series hiding, and can we recover the
exact answer from it?* One cubic potential from turning points to the Borel plane, the Stokes
graph, a pentagon triangulation, the ``A_2`` quiver, the spectrum as a maximal green sequence,
wall-crossing as a cluster mutation measured to seven digits, and finally the same periods
recomputed from the spectrum alone through integral equations - agreeing to six digits with the
resummation, and producing the golden ratio in the infrared.

**[The double well](@ref)** - *how do you compute a quantity that is zero to all orders of the
only expansion you have?* The tunnelling splitting of a symmetric double well: the instanton as
a visible edge of the Stokes graph, a gap of ``6\times10^{-6}`` from a resummed quantization
condition, agreement with the one-instanton estimate, and the exact relation that determines
the instanton series from perturbation theory, verified to twenty digits.

**[Seiberg-Witten SU(2)](@ref)** - *which particles exist in a strongly coupled gauge theory,
and how does the list change with the vacuum?* Periods on the whole moduli space checked
against a three-instanton expansion, monodromy measured by transporting a differential equation
and coming out as exact integers, the wall of marginal stability, the two spectra it separates,
and a Kontsevich-Soibelman identity closing to exactly zero in rational arithmetic.

## Manual

| Page | Covers |
|------|--------|
| [Problems and turning points](problems.md) | `SchrodingerProblem`, energies, roots of ``Q`` and their multiplicities |
| [The WKB expansion](wkb.md) | the Riccati recursion, period integrals with branch tracking, Voros symbols |
| [Stokes graphs](stokes.md) | the augmented-state tracer, topology signatures, saddle connections |
| [Spectra](spectra.md) | energy-dependent cycles, exact quantization, eigenvalues, double-well splittings |
| [Wall-crossing](ddp.md) | lateral Borel sums, the Delabaere-Dillinger-Pham jump, wall-crossing as ``y``-mutation |
| [The cluster bridge](bridge.md) | ideal triangulations, quivers, charge lattices, signed frames, spectra of stable states |
| [Thermodynamic Bethe ansatz](tba.md) | the conformal-limit Gaiotto-Moore-Neitzke equations: quantum periods from the spectrum alone |
| [Seiberg-Witten geometry](seiberg_witten.md) | ``SU(2)`` periods on the ``u``-plane, monodromy, the wall, wall-crossing |
| [Plotting](plotting.md) | the Makie extension |

The [API index](@ref) lists every exported name.

## Related packages

This is the bridge of a family of Julia packages for **exact and asymptotic methods** - where a
discrete, exactly computable structure controls a continuous, only-asymptotically-defined one.
[Resurgence.jl](https://github.com/benedikt-nagler/Resurgence.jl) is the continuous foundation,
providing the Borel-Padé summation, transseries and alien calculus used here;
[ClusterAlgebras.jl](https://github.com/benedikt-nagler/ClusterAlgebras.jl) is the discrete one,
providing quivers, seed mutation, green sequences and Donaldson-Thomas invariants. Neither
depends on the other, or on this package.

Names from those packages are **not** re-exported: use `import Resurgence` /
`import ClusterAlgebras` and qualify. Public types here expose plain numbers and
`Resurgence.FormalSeries`; AbstractAlgebra ring elements never leak out of the recursion.
