```@meta
CurrentModule = ExactWKB
```

# Seiberg-Witten geometry

Pure ``SU(2)`` ``\mathcal N = 2`` gauge theory is the exact Wentzel-Kramers-Brillouin problem
that made the subject famous. Its Nekrasov-Shatashvili limit is the Mathieu equation, which in
this package's
normalization ``\hbar^2\psi'' = (V-E)\psi`` reads

```math
V(x) = 2\Lambda^2 \cos 2x, \qquad E = u,
```

with ``u`` the Coulomb modulus. The two periods of the Seiberg–Witten curve are the classical
actions of that problem: the electric ``a(u)`` over a real cycle, the magnetic ``a_D(u)`` over
the tunnelling cycle through the barrier.

This layer computes them directly as period integrals rather than through the polynomial
recursion engine - ``\cos 2x`` is not a polynomial. It is deliberately firewalled:
`src/sw_curve.jl` imports nothing from the package beyond `Resurgence.FormalSeries` and
`QuadGK`, which keeps it extractable as a standalone period module.

A worked walk through this layer, with the physics spelled out, is
[Seiberg-Witten SU(2)](@ref).

## The u-plane

[`sw_periods`](@ref) gives both periods in closed form anywhere on the ``u``-plane, via
complete elliptic integrals evaluated through Carlson symmetric forms (so complex parameters
work, unlike real-only elliptic libraries):

```math
a(u) = \frac{2}{\pi}\sqrt{u+2\Lambda^2}\,E(\tilde m), \qquad
a_D(u) = -\frac{8i\Lambda}{\pi}\bigl[E(m) - (1-m)K(m)\bigr],
```

with ``\tilde m = 4\Lambda^2/(u+2\Lambda^2)`` and ``m = (2\Lambda^2-u)/4\Lambda^2``. The direct
quadratures are kept in the source as the oracle that pins these forms.

The curve degenerates at [`sw_singularities`](@ref) - the monopole point ``u = +2\Lambda^2``
where ``a_D \to 0``, and the dyon point ``u = -2\Lambda^2`` - plus the weak-coupling puncture
at infinity. Encircling any of them permutes the periods by an integer symplectic matrix
([`sw_monodromy`](@ref)); [`continue_periods`](@ref) transports them along an arbitrary path by
the Picard–Fuchs equation ``(u^2-4\Lambda^4)\Pi'' + \Pi/4 = 0``, which is a genuine analytic
continuation across the branch cuts and reproduces those integer matrices when the path closes.

The normalization of ``a_D`` (a factor 4 relative to the naive tunnelling action) is not a
convention one can choose freely: it is pinned simultaneously by integrality of the monodromy,
by the massless dyon at ``u = -2\Lambda^2`` carrying a lattice charge, and by metric positivity
``\mathrm{Im}\,\tau > 0``.

## Quantum periods

[`quantum_sw_periods`](@ref) adds the Nekrasov–Shatashvili ``\hbar^2`` corrections to *both*
periods. The magnetic cycle encloses turning points, where the standard Dunham integrand
``S_2 \sim Q^{-5/2}`` is not pointwise integrable - but no Voros regularization is needed.
Dropping a total derivative and using ``Q'^2 = 16\Lambda^4 - 4(Q+u)^2`` reduces everything to
``u``-derivatives of the classical action, and the Picard–Fuchs relation eliminates the third
derivative, leaving a universal first-order operator on *any* period ``\Pi``:

```math
a_2(\Pi) = \tfrac16 \Pi'(u) - \frac{u\,\Pi(u)}{12(u^2-4\Lambda^4)}.
```

Applied to ``a`` it reproduces the direct Dunham quadrature and its ``-\Lambda^4/(4u^{5/2})``
tail; applied to ``a_D`` it *is* the magnetic quantum period.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/sw_curve.jl"]
```

## States and wall-crossing

The Bogomolny-Prasad-Sommerfield (BPS) quiver of pure ``SU(2)`` is the Kronecker quiver - affine
``\tilde A_1``, hence cluster-*infinite*, which is the combinatorial reason the weak-coupling
spectrum is an infinite tower rather than a finite list. [`su2_bps_states`](@ref) returns the spectrum in either
chamber, resolved automatically from the wall of marginal stability
([`ms_wall`](@ref), [`sw_chamber`](@ref)):

- **strong coupling** (inside the wall): two states, monopole ``(1,0)`` and dyon ``(-1,2)``;
- **weak coupling** (outside): the W-boson ``(0,2)`` with ``\Omega = -2``, plus the
  preprojective and preinjective dyon towers ``(1,2n)`` and ``(-1,2n+2)``.

[`verify_su2_wall_crossing`](@ref) checks the Kontsevich–Soibelman identity relating the two as
an exact identity of truncated power series - the ``m=2`` analogue of the ``A_2`` pentagon. It
is self-validating: closure to a given degree simultaneously pins the spectrum content, the
``\Omega`` values including the vector multiplet's ``-2``, and the ordering. No external
"answer" is supplied.

```@autodocs
Modules = [ExactWKB]
Private = false
Pages = ["src/sw_bps.jl"]
```
