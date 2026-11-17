# Performance & computational bottlenecks - ExactWKB.jl

A running record of where the compute goes, why, what has already been mitigated, and
what to try next. **Do not optimize speculatively** - every entry here is something
that actually bit a real calculation (most of them the M5 double-well spectroscopy,
2026-07-17). Reach for these when a run is slow, not before.

Rule of thumb from M5: the exact-WKB spectral pipeline is
`wkb_eigenvalue → (Newton) → quantization_condition → voros_symbol → wkb_period
(quadrature) + voros_value (Borel–Padé–Laplace)`, all in BigFloat. Cost multiplies
along that chain, so a fix low in the chain (a single quadrature) pays off ~100×
upstream (a full splitting solve).

## Ranked bottlenecks

### 1. Even-period quadrature runaway on barrier cycles - FIXED, keep in mind

**Symptom (M5):** one barrier `voros_symbol` took **26 s** at order 6 / 120-bit while
the well symbol took 1 s; a single order-6 double-well eigenvalue was minutes.

**Root cause:** `voros_symbol` computes the even periods `v₂, v₄, …` only to assert
they are ≈ 0 (a branch-error sanity check). The high-order integrand `pₘ(z)·√Q/Qᵏ`
has enormous, sharp peaks near the turning points (Q → 0), magnitude O(1e6+), while the
integral itself cancels to ≈ 0. With the atol tied to the (small) classical period, the
adaptive Gauss–Kronrod quadrature kept subdividing the peaks trying to resolve zero to
an impossible relative precision, exhausting `maxevals = 10⁶`.

**Mitigation applied** (`src/voros.jl`): compute the odd periods first, track their
magnitude `oddscale`; give the even periods a coarse tolerance (`atol = 1e-3·(1+oddscale)`,
`rtol = 1e-2`) **and a capped `maxevals = 10⁴`**. Barrier symbol 26 s → 1.4 s; harmonic
exactness preserved (tail < 1e-37). A branch error still shows up as an O(oddscale) even
period, so the check still bites.

**Future:** the even-period check is pure insurance. Consider skipping it entirely above
a trusted order, or replacing the quadrature check with the algebraic identity
`S_even = −½ ∂_z log S_odd` evaluated at one point (already available as
`even_odd_residual`) - no contour integral at all.

### 2. Adaptive quadrature order vs. BigFloat precision - PARTLY FIXED

**Symptom (M5):** a stiff order-15 quartic period took **7.8 s** at the default
Gauss–Kronrod order 7; the same integral at order 21 took **0.22 s** (35×), identical
value.

**Root cause:** the near-turning-point `1/Qᵏ` integrands are analytic on the contour but
sharply peaked; low-order adaptive quadrature pays for that with deep bisection, and at
BigFloat precision each function evaluation is expensive, so bisection depth dominates.
High fixed order converges exponentially and avoids the subdivision.

**Mitigation applied** (`src/periods.jl` `_quad_order`): Gauss order scales with
precision - 7 for Float64, `clamp(precision/8, 21, 61)` for BigFloat.

**Future:** the order heuristic is coarse. The true cost driver is the closest approach
of the contour to a turning point (peak sharpness ∝ 1/distance). A contour-aware order
(or a turning-point-avoiding contour shape / a change of variables that regularises the
`√Q` endpoint behaviour, e.g. a trig substitution around each turning point) would beat a
flat order. This is the single highest-leverage quadrature improvement.

### 3. Every Newton iterate re-runs the whole WKB pipeline - OPEN (flagged in the M5 plan)

**Symptom (M5):** a double-well splitting at small g (≈ 0.03) cost **70–90 s** - two
parity eigenvalue solves, each a finite-difference Newton of ~10–15 iterations.

**Root cause:** `wkb_eigenvalue` does finite-difference Newton, so each iteration
evaluates `quantization_condition` at **three** energies (`E`, `E±h`), and each of those
rebuilds the WKB expansion and re-integrates every cycle's Voros symbol from scratch at a
new `E`. Nothing is cached between iterates, and the ± evaluations double the count. For
the double well each condition eval also sums two cycles (well + barrier). So one
splitting ≈ (2 parities) × (~12 iters) × (3 evals) × (well + barrier Voros + median sum).

**Mitigations already in place:** within a *single* `quantization_condition` call,
`_condition_data` builds one WKB expansion and reuses it for the well and barrier
(`src/quantization.jl`). The Newton step is floor-clamped and backtracks out of invalid
layouts (robustness, not speed).

**Future (highest structural leverage):**
- **Analytic-ish E-derivative.** The condition is `Im log V_B + …`; its E-derivative is a
  period whose integrand is the E-derivative of `S_odd`. Computing `∂/∂E` by differentiating
  the integrand (one extra period per cycle) instead of two extra full condition evals cuts
  Newton cost by ~⅓ and removes the finite-difference step-size fragility.
- **Cache/warm-start across iterates.** Newton energies cluster; the WKB expansion depends
  smoothly on E. A per-(prob) expansion cache keyed on E, or a low-order predictor for the
  Borel–Padé data, would avoid rebuilding from scratch each step.
- **Better seed.** `_seed_energy` uses the leading (classical) Bohr–Sommerfeld action.
  Seeding with the median-summed `B = n + ½` root instead would land Newton within one or
  two steps of the answer and avoid the oscillatory-condition backtracking that small-ħ
  double wells trigger.

### 4. Borel–Padé–Laplace median sum degrades as ħ → 0 - OPEN, intrinsic to small ħ

**Symptom (M5):** the ZJJ one-instanton splitting law needs the ground state *below the
barrier top*, i.e. small g / small ħ (g = 0.08 already fails: E₀ > barrier). But smaller
ħ makes every `voros_value` (the lateral/median Laplace sum) more expensive, because the
Laplace integral must be carried further out along the Borel ray and the summation's
adaptive quadrature refines more. This is the deep reason double-well spectroscopy is
costly: the physically interesting regime (deep tunneling) is the numerically hard one.

**Root cause:** exponential accuracy of the resummation at small ħ costs proportionally
more Borel-plane work; this lives in `Resurgence.lateral_sum`, not here.

**Future:** this is the designated territory of **Resurgence M8** (Arb rigor, fast
summation kernels) and the conformal/hyperasymptotic tools of Resurgence M5. Do not
optimise the sum inside ExactWKB - feed it better (fewer, cheaper calls per §3) and lean
on Resurgence for the kernel. A hyperasymptotic (optimal-truncation) shortcut for the
leading small-ħ behaviour may sidestep the full Laplace integral where only the
one-instanton prefactor is wanted.

### 5. BigFloat WKB (Riccati) recursion at high order - OPEN (known, ecosystem-level)

**Symptom:** order-14–16 BigFloat expansions are a noticeable fraction of each condition
eval; the exact-`QQ` ring blows up faster still (numerator degrees grow with order).

**Root cause:** documented in `ARCHITECTURE.md` (“The Riccati representation”): the
numerator polynomial degree grows with order; exact arithmetic compounds it. This is the
ecosystem's known cost driver, with **Resurgence M8 fast kernels** the planned relief.

**Future:** for spectral work the median sum tolerates modest orders - prefer *lower* WKB
order + accurate summation over high order (the M5 plan's explicit guidance: “drop the
recursion order and lean on median-summation accuracy before reaching for performance
work”). Order 10 already gave 1e-10 double-well eigenvalues; order 14–16 was overkill and
was walked back in the tests.

### 6. Redundant work in the seed bisection and the verify_zjj stencil - MINOR

- `_seed_energy` runs a growth loop + up to 80 bisection steps, each calling
  `turning_points` (BigFloat root-finding) and a `period_integral`. Cheap per call
  (~0.5 s total) and dwarfed by §3, but it recomputes turning points it could cache.
- `verify_zjj` evaluates `B` and `A` on a **5-point** E-stencil (5 WKB expansions + 10
  Voros symbols) for a 4th-order finite difference. Fine now that §1/§2 are fixed, but the
  analytic E-derivative of §3 would replace the stencil with a single evaluation.

## Cross-cutting principles (for future perf work here)

1. **Fix low in the chain.** A quadrature that is 10× faster makes a splitting solve ~10×
   faster; the opposite direction (optimising the Newton loop) helps far less if each eval
   is still slow. §1 and §2 were the right first targets for that reason.
2. **The turning points are the enemy.** Every quadrature bottleneck here traces to the
   `1/Qᵏ` peak at a turning point. Contour shaping and endpoint-regularising substitutions
   attack all of them at once.
3. **Prefer accurate summation over high WKB order.** Cost grows fast in order; the
   median sum is forgiving. This is a genuine free lunch for spectra.
4. **Small ħ is intrinsically hard** and is Resurgence's problem, not ours - don't
   reimplement summation kernels in ExactWKB.
5. **Measurement note (workflow, not code):** Julia buffers stdout to a pipe, so
   background-run progress is invisible until completion. When profiling, run in the
   foreground or `flush(stdout)` after each print - this cost real debugging time in M5.

## Concrete “when we come back to this” shortlist

In rough order of leverage for double-well / spectral runs:
1. Analytic (integrand-differentiated) E-derivative in `wkb_eigenvalue`/`verify_zjj` (§3, §6).
2. Turning-point-aware contour shape or endpoint substitution for the period quadratures (§2).
3. Median-summed `B = n+½` seed for the Newton solve (§3).
4. WKB-expansion cache / warm start across Newton iterates (§3).
5. Skip or algebraically replace the even-period sanity check above a trusted order (§1).
