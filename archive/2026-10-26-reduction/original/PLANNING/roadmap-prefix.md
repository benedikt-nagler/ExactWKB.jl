# roadmap - ExactWKB.jl

Ecosystem context: `../../PLANNING/roadmap.md` (current priorities) and
`../../PLANNING/archive/research-program-roadmap.md` (historical design).

## Current scope and adoption gate (2026-09-19)

The supported core is a bounded Schrödinger/WKB workflow: potentials, expansions, periods,
Stokes geometry and quantization, extended by the cluster dictionary and wall-crossing.
Seiberg-Witten SU(2), refined identities and finite-radius GMN metrics remain specialized
applications with explicit validity conditions. They do not establish a general BPS or
hyperkähler platform. No repository split is scheduled.

**Coverage to make explicit before v0.1:** the generic Riccati/quantization pipeline is
polynomial; rational potentials have graph/topology support and specialized routes.
Regular punctures can appear in triangulations, but `charge_basis` and the downstream
signed-frame/BPS route exclude them. Exact quantization is limited to the implemented
problem classes and conditions. Enumerating these boundaries is not approval to lift them.

Choose one external researcher's calculation within the supported domain. Complete that
workflow before widening the theory. Replace private cross-package calls, separate routine
from extended numerical tests, and supply a reproduction against an independent calculation.
The support review must distinguish formal expansions, numerical agreement and certification.

**Third-order WKB and higher-rank GMN expansion are deferred.** Reopen only for a concrete
external calculation after the current workflows are supported. Require a bounded target
and an independent oracle, including a published reference calculation where applicable.
Neither a ClusterML preprint nor a new theorem is a software prerequisite. Research novelty
claims remain subject to the separate novelty gate.

**Exit criterion:** an agreed support list, tagged v0.1, clean-environment installation,
routine CI within ten minutes, public cross-package interfaces, and an external user who
can reproduce and interpret the selected calculation.

## Geometry integration boundary (2026-09-19)

The [geometry track](../../PLANNING/decisions/geometry-computation.md) develops general
local periods outside ExactWKB.jl. Preserve the specialized SU(2) calculations as reference
cases. Its known monodromy matrices are not a general reconstruction algorithm.

When a concrete geometry consumer needs it, design a bounded adapter for externally supplied
periods or central charges. Match bases, pairings and normalizations explicitly. Require
agreement with existing SU(2) results before replacing an internal path. Preserve the public
representation boundary and avoid a new dependency until the adapter design justifies it.
Supplied central charges do not alone authorize arbitrary BPS-spectrum inference.

General period work does not wait for the WKB external trial or v0.1 release. Those remain
parallel support obligations. Further WKB/BPS theory expansion stays behind the new geometry
priority unless a separate concrete calculation justifies it. No code extraction is scheduled.

## Milestones

### M3 - Schrödinger / Stokes-graph core  ⟵ DONE (2026-07-12)

Scaffold + the exact-WKB layer, no cluster bridge yet. Files, in implementation order
(all **done**):

1. scaffold + Aqua-only suite
2. `errors.jl` + `potentials.jl`
3. `turning_points.jl`
4. `wkb_recursion.jl`
5. `periods.jl`
6. `voros.jl`
7. `stokes_graph.jl`
8. `saddles.jl`
9. `show.jl` + Makie extension
10. precompile, `examples/m3_flagship.jl` end-to-end, full suite

Full suite green (Aqua + all `test_*.jl`, incl. headless CairoMakie). The flagship
double-well chain (turning points → WKB → Voros → Resurgence Borel/Padé → Stokes graph
with the instanton saddle → PNG) runs from `examples/m3_flagship.jl`. See
`RECENT_CHANGES.md` for the per-file log until it is folded in here.

### DDP layer - the convention-pinned connection formula  ⟵ DONE (2026-07-14)

Ecosystem step 0 item 4 (`src/ddp.jl` + `test/test_ddp.jl` + `examples/ddp_demo.jl`):
Voros wall-crossing = cluster y-mutation, every sign pinned by an exact oracle on the
cubic A₂ potential (see the convention ledger in `src/ddp.jl` and its summary in
`ARCHITECTURE.md`). `voros_value` (lateral/median sums via Resurgence),
`verify_ddp` (measured jump = ⟨γ,γ'⟩·log(1+V'), Stokes constant an exact integer),
`intersection_pairing` (canonical-lift crossing count), `ddp_seed` (B = −P),
`verify_ddp_mutation` (jump = FZ μ_k, opposite direction fails by ≥ 10×), pentagon /
DT closure against ClusterAlgebras.jl's `y_system`/`dt_transformation`/`omega`. The `ClusterAlgebras`
dependency is live (Aqua stale-dep ignore removed).

### M4 - the cluster bridge  ⟵ DONE (2026-07-16)

The Iwaki–Nakanishi dictionary, fully automatic (see `RECENT_CHANGES.md` for the
per-file log until folded in here): `src/triangulation.jl` (generic Stokes graph →
ideal triangulation of the (d+2)-gon by a combinatorial planar face walk →
`triangulation_quiver`, finite type A_{d−1} in every chamber; flip = mutation
oracles), `src/charge_lattice.jl` (`charge_basis`: validated decay-oriented cycles,
central charges, intersection pairing; `bridge_seed`), `src/bps.jl` (`bps_spectrum`:
greedy phase-ordered maximal green sequence, gated by the confirmed-saddle
cross-check), `plot_triangulation`, `examples/m4_flagship.jl` (cubic A₂ + quartic A₃
+ quintic A₄ stress test; Borel poles on the ±Z_γ lattice; integer Stokes constants;
no hand-built input). Full suite 542 green.

**Scope boundary (record; CORRECTED 2026-07-31 by M7a):** degree-d polynomial
potentials with simple turning points triangulate a (d+2)-gon, so the reachable cluster
types are exactly A_{d−1}. Non-A types need **punctures**, i.e. *poles* of Q - a double
pole is a puncture (D-type, the (A₁,D_n) theories), a simple pole an orbifold point
([IN16]). This entry used to add "or higher-order turning points (the deferred Weber
local models)", which is **wrong**: the marked-point count is d+2 whatever the zero
multiplicities are, so the surface stays a disk, and an order-m turning point's dual
cell is an (m+2)-gon - a *degenerate* triangulation of the same (d+2)-gon, a
positive-dimensional face of the A_{d−1} associahedron, not a new type. The route to
non-A types is the rational-potential item alone.

**M4 findings (2026-07-16) - CLOSED by M4b, see below.** The per-chamber
cycle-orientation rule beyond the DDP all-decay frame was not implemented, so
(1) *cyclic chambers* (a triangulation triangle with three internal edges, generic
for d ≥ 5) had no all-decay seed and threw `ChamberError`; (2) some tree chambers
(quintic real-axis fan) swept a cluster-consistent but unphysical green sequence,
caught only by the saddle gate; and `bps_spectrum` had to *search* chambers for one
that worked (empirically the top chamber).

### M4b - the signed-flip layer  ⟵ DONE (2026-07-19)

`src/signed_frame.jl` + `src/triangulation.jl`'s new `flip`/`canonical_reorder`.
The decay rule fixes a *representative* of ±γ_e; the physical cycle is
`γ_e^phys = ε_e·γ_e`, and ε is solved by inverting the **tightened keystone**
`P[i,j] == −ε_i ε_j·B[i,j]` (M4 checked only `|P| == |B|` - that `abs` was the slack
hiding the layer). Solving is a spanning-tree propagation over the quiver graph;
solvability needs the sign cocycle to close around every cycle, which is vacuous on
tree chambers and a *genuine* condition on cyclic ones - the one M4 could not meet.
`ChargeBasis` gained `signs`/`physical_charges`/`signed_pairing`; `contours` and
`central_charges` stay in the decay frame, so ledger item 4 and the M5 orientation
ledger are untouched.

A second, independent defect was fixed at the same time: **ledger item 6 is
chamber-relative**. `_sweep_chamber` ordered walls by *absolute* phase, which is only
right at the top of `[0,π)`; the general rule is increasing `mod(θ₀ − θ_c, π)`, of
which the absolute form is the θ₀ → π⁻ case. This alone explains why M4 enumerated
correctly in the top chamber and misordered elsewhere.

Results: every chamber of cubic/quartic/quintic now yields the physical spectrum -
7 states in all six quintic chambers, including the four *cyclic* ones and the fan
chamber (previously 6 unphysical states). The chamber *search* in `bps_spectrum` is
gone (default θ = widest wall-free gap, one attempt); `_matches_saddles` changed
meaning from chamber filter to verification of a frame correct by construction.

**M4b findings (record):**
- **`flip` is handed.** Not a free binary labelling choice: a fixed convention matched
  the tracer 6/6 crossing `lo→hi` and 0/6 crossing `hi→lo`, and `flip∘flip ≠ id`.
  `flip(t, k; direction = ±1)` - the Stokes graph rotates with θ, so the collapsing
  strip reconnects oppositely by crossing sense.
- **ε is chamber-local; there is no transport theorem.** It is *not* continuous across
  a wall even on diagonals whose turning-point pair - hence whose cycle - persists
  (cubic θ ≈ 0.785): the cycle persists but the seed is a different object per
  chamber. The residual gauge (one global sign per component) is physically inert and
  set by convention. What the flip walk proves, and all `verify_signed_frame`
  asserts, is that it lands on exactly the traced triangulation (12/12 chambers).
- **Wall phases live mod π.** The quintic's saddles at θ_c = 0 `(1,2)` and θ_c = π
  `(3,4)` are one degenerate wall; not reducing makes the top chamber look 4.7e-11
  wide.
- **Reference chamber = a pin, not a derivation.** The top chamber `(θ_max, π)` has
  ε ≡ +1 (verified in all three fixtures); it is the gauge in which ledger items 4/5/6
  were originally fixed. A potential whose top saddle sits at π⁻ has no traceable
  reference chamber - `reference_theta` throws (recorded scope boundary).
