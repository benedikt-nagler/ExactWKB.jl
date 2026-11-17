# Archived, experimental, not validated

This directory preserves the working tree before the agreed reduction on 2026-10-26.
`original/` contains the original source, tests, examples, documentation. Planning recovery limits are recorded in `planning-recovery.md`.
The head revision, local patch, status and non-notebook SHA-256 digests record provenance.
Notebooks were copied without inspection. No archived file is loaded by the active package,
routine tests or documentation build. This is a historical snapshot, not a supported release.

## Removed capabilities

SU(2), its chambers and deformed walls, finite-radius GMN, hyperkahler metrics,
implicit GMN derivatives, motivic refinements, TBA, automatic BPS spectra,
triangulations and signed charge frames were removed. Double-well and uniform
quantization, spectral sensitivities and general merging-pair Weber mappings were removed.
Mixed files are preserved in full; the active copies retain only the agreed core.

## Reproduced faults

- At u = -1.999 + 0.001i the SU(2) spectrum reports masses 2.5444893706 and
  5.0894791389 while the charge (1,2) has mass 0.0007071289 in the same period frame.
  The hard-coded dyon charge is inconsistent with this branch.
- `su2_torus(sw,10,[0.5,0.9];chamber=:weak,tower=1)` assigns sigma = -1 to (0,2).
  This violates the quadratic-refinement identity on the declared integral lattice.
- The original square-contour period loses a sheet winding; root clustering also
  merges the two simple zeros of z^2-1e-10. These shared-core faults are repaired in
  active code, not retroactively changed in this snapshot.

## Limits

A repeated mutation word or a quantum-dilogarithm identity does not independently
validate a physical BPS spectrum. Fixed-point changes are not a total error bound.
Projected real symmetric metrics do not independently establish their raw reality or symmetry.
The original double-well conditions and general Weber transfer have not received a
complete derivation with a verified validity range. Existing numerical matches are retained
as evidence for those tested cases only. No complete revalidation of the removed theory was done.
