# References

Citation keys used in `src/` comments and in the other `PLANNING/` files. Cite as
`[KEY] §n` where the mathematics is not textbook-standard; the annotation in
parentheses records what this package actually takes from the source.

Conventions that are *ours* rather than taken from a source live in the convention
ledger at the top of the owning `src/` file (`sw_curve.jl`, `double_well.jl`,
`quantization.jl`, `sw_bps.jl`, `signed_frame.jl`), each item naming the oracle that
pins it.

## Exact WKB

- **[Voros83]** A. Voros. *The return of the quartic oscillator: the complex WKB
  method.* Ann. Inst. H. Poincaré A 39 (1983) 211-338. (Voros symbols / quantum
  periods, the exact quantization condition, the quartic as the model problem)
- **[DDP97]** E. Delabaere, H. Dillinger, F. Pham. *Exact semiclassical expansions
  for one-dimensional quantum oscillators.* J. Math. Phys. 38 (1997) 6126-6184.
  (the DDP discontinuity formula for Voros symbols - the jump rule implemented in
  the DDP layer; median summation of the Voros multiplier)
- **[KawaiTakei]** T. Kawai, Y. Takei. *Algebraic analysis of singular perturbation
  theory.* AMS, 2005. (Stokes graphs, turning points, the odd/even split of the
  Riccati solution used by `wkb_recursion.jl`)
- **[BNR82]** H. L. Berk, W. M. Nevins, K. V. Roberts. *New Stokes' lines in WKB
  theory.* J. Math. Phys. 23 (1982) 988-1002. (higher-order / merging turning points
  and the new Stokes curves they emit - the setting of `weber.jl` and of the `m+2`-ray
  tracer in `stokes_graph.jl`)
- **[Takei08]** Y. Takei. *Sato's conjecture for the Weber equation and transformation
  theory for Schrödinger equations with a merging pair of turning points.* RIMS
  Kôkyûroku Bessatsu B10 (2008) 205-224. (**the** reference for `weber.jl`: the Voros
  coefficient of the Weber equation and its identification with the Stirling correction
  of `log Γ(½+ν)`)
- **[KoikeTakei11]** T. Koike, Y. Takei. *On the Voros coefficient for the Whittaker
  equation with a large parameter - some progress around Sato's conjecture in exact WKB
  analysis.* Publ. RIMS 47 (2011) 375-395. (the companion computation; the Bernoulli
  closed form our `weber_voros_series` implements)
- **[DLMF-Gamma]** *NIST Digital Library of Mathematical Functions*, §5.11
  (Gamma function asymptotic expansions), dlmf.nist.gov. (the expansion of
  `log Γ(z+½)`, whose `(2^{1−2n} − 1)` factor distinguishes it from the plain Stirling
  series; the functional equation and `Γ(½) = √π` are the `weber_log_gamma` oracles)

## The cluster bridge

- **[IN14]** K. Iwaki, T. Nakanishi. *Exact WKB analysis and cluster algebras.*
  J. Phys. A 47 (2014) 474009. arXiv:1401.7094. (**the** bridge reference: the
  equation normalization `ħ²ψ″ = Q(z)ψ` with `Q = V − E` that this package adopts
  package-wide, the Stokes-graph → ideal-triangulation dictionary, Voros jumps as
  y-mutations)
- **[IN16]** K. Iwaki, T. Nakanishi. *Exact WKB analysis and cluster algebras II:
  simple poles, orbifold points, and jump trajectories.* J. Phys. A 49 (2016) 175201.
  arXiv:1409.4641. (the extended dictionary: poles of `Q` give punctured/orbifold
  surfaces rather than disks - the classification `src/rational_potentials.jl`
  implements for irregular poles and refuses for orders 1 and 2)
- **[FST08]** S. Fomin, M. Shapiro, D. Thurston. *Cluster algebras and triangulated
  surfaces I: cluster complexes.* Acta Math. 201 (2008) 83-146. arXiv:math/0608367.
  (ideal/tagged triangulations, flips, the signed adjacency matrix - the target of
  the `triangulation.jl` construction)

## BPS spectra and wall-crossing

- **[KS08]** M. Kontsevich, Y. Soibelman. *Stability structures, motivic
  Donaldson-Thomas invariants and cluster transformations.* arXiv:0811.2435.
  (the KS wall-crossing formula as an identity of quantum-dilogarithm / cluster
  transformations; `verify_ks_*` in `sw_bps.jl` closes it order by order)
- **[DG09]** T. Dimofte, S. Gukov. *Refined, motivic, and quantum.* Lett. Math. Phys.
  91 (2010) 1-27. arXiv:0904.1420. (the refined/motivic wall-crossing formula and the
  refined index `Ω(γ, y)`: a hypermultiplet contributes `1`, a vector multiplet
  `−y − y^{-1}` - the reading `refined.jl` and the SU(2) refined identity implement,
  with the quadratic-refinement sign measured rather than imported)
- **[Kel11]** B. Keller. *On cluster theory and quantum dilogarithm identities.* In
  Representations of Algebras and Related Topics, EMS (2011) 85-116. arXiv:1102.4148.
  (the quantum-dilogarithm identity attached to a maximal green sequence, and its
  independence of the sequence - the cluster-side statement the refined bridge lands on)
- **[Rei10]** M. Reineke. *Poisson automorphisms and quiver moduli.* J. Inst. Math.
  Jussieu 9 (2010) 653-667. arXiv:0804.3214. (wall-crossing invariance of the ordered
  dilogarithm product - the identity behind the quantum pentagon reached from two
  potentials in `verify_refined_wall_crossing`, and behind the Kronecker case)
- **[GMN10]** D. Gaiotto, G. Moore, A. Neitzke. *Four-dimensional wall-crossing via
  three-dimensional field theory.* Comm. Math. Phys. 299 (2010) 163-224.
  arXiv:0807.4723. (the physics derivation of the KS formula; the charge lattice and
  its Dirac pairing; **and** the finite-radius twistor construction implemented in
  `hyperkahler.jl` - the Darboux coordinates `𝒳_γ(ζ)`, the Laurent form of `ϖ(ζ)`,
  the hyperkähler metric on the moduli space of the theory on ℝ³ × S¹)
- **[GMN-TBA]** D. Gaiotto, G. Moore, A. Neitzke. *Wall-crossing, Hitchin systems,
  and the WKB approximation.* Adv. Math. 234 (2013) 239-403. arXiv:0907.3987.
  (BPS states from finite WKB trajectories, and the TBA integral equations solved in
  `tba.jl`: exact quantum periods from BPS data alone, chamber independence; §6 also
  carries the quadratic-differential ↔ marked-surface dictionary that
  `rational_potentials.jl` and the general `ideal_triangulation` implement - a pole of
  order `p ≥ 3` is a boundary circle with `p − 2` marked points)
- **[OV96]** H. Ooguri, C. Vafa. *Summing up D-instantons.* Phys. Rev. Lett. 77
  (1996) 3296-3298. arXiv:hep-th/9608079. (the one-hypermultiplet hyperkähler
  geometry whose instanton sum `hyperkahler.jl` reproduces in closed form - the
  Bessel-K series `ooguri_vafa_instantons`)
- **[GMN-SN]** D. Gaiotto, G. Moore, A. Neitzke. *Spectral networks.* Ann. Henri
  Poincaré 14 (2013) 1643-1731. arXiv:1204.4824. (spectral networks; the
  higher-rank generalization of the Stokes graphs used here)
- **[BDP14]** T. Brüstle, G. Dupont, M. Pérotin. *On maximal green sequences.*
  Int. Math. Res. Not. 2014 (2014) 4547-4586. arXiv:1205.2050. (maximal green
  sequences = the BPS spectrum ordering; the enumeration lives in the cluster core)

## Seiberg-Witten SU(2)

- **[SW94]** N. Seiberg, E. Witten. *Electric-magnetic duality, monopole
  condensation, and confinement in N = 2 supersymmetric Yang-Mills theory.*
  Nucl. Phys. B 426 (1994) 19-52. arXiv:hep-th/9407087. (the SU(2) curve, the
  periods `(a, a_D)`, the singularities at `u = ±2Λ²`, the monodromy group)
- **[NS09]** N. Nekrasov, S. Shatashvili. *Quantization of integrable systems and
  four dimensional gauge theories.* arXiv:0908.4052. (the NS limit; why the pure
  SU(2) quantum periods are the Mathieu problem's exact-WKB periods)
- **[BilalFerrari]** A. Bilal, F. Ferrari. *Curves of marginal stability, and weak
  and strong-coupling BPS spectra in N = 2 supersymmetric QCD.* Nucl. Phys. B 480
  (1996) 589-622. (the marginal-stability wall; its crossing of the imaginary axis
  at `0.86·2Λ²` is our `ms_wall` oracle)
- **[DLMF-Mathieu]** *NIST Digital Library of Mathematical Functions*, ch. 28
  (Mathieu functions), dlmf.nist.gov. (the `a_ν(q)` characteristic-value expansions
  used as the rung-1 period oracle; §28.2 for the Floquet exponent `ν` and the fact
  that both `a_r(q)` and `b_r(q)` sit at integer `ν = r` - the band-edge pair whose
  midpoint the perturbative NS series predicts, our rung-5 ħ⁴ oracle)
- **[Carlson95]** B. C. Carlson. *Numerical computation of real or complex elliptic
  integrals.* Numer. Algorithms 10 (1995) 13-26. (the symmetric forms `R_F`, `R_D`
  implemented in-file in `sw_curve.jl` so `K`, `E` are available at complex
  parameter; parameter convention `K(m) = ∫₀^{π/2} dθ/√(1 − m sin²θ)`)

## Double well and exact quantization

- **[JZJ04]** U. Jentschura, J. Zinn-Justin. *Multi-instantons and exact results I:
  conjectures, WKB expansions, and instanton interactions.* Ann. Phys. 313 (2004)
  197-267. (the ZJJ normalization of the two characteristic functions, eqs.
  2.23-2.25 - the convention `double_well.jl` adopts; the exact `g = 0.08`
  eigenvalues quoted in fig. 2.1 are an oracle)
- **[DU14]** G. Dunne, M. Ünsal. *Uniform WKB, multi-instantons, and resurgent
  trans-series.* Phys. Rev. D 89 (2014) 105009. arXiv:1401.5202. (the P/NP relation
  between the two quantum periods, verified in `double_well.jl`)
- **[GT15]** I. Gahramanov, K. Tezgin. arXiv:1512.08466, eq. 4. (the explicit form of
  the Dunne-Ünsal relation transcribed into ledger item 2 of `double_well.jl`)
