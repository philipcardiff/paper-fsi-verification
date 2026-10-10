# Manuscript update, October 2026

Branch `manuscript/evidence-update-2026-10`, based on
`origin/womersley-temporal-diagnosis` (PR #1, `0b08caa`). This file lists every
changed claim (old → new) with its evidence source, then the TODOs that
remain. The evidence branches were read but not merged, so this branch has no
evidence files. The sources are:

- **3DT**: `paper-data/3dtube_level3_findings.md` on `evidence/3dtube-level3`
  (PR #3). Section 12.6 and later are used; sections 1–11 are superseded.
- **HT**: `paper-data/hronturek_fsi3_4x_findings.md` on
  `evidence/hronturek-fsi3-4x` (PR #4). Sections 1–16 are used; the claim in
  section 12 that "the solid explains" the displacement changes is retracted
  in section 13.
- **HT-EB**: section 16 of the same file on `evidence/hronturek-energy-balance`
  (PR #546 exposure, energy balance, `u_x` ~ A²).
- **HT-RV**: solids4foam `verification/hronturek-replay-variants` `d700b64c8`,
  `HronTurek/verification/README.md`, section "Fluid-side variants". This
  includes the coordinator correction that withdraws `gcl`.
- **HV**: solids4foam `verification/hessenthaler-phaseI` `9f33c42a2`,
  `hessenthalerFsi/verification/README.md`, and solids4foam issue #552.
- **BEAM**: `paper-data/beam_crossflow_graded_mesh_findings.md` §9 on
  `evidence/beam-crossflow-graded-mesh` (PR #2), plus
  `status_review_2026-10-09.md` §6.
- **PR546**: solids4foam PR #546 description, including the audit section.
- **SRC**: solids4foam `development` source and tutorial dictionaries.

## Front matter and methods (commit 1)

| Location | Old | New | Source |
|---|---|---|---|
| article.tex title | "A quasi-monolithic cell-centred finite volume approach for fluid-solid interaction" | "Benchmark problems are not benchmark solutions: an open verification suite for continuum fluid–solid interaction with quantified numerical uncertainty" | brief |
| abstract | Quasi-monolithic JFNK method; "second-order spatial and temporal accuracy, agreement with partitioned reference solutions" | Verification thesis, A/B/C/C-weak categories, the partitioned framework, and the FSI3, 3dTube and convergence-test findings | all |
| keywords | Monolithic, quasi-monolithic, JFNK | Code verification, benchmark, reference solution, numerical uncertainty, partitioned coupling | – |
| article.tex top `\todo` | "title… not rewritten" | Notes what was rewritten; Section 2 still needs a consistency check | – |
| Data availability | `feature-petsc-snes` branch; solid-benchmarks repo | solids4foam tutorials and `Allverify` drivers; `\todo` for tag/DOI | SRC |
| Intro, opening paragraph | "This paper presents a quasi-monolithic … method" | Benchmark problems vs benchmark solutions | brief |
| Intro, monolithic FE/FV gap paragraphs (old lines 28–51) | Motivation for the quasi-monolithic solver | Removed. Replaced by a paragraph on the verification literature (Turek2006, Richter2012/2017, Heil, Eken, Lozovskiy, Tuković) | – |
| Intro, contributions 1–4 | Quasi-monolithic formulation, JFNK preconditioner, PETSc implementation, benchmark study | Reference classification; spatial/temporal/iterative studies; defects and misleading agreement found by verification; executable open suite | brief |
| Intro, Robin citation | `\citep{Tukovic2019}` (a Simo–Reissner beam paper) | `\citep{Tukovic2018b}` (Robin pressure condition) | bib |
| sectionNumericalModel.tex | Quasi-monolithic method, Algorithm 1 (quasi-monolithic) | Partitioned framework: pimpleFluid ALE fluid (PIMPLE, Rhie–Chow, BDF2, BDF2-consistent mesh flux), velocityLaplacian mesh motion, small-strain and total-Lagrangian FV solids (StVK, neo-Hookean, linear elastic; standard and high-order MLS variants; PETSc SNES), directMap/AMI transfer, Dirichlet–Neumann (fixed, Aitken, IQN-ILS, predictor), Robin–Neumann (`elasticWallPressure`, unrelaxed), convergence test r1/r2 with legacy min and strict max, Algorithm 1 (partitioned time step), toolchain | SRC |
| sectionNumericalModelQuasiMonolithic.tex | – | New file, **not input**: the old section, with labels prefixed `qm:` | – |

## 3dTube (commit 2)

| Location | Old | New | Source |
|---|---|---|---|
| BS 3dTube "Reference" | "independent monolithic finite element codes, both implicit Euler" | Lozovskiy is FE. Eken has a side-centred FV fluid with an FE solid. Both use implicit Euler at Δt 1e-4 with no time-step study | 3DT §3 |
| BS "Evidence" | "Two mesh levels … a third has not been run" | Three levels (16 000 / 128 000 / 1 024 000 fluid cells; Δt 2.5e-5 / 1.25e-5 / 6.25e-6). Mean iterations 4.8 / 3.9 / 3.2 | 3DT §12.6, §3 |
| `tab:3dtube` | Two levels: u_r 0.15984/0.15976; u_z −0.08823/−0.08660; t_arr 5.874/5.854; c_p 4.771/4.695 | Three corrected levels: u_r 0.15978/0.15946/0.15851 (−0.20%, −0.60%, no order); u_z −0.08819/−0.08659/−0.08571 (order 0.87); t_arr 5.873/5.843/5.825 (0.72); c_p 4.654/4.635/4.594 (no order) | 3DT §12.6 |
| Space bullet | "u_r,max changes by only 0.05%" | u_z and t_arr monotone, orders 0.87/0.72 (fixed Δt 0.98/0.65); u_r,max changes grow, uncertainty ≥0.6%, spatial convergence not established | 3DT §12.6–12.7 |
| Time bullet | L1 halvings 1.74/0.83/0.016% | Fine levels 0.1–0.2% (L2 +0.16/+0.18/+0.10%; L3 +0.10/+0.20/+0.17%); the path is mesh-dominated | 3DT §12.6 |
| Coupling | RN/IQN-ILS ≤0.08% (history 0.27%); 4.7 vs 15.5 iterations | 0.07% (history 0.26%) | 3DT §12.6 |
| New bullet | – | Defect (a), −O3 tmp-reuse aliasing (axial wall shear up to ≈2.5× on step 1; root cause in OpenFOAM tmp reuse with `__restrict__`; the audit fixed ~15 more sites). Defect (b), scale-dependent barycentric degeneracy (1/3, 1/3, 1/3) on L2–L3 only. Both fixed in PR #546. Platform agreement <0.01% | 3DT §12.5–12.6; PR546 (audit: 2 + 13 sites) |
| Literature bullet | "FE/FV difference … +2.5% / −1.4%" | Implicit Euler damping 0.036–0.038 mm = 92–115% of the gap (86–126%); within +2.5/−1.4% (L1), +1.3/−2.6% (L2); Euler still 4.4% low at 1.25e-5. "Much of the discrepancy is reproduced by the difference in time integration"; not FE vs FV; L3 is 1.0% above Tuković | 3DT §6, §12.6 |
| BDF2 comparison bullet | "L2 1.8% above Tuković; c_p 2.4% below estimate" | Removed. Replaced by "L3 1.0% above Tuković" | 3DT §12.6 |
| Lesson | "u_r,max converged, u_z,min not" | Case acts as a defect detector; u_z and t_arr sub-nominal; u_r,max not converging; C-weak | 3DT §12.7 |
| `\todo` | ESSENTIAL: run level 3 | Single-platform L3, L4 or independent solution, upstream report, small-Δt instability | 3DT §12.8 |
| Tables (overview, orders, QoI, published, cost, tiers) and the CB text on the time/space balance | Two levels; "u_r,max 0.05%"; FE "21–24%"; cost "L1 425 s, L2 7 197 s on 4, L3 >10 h est."; "an order of magnitude below" | Three levels; orders as above; "20–25%, 92–115% reproduced"; cost L1 504 s/4, L2 3 641 s/8, **L3 22 502 s (6.25 h) on 32 ranks (720 core-h)**; "several times smaller (0.1–0.2% vs 0.3–1.0%)" | 3DT §2, §12.6 |

## Hron–Turek (commit 3)

| Location | Old | New | Source |
|---|---|---|---|
| BS Reference | FSI2/FSI3 "uncertain by ≈1%" | FSI2 L3→L4 0.9–2.5%. FSI3 L3→L4 3.8% (u_x mean), 4.0% (u_x amp), 4.5% (drag amp), 2.6% (lift), 1.6% (u_y): 3–5% uncertainty for u_x and force amplitudes, ~2% for u_y, ≲0.3% for drag mean and frequency | HT §3 |
| BS Reference | Summary values (no provenance) | Level 4 of the 2006 proceedings at Δt 0.0005, kept on Featflow only as old values | HT §3 |
| Driver check | 0.7% (1.1% drag amp) | 0.6% (0.8% drag amp; frequencies 0.3%) | HT §3 |
| Evidence | "two mesh levels; the 4x level has not been run" | FSI3 three levels (5 336 / 21 344 / 85 376 fluid; 630 / 2 520 / 10 080 solid); every step converged; 7.0/7.0/6.7 iterations | HT §2 |
| `tab:hronturek` | FSI3 + FSI2 errors at 1x/2x; FSI3 apparent orders 1.1–3.8 | FSI3 1x/2x/4x values, Featflow, error at 4x, observed order (u_x −2.232/−2.788/−3.085, 0.91; u_y 29.85/34.18/36.23, 1.07; …); smoothed force amplitudes | HT §4 |
| `tab:hronturek_fsi2` | (inside `tab:hronturek`) | New table with the unchanged FSI2 1x/2x errors | inventory |
| FSI2/FSI3 bullet | "within 1.2–3.1% … apparent orders 1.1 to 3.8 … 1x outside asymptotic range" | FSI3 displacements of order ≈0.9–1.1; the 2x agreement was a crossing (4x u_y +3.5%, u_x −7.1%); 2.6–3.8 withdrawn as artefacts; u_x mean ∝ A^1.88–1.93 so not independent | HT §1, §4; HT-EB "u_x mean" |
| Lift bullet | FSI3 "48.9% → 13.3%" | 4x raw +15.0% (lift) and +18.1% (drag) are noise-inflated; smoothed lift 11.9% → 7.6%, drag amplitude still rising | HT §4, §7 |
| New bullet | – | 2x Δt/tolerance ≤1.7% of the displacement; IQN-ILS stalls at 5–9e-6 with tolerance 1e-6 at Δt 0.00025 | HT §5, §7 |
| New paragraph "Diagnosis" | – | CSM2/CSM3 orders 1.4–1.6 → 1.8–1.9; solid-only refinement 5% of the u_y change (+0.33%), lift −12.3%, drag −5.5%; CFD3 4x within 0.16% (drag mean) and 0.05% (lift) of Featflow, drag amp +2.0%, frequency +1.2%; PR #546 output identical to every written digit; replay Q_in 50.3/72.2/81.0 (order 1.3), 1.1–1.3 for every valid three-level variant, spread 39/21/12 N/m; W change below resolution; one-mode balance 50–98% (hypothesis); geometric-singularity **hypothesis** | HT §12–16; HT-EB; HT-RV |
| `\todo` ESSENTIAL run FSI3 4x | – | Replaced by follow-ups: singularity test, 8x, FSI1 4x, FSI2 4x, windows; `gcl` noted as withdrawn and not used | HT §10; HT-RV |
| Lesson | "references uncertain by ≈1% … FSI2/FSI3 not demonstrably asymptotic" | 3–5% reference uncertainty; the two-level crossing; first-order sequence not caused by the solid, the static fluid or a defect; "disagrees" vs "unresolved"; category C | HT; brief |
| Tables: overview, orders, QoI, published, cost, tiers; RS extended-suite text | FSI3 2 levels; "1.1–3.8"; "2x displacements 2.1–3.1%, lift 13.3%"; "~1%"; FSI2/3 cost "2x 4.0–4.1 h on 7–8 ranks"; "needs FSI3 4x"; "≈1% uncertainty" | 3 levels; 0.9–1.1; 4x statement; 3–5%; **FSI3 1x 1.4 h serial, 2x 5.5 h on 8 ranks, 4x 40.8 h on 32 ranks** (the FSI2 row is unchanged); three levels, order ≈1; 3–5% | HT §2–4 |
| CB "two levels" paragraph | "FSI1–3 and 3dTube have two levels; FSI3 apparent orders 1.1–3.8" | Third levels changed both conclusions (crossing; changing interpolation operator). Super-nominal orders are listed as pre-asymptotic | HT, 3DT |

## Methodology and iterative error (commit 4)

| Location | Old | New | Source |
|---|---|---|---|
| VM `\todo` on the quasi-monolithic Section 3 | "Section 3 currently describes the quasi-monolithic solver …" | A sentence stating that every result is partitioned and that validation is separate | – |
| VM Table 1 | C: "beamInCrossFlow, original form"; C-weak: "beamInCrossFlow modified form" | C: "beamInCrossFlow (Richter 2012; in progress)"; the modified form is removed | BEAM |
| VM Iterative paragraph | Tolerance sweeps only | Adds: what the convergence test measures, that a steadiness check is blind to accumulated error, and the strict test plus no restarts from weak-test states, given as a **recommendation**; states that the completed studies used the legacy test; `\todo` to audit r1 at acceptance | HV; SRC |
| VM new paragraph "Code, build and platform" | – | The toolchain is part of the verification configuration; cross-platform and symmetry/decomposition checks | 3DT §12 |
| CB intro | "the one such failure found, the womersleyTube temporal order …" | Adds sub-nominal convergence and points to the new defects subsection | – |
| CB new observation "Coupled limit-cycle QoIs …" | – | FSI3 components behave well but the coupled order is ≈1; "disagrees" vs "unresolved" | HT, HT-RV |
| CB iterative errors | "Only FSI2 contains a coupling-tolerance sweep" | Adds the FSI3 force noise (7.3 vs 1.7 N/m rms) and the IQN-ILS stall | HT §7 |
| CB new paragraph "A converged-looking QoI …" | – | Hessenthaler legacy test: ~3 iterations per step; after 1 s every step has r1 > 1e-4 (median ~6.8e-2; 87–91% > 1e-2; r2 median 3.5e-5); lag ~0.03 mm; creep +0.003–0.005 mm/s (16.84/16.90/16.98 mm at 30/60/80 s); strict test 10–11 iterations, plateau 17.0195 ± 0.0002 mm over 7.5 s; restarts retain history (17.02 vs 16.90 mm); issue #552 | HV; #552 |
| CB new subsection "Defects exposed by verification" (`tab:defects`) | – | Womersley end flux, 3dTube aliasing, 3dTube barycentric weights, Hessenthaler convergence test; toolchain in the verification matrix | PR #1, 3DT, HV |
| BS beamInCrossFlow | "[PROVISIONAL]"; Richter **2015** values 5.95e-5 m and 1.33 N; uniform refinement +0.07%/−1.4%; modified form vs Tuković | "[IN PROGRESS]"; Richter **2012** (bib entry already present), extrapolated 5.924e-5 m and 1.327 N (≤1%); the quoted values are rounded; earlier runs solved a different problem (beam position, point A, inflow, start-up), so the comparison is withdrawn; no new results | BEAM §9; status §6 |
| BS Hessenthaler | "[IN DEVELOPMENT] … no case exists" | "[IN PROGRESS]"; primarily validation; separate verification campaign with a fixed model; status is the convergence-test finding only; `\todo` on the PR #546 A/B check and μ* | HV; status §3, §7 |
| Discussion `\todo` | – | Adds three themes: toolchain, limit-cycle QoIs, accumulated iterative error | – |

## Cleanup (commit 5)

| Location | Old | New | Source |
|---|---|---|---|
| BS blobInTreacle | "supports at least second-order behaviour"; "Richardson estimate 0.12026 m, finest within 0.1%" | Reported as observed self-convergence behaviour only; super-nominal; no extrapolated value | status §8 |
| BS mokLidDrivenCavity | "suggesting a limit near 0.285 m"; lesson "≈1% in the peak"; todo "Richardson estimates" | No limit extrapolated; "last change 0.78% in the peak"; no Richardson estimate is justified | status §8 |
| BS cavityFlexibleBottom | "converges, at approximately second order" | One triple (net 1.96), level 4 diverged: second order not established | status §8 |
| CB orders table | cavity "1.96 (net)"; blob "2.40" | "(one triple, level 4 diverged)"; "(super-nominal)" | – |
| BS overview, CB orders/QoI/published/cost, RS tiers | beamInCrossFlow "4/5 levels, not asymptotic, 8x 11 h on 64 cores"; Hessenthaler "periodic, in development" | "in progress" | status §6–7 |
| CB QoI type (iii) | Near-wall force resolution in beamInCrossFlow (hypothesis) | The in-phase fluid force on the Hron–Turek flag | HT-RV |
| CB cost/dimensionality text | 3dTube two levels; beam near-wall limit | 3dTube three levels; FSI3 4x cost; collapsibleChannel solid vs FSI3 solid share | 3DT, HT |
| RS core/extended text and infrastructure list | beam "near-body forces"; no convergence-test or toolchain reporting | FSI3 component studies; report the convergence test, the residual reduction at acceptance, and the compiler/flags | – |
| VM:115 version `\todo` | "3dTube and Hron–Turek on v2412" | 3dTube L1–2 v2512 (confirmed on v2412), L3 v2412 only; FSI3 v2512; Hessenthaler v2412; which results predate PR #546 still has to be stated | 3DT §12.6; HT §2; HT-EB |

## Remaining TODOs (in the .tex unless noted)

1. Check that Section 2 (math model) is consistent with the partitioned
   method, e.g. velocity- vs displacement-based solid unknowns.
2. Give the solids4foam commit, OpenFOAM release, compiler and flags, PETSc
   and platform for every result (VM `\todo`), and say which results predate
   PR #546. Check that the `predictor` description matches each commit
   (changed in PR #492).
3. Audit r1 at acceptance for every completed study run under the legacy test
   (issue #552 list: FSI1, cavityFlexibleBottom, beamInCrossFlow).
4. 3dTube: a second-platform L3 replicate, an L4 or an independent BDF2-type
   solution for u_r,max; an upstream OpenCFD report on tmp reuse with
   `__restrict__`; the small-Δt backward instability.
5. Hron–Turek: test the geometric-singularity hypothesis; 8x FSI3; FSI1 4x;
   FSI2 4x; longer windows; cite the Featflow FSI/CSM/CFD tables precisely.
6. beamInCrossFlow: finish the Richter-definition campaign. Settle the inflow
   (peak 0.3 vs 0.2 m/s) and point A from the text of Richter (2012). Run a
   PR #546 A/B check. Verify the Richter (2012) Table 5/6 values quoted from
   BEAM §9.
7. Hessenthaler: strict-test verification runs from t = 0; the PR #546 A/B
   check; record the justification of μ*; add the Hessenthaler et al. (2017)
   bib entry (still `\citetodo`).
8. Essential tolerance sweeps (strict test) for ringAddedMass and
   collapsibleChannel.
9. Discussion and Conclusions are not drafted. Pre-existing `\citetodo`s
   remain (Roache, Oberkampf & Roy, Womersley, Formaggia, etc.).
10. Release tag/commit and Zenodo DOI (Data availability).

## Evidence conflicts to resolve

- **Replay-variant order range.** The brief gives "~1.1–1.5". The `1.5` is
  the withdrawn `gcl` variant. The valid variants run on three levels give
  baseline 1.31, pwall 1.15 and gradS4f 1.10, and the coordinator correction
  says "about 1.1–1.3". The manuscript uses **1.1–1.3**. The wall-velocity
  and mesh-diffusivity variants ran on two levels only (mesh_inv was a
  screening, not valid for restarts), so no order exists for them.
- **"Accepted after ~1 iteration".** The README and issue #552 give a mean
  of ~3.1 iterations per step after 1 s; r2 can fall below the tolerance
  after a single relaxed iteration. The manuscript says "after about three
  iterations".
- **3dTube "now v2512".** Corrected L1/L2 are v2512 on xenosim, reproduced
  exactly on v2412/MeluXina. Corrected L3 ran on **v2412 (MeluXina, −O2,
  128 ranks)** only. The 6.25 h on 32 ranks figure is the xenosim L3 run
  from before the fixes. The cost is representative but is not the corrected
  run.
- **Hron–Turek 2x cost.** 5.5 h ran on **8** ranks, not 32. Only 4x
  (40.8 h) used 32 ranks.
- **Featflow FSI3 L3→L4.** The drag amplitude changes by 4.5%, not "~4%",
  and u_x by 3.8% (mean) and 4.0% (amplitude). The manuscript quotes the
  exact values.
- **Hessenthaler and beam campaigns** ran on trees without PR #546
  (v2412, −O2, AMI for Hessenthaler). The aliasing defect is −O3-only, so
  −O2 builds should be immune to defect (a); defect (b) can still apply.
  The A/B check is outstanding.
