# Verification programme status review — 2026-10-09

A read-only review of the live state of the paper repository, solids4foam
PRs/branches, and XenoSim/MeluXina jobs as of 2026-10-09 ~23:50 CEST.
It is a snapshot to guide the next work. It is not paper text.

## 1. Paper repository

- `main` is at `5938e09` (2026-10-03). `compress-sections-4-7` is already
  merged.
- Open PRs, all with no reviews and all mergeable:

| PR | Branch | Head | Contents |
|---|---|---|---|
| #1 Womersley | `womersley-temporal-diagnosis` | `0b08caa` | The only PR that edits `.tex`. Depends on solids4foam #533 being merged. |
| #2 beamInCrossFlow | `evidence/beam-crossflow-graded-mesh` | `8a2b4fe` | OLD graded-mesh evidence. Contains no Richter-rebuild evidence. |
| #3 3dTube | `evidence/3dtube-level3` | `30ec6cf` | Corrected 3-level evidence. Its sections 1–11 are superseded by section 12.6. |
| #4 Hron–Turek | `evidence/hronturek-fsi3-4x` | `c2c1510` | Sections 1–15: CSM3, solid refinement, CFD3, prescribed-motion ALE. Section 12 is retracted by section 13. |

- Suggested order: merge #1 first, then rewrite the 3dTube and Hron–Turek text
  on top of it. Otherwise there are about 8 hunk conflicts at VM:115, CB:7–111
  and RS:22.

## 2. solids4foam

| PR | State |
|---|---|
| #546 `fix/fsi-tmp-alias-and-barycentric-weights` `28751407c` | CI green on all three regression arms; mergeable. The -O3 reproducer is only on `verification/3dtube-level3`. Other `tmp & tensorField` call sites have not been audited. |
| #533 Womersley `a507bfd41` | Builds and regression arms pass; markdown-lint fails. |
| #544 beam `4c3f8b9d7` | The branch is now the Richter rebuild (beam x = 0.4–0.5, E = 1.4 MPa). The PR title and body are stale. All 3 regression arms FAIL, most likely because the reference values were not regenerated. |
| #535 Hessenthaler `4f549b500` | Stale and conflicting. The real work is on `verification/hessenthaler-phaseI` `683b0861` (`Allverify` plus about 2,000 lines of scripts), which has no PR. |
| #539 channelLeaflet | Conflicting; fails check-whitespace. |

Hron–Turek diagnostics are five stacked branches with no PR:
`fsi3-4x → csm3 → fsi3-solid-refinement → cfd3 → ale-prescribed`
(`27c9bfdf`). The replay scripts on XenoSim are uncommitted.

## 3. Code-defect exposure of running campaigns (HIGH PRIORITY)

Neither MeluXina tree contains #546 (checked with the GitHub compare API):

- beam: `/project/home/p201287/philip/solids4foam-beam-crossflow` at
  `d6679559`, with uncommitted edits;
- Hessenthaler: `/project/home/p201287/philip/s4f-verif-v2412/solids4foam` at
  `aee0c8e3`.

Both are built with OpenFOAM v2412 and foss-2024a (GCC), the MeluXina toolchain
on which the 3dTube cross-platform difference was traced to these defects.
Hessenthaler uses AMI interface transfer. Until an A/B comparison shows the
effect is negligible, treat these campaigns as possibly contaminated. The
owning agents were asked on 2026-10-09 to confirm, rebuild with #546, run an
A/B comparison and restart.

## 4. Case readiness

| Case | Class | Status |
|---|---|---|
| ringAddedMass | A | Ready. Minor open item: IQN amplitude loss (≤2.1%/period, decreasing with mesh). |
| womersleyTube | A | Ready once #533 and paper #1 are merged. |
| collapsibleChannel | B | Ready. Fine levels reach the reference/model-form floor. |
| Hron–Turek FSI3 | C | Diagnosis nearly closed, see §5. No 8x run. |
| 3dTube | C-weak | Ready as a defects plus time-integration case. u_r,max spatial convergence not established. No L4. |
| beamInCrossFlow (Richter) | – | Not ready. See §6. |
| Hessenthaler Phase I | validation; verification in progress | Not ready. Time/steady-state error dominates, see §7. |
| Mok, cavityFlexibleBottom, blobInTreacle | C-weak | Supplementary. Remove the unsupported extrapolations. |
| channelLeaflet | B-provisional | Paused. |

## 5. Hron–Turek FSI3: what is now known

- Rigid CFD3 at 4x is within 0.05–0.16% of Featflow. On the matched path the
  drag mean is 446.01, 441.81, 440.15 N/m and the lift amplitude is 296.0,
  424.4, 438.0 N/m. The fluid operator alone is not the source of the
  first-order behaviour.
- Under prescribed flag motion (ALE), forces change by about 0.1–3% from 2x to
  4x and are non-monotone at the noise level.
- Solid-only refinement at 2x fluid explains about 5% of the u_y-amplitude
  change and none of the u_x-mean change.
- Coupled u_y amplitude: 29.85, 34.18, 36.23 mm (+6% from 2x to 4x).
- Working hypothesis, not yet tested: the limit-cycle amplitude is set by a
  small net energy balance per cycle, so small force-phase errors are amplified
  into large amplitude changes.
  - Test: net work ∮F·v over the replayed trajectory on 1x/2x/4x (XenoSim job
    10728 and the finished replay runs).
  - The u_x mean is roughly quadratic in the u_y amplitude, so it is not an
    independent QoI.
- The force-noise/tolerance-floor study was never run.

## 6. beamInCrossFlow (Richter)

- The committed case geometry and material match Richter.
- The inlet `maxVelocity` 0.3 → 0.2 edit used by L0–L4 is uncommitted. The
  committed README ("peak 0.3, mean 0.2") is self-inconsistent: a
  bi-parabolic profile has mean = 4/9 × peak.
- The inflow must be settled from the text of Richter (2012), not from
  agreement with Richter's values.
- L0–L3 (v2412, dt = 0.00625, t = 8 s):

| Level | u_x(A) | F_x |
|---|---|---|
| L0 | 4.026e-5 m | 1.2154 N |
| L1 | 5.141e-5 m | 1.2942 N |
| L2 | 5.711e-5 m | 1.3223 N |
| L3 | 5.915e-5 m | 1.3311 N |

  - Observed u_x orders are about 1.0 and then 1.5, so the sequence is not
    asymptotic.
  - The L3 agreement with Richter (5.92e-5 m) is not evidence of convergence.
- L4 (32.1M cells, dt/4):
  - Previous attempts diverged.
  - The running chain at about 270 s/step needs about 140 h.
  - writeInterval of about 24 h risks losing most work at each 48 h walltime.
- The approach from below suggests solid over-stiffness (compare the
  Hron–Turek CSM2 result). A solid-only diagnostic should come before L4.

## 7. Hessenthaler

- Fixed μ* = 64.2 kPa; its justification still has to be recorded.
- F2S2 tip y keeps drifting with run length: 16.84, 16.90, 16.98 mm at
  T = 30, 60, 80 s.
- dt = 0.001 vs 0.002 differs by 0.29 mm at T = 60. That exceeds the spatial
  differences among F3S2, F3S3 and F4S2.
- The steady-state/temporal error must be controlled before spatial families
  are meaningful.
- The pseudo-transient runs hv_PTcd and hv_PTe cannot reach their end times,
  and two of them duplicate each other.

## 8. Manuscript (`main`): stale or contradicted claims

- **Monolithic-solver framing.** The title, abstract and contribution 4
  (article.tex:96, 121; Intro:53–60) present the quasi-monolithic solver as
  the source of the evidence. This must be rewritten. The numerical-method
  section describes only that solver; the partitioned procedures behind every
  result are missing.
- **3dTube** (BS:33, 309–346; CB:32–180; RS:29; VM:115):
  - pre-defect two-level numbers;
  - "u_r,max converged, u_z,min not" is wrong in both halves;
  - the FE/FV explanation should be time integration (implicit Euler vs BDF2);
  - neither defect is mentioned.
- **Hron–Turek** (BS:227–290; CB:31–170; RS:28,57):
  - two-level apparent orders of 2.8–3.8 are artefacts;
  - "within 1.2–3.1% of Featflow" and "≈1% reference uncertainty" should be
    3–5% for the amplitudes and u_x;
  - says 4x has not been run, but it has.
- **Womersley:** stale on `main`; fixed in PR #1.
- **Unsupported extrapolations:**
  - blobInTreacle Richardson value 0.12026 m from p = 2.4 (BS:497);
  - Mok limit of 0.285 m (BS:398);
  - cavity "second order" with level 4 diverged (BS:433).
- **Beam:** cites Richter 2015 for values from Richter 2012, and compares
  against a non-matching problem (BS:444–464).
- **Open items:** 24 `\todo`, 15 `\citetodo`; Discussion and Conclusions are
  undrafted.

## 9. Recommended narrative additions

1. **The toolchain is part of the verification matrix.** An -O3 aliasing
   defect produced a cross-platform discrepancy that looked like a platform
   effect.
2. **Coupled limit-cycle QoIs can amplify small operator errors.** If §5 is
   confirmed, a sub-nominal displacement order can be a property of the QoI
   rather than of a single operator.
3. **Classification of the 3-D cases.** The beam and Hessenthaler cases are
   likely to end as "partially verified, with stated uncertainty". The three
   2-D cases remain the backbone.
