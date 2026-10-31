# FEBio v4.13 collapsibleChannel: static-wall and fixed-shape fluid evidence, applied-load discrepancy and coupled comparison not established

**Status: DRAFT diagnostic evidence package. The coupled FEBio-vs-oomph-lib comparison is NOT ESTABLISHED (no full 140-step run was made). Nothing here is a successful verification.** Stock, unmodified FEBio v4.13 (tag v4.13, clean; binary sha256 `f5ae1ed6f347f875c392f93d847f2369f437345a9b785296d46b0bd6fbb72f02`), solids4foam `development` commit `b200718f4f8390c90f9c7e3ded8b74f3e13941cc` (case inputs and reference files), 6 Slurm jobs, 0.195 allocated CPU-h. No stock FEBio source change is proposed and no manuscript change is requested.

Evidence labels: PROVEN / SUPPORTED / UNRESOLVED / FAILED / NOT TESTED. Classes: VERIFIED AGREEMENT / QUALIFIED AGREEMENT / DIAGNOSTIC DISCREPANCY / NOT ESTABLISHED.

## Executive summary
- **Static wall (CC-S1, solid only).** Stock hex8 at the tutorial's wall mesh (80x8, aspect 20) is 11.6 % too stiff against the oomph-lib beam value -0.143628 m: **DIAGNOSTIC DISCREPANCY** (PROVEN). Errors 80x8 -11.63 %, 160x4 -3.054 % (missed the pre-registered 3 % gate by 0.054 pp), 160x8 -3.000 %, 320x8 -0.757 %; second-order convergence (ratios 3.9, 4.0; Richardson limit within 0.01 % of the beam value). 320x8: **QUALIFIED AGREEMENT** (not verified). Amendment A1 (switch to the 320x8 wall and a 7,040-cell fluid mesh) is reviewed in `docs/PLAN_AMENDMENTS_AND_ACCEPTANCE.md`: acceptable as a documented, conservative deviation, with stated qualifications.
- **Fixed collapsed-shape fluid (CC-F1).** Mass balance 6e-7, steady drift 1e-11, streamwise parity 7e-5: pass their own gates (**QUALIFIED AGREEMENT**, one mesh, no external reference). The pre-registered force gate **failed**: the nodal wall pressure that FEBio applies to the wall differs from the element-mean pressure by -6.99 % (F_y) and -8.63 % (first-mode load) through a **wall-normal row-parity mode** (wall row = row 2 to 4e-3 Pa; row 1 differs by about 1 Pa).
- **Dynamic diagnostic (D1, 20 coupled steps through the first collapse).** The hypothesis that this is negligible in the dynamic regime is **refuted**: applied nodal wall load exceeds the element mean by +19 % to +65 % during loading (193 vs 77 Pa on the wall row vs row 1 at x = 10, t = 0.25 s). A free-body check reproduces the applied nodal load to 0.1-1 % for t = 0.225-0.35 s (valid only to |u| <= 0.1425 m). The wall lags oomph-lib by up to 12.8 % of the peak (t = 0.425 s). The lag matches the analytical load-balance law in sign and size; this is **not proof of a sole cause** (single mesh and time step, rhoi 0.9, E_mesh 1, wall 0.76 % too stiff).
- **Not established:** first trough value and time (reference -0.2160 m at 0.4685 s), trough-to-trough period (0.65 s; the reference README's "0.8 s" is not reproduced from the CSV), rebound, second trough, and the coupled comparison. A meaningful coupled comparison needs a fluid operator without wall-normal pressure parity, which stock v4.13 does not provide.

## Compact pass/fail matrix
| Quantity | Reference | FEBio result | Pre-registered gate | Classification | Evidence |
|---|---|---|---|---|---|
| Mesh counts: tutorial 3,200 fluid / 640 solid; pilot 7,040 fluid + 2,560 wall | tutorial inputs | counted from regenerated inputs | n/a | consistent | PROVEN |
| Static wall, tutorial mesh 80x8 | -0.143628 m | -11.63 % | 3 % | **DIAGNOSTIC DISCREPANCY** | PROVEN |
| Static wall 160x4 | same | -3.054 % | 3 % | gate missed by 0.054 pp (A1) | PROVEN |
| Static wall 320x8 | same | -0.757 % (2nd order, 3 meshes) | 1 % | **QUALIFIED AGREEMENT** | SUPPORTED |
| Static reference re-derivation (own beam code) | -0.143628 m | -0.14359 m extrapolated (-0.03 %) | n/a | same model class, about 0.05 % | SUPPORTED |
| Fixed-shape fluid: mass, steadiness, streamwise parity | own gates | 6.0e-7; 1.2e-11; 6.8e-5 of V | 1e-3; 1e-3; 5 % | **QUALIFIED AGREEMENT** (one mesh) | PROVEN (data) |
| Fixed-shape applied wall load F1 vs F2 | n/a | F_y -6.99 %, first-mode -8.63 % | 3 % (stop above 8 %) | **DIAGNOSTIC DISCREPANCY** | PROVEN (data); cause SUPPORTED |
| Dynamic applied wall load (D1, t 0.10-0.30 s) | n/a | +19 to +65 % | H_D: <= 3 % | **DIAGNOSTIC DISCREPANCY** (H_D refuted) | PROVEN (data) |
| Free-body check of the applied load | n/a | equals nodal load 0.1-1 % (t 0.225-0.35 s) | n/a | supports applied = nodal | SUPPORTED |
| Wall midpoint history vs oomph-lib (20 steps) | converged oomph-lib | max -12.8 % of peak at t = 0.425 s; trough not reached | informational | diagnostic observation only | not a validation |
| First trough value / time | -0.2160 m / 0.4685 s | not reached at 0.5 s | within 8 % | **NOT ESTABLISHED** | NOT TESTED |
| Trough-to-trough period, rebound, second trough | 0.65 s; -0.0404 m; -0.2037 m | no full run | n/a | **NOT ESTABLISHED** | NOT TESTED |
| Coupled FEBio vs oomph-lib | n/a | no full run (stop rule) | G1-G3 | **NOT ESTABLISHED** | NOT TESTED |

## Contents
- `docs/` SCIENTIFIC_SUMMARY; EQUIVALENCE_AUDIT; PLAN_AMENDMENTS_AND_ACCEPTANCE (pre-registration, A1 and A2 with my review and the acceptance criteria); RESULTS (CC-S1, CC-F1, D1); PARITY_DEFINITIONS; LOAD_BALANCE_ANALYSIS; SLURM_ACCOUNTING; LIMITATIONS_AND_UNRESOLVED; RECOMMENDATION_AND_REOPENING (incl. relation to PRs #27/#28); `archive/PILOT_RESULTS_original.md` (the original report, superseded where `provenance/INCONSISTENCIES.md` differs).
- `results/` CSV/JSON tables (S1, F1, D1, free-body, oomph-lib reference features, independent re-verification of 141 quantities, beam re-derivation); `figures/` five figures regenerated from the tables and inspected; `inputs/` the four S1 inputs and the frozen-shape file (the large F1 and D1 inputs, 1.7-2.0 MB, are regenerated and hash-checked by `reproduce.sh`); `raw/` compact raw outputs (final-step and per-step node logs, gzip; logs) and the original analysis results; `reference/` the oomph-lib CSVs with `ORIGIN.md` (origin and hashes); `scripts/` generator, analysis, reference and verification code; `provenance/` hashes, Slurm accounting, ledger extract, retained-HPC-output list (relative path, size, sha256), `REPRODUCTION_RECORD.md`, `INCONSISTENCIES.md`; `reproduce.sh`; `SHA256SUMS`.
- Not shipped: `.xplt` files, F1/D1 `.feb` inputs and `mesh2d.json` (regenerable), empty Slurm output files; their sizes and sha256 are in `provenance/retained_hpc_outputs.csv`.

## Reproduce (no simulation)
`./reproduce.sh` (Python 3 with numpy, scipy, matplotlib; about 30 s on one core; `FULL=1` adds the 40-element beam run): regenerates the six inputs and checks their hashes against the shipped inputs and the ledger, rebuilds the shape file, re-runs the three original analysis scripts and the independent re-computation on the stored raw data, rebuilds the tables and figures. Result: `provenance/REPRODUCTION_RECORD.md`. Re-running FEBio itself needs v4.13 with MKL; commands are in `provenance/REPRODUCTION_RECORD.md`.

## What this package does not claim
No successful verification of collapsibleChannel by FEBio; no statement about solids4foam's results; no statement that F1 or F2 is the physically correct load; no claim that parity is the sole cause of the wall lag; the 20-step diagnostic is never presented as agreement.
