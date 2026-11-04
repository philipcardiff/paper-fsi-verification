# Slurm accounting, job ledger, failed and diagnostic attempts, corrections

## 1. Accounting (sacct export, `provenance/slurm_jobs.csv`; ledger extract `provenance/ledger_extract.csv`)
Allocated CPU-h = allocated cores x ElapsedRaw / 3600 (the project convention). 19 jobs, job ids 14165-14183 (consecutive; no other job of the workspace falls in this range or in the tranche), one job at a time (2 cores for stage B, 4 for stage C), partition `main`, stock binary sha256 f5ae1ed6... in all 19 `run_meta_start.json` files. **Total 1864 CPU-s = 0.5178 allocated CPU-h** (budget of the authorised tranche 3.0 CPU-h). sacct ElapsedRaw has one-second resolution, so the very short stage-B jobs carry up to 1 s x 2 cores of rounding each (total effect below 0.01 CPU-h).

| Slurm job | case | role | cores | elapsed (s) | CPU-s | allocated CPU-h | log termination |
|---|---|---|---|---|---|---|---|
| 14165 | B_CSM1_h8_nt4 | diagnostic (generator v1, body-force sign check) | 2 | 1 | 2 | 0.00056 | UNKNOWN |
| 14166 | B_CSM1_h8_nt4 | stage B | 2 | 1 | 2 | 0.00056 | NORMAL |
| 14167 | B_CSM1_h8_nt2 | stage B | 2 | 1 | 2 | 0.00056 | NORMAL |
| 14168 | B_CSM1_h8_nt8 | stage B | 2 | 2 | 4 | 0.00111 | NORMAL |
| 14169 | B_CSM1_h8_nt16 | stage B | 2 | 9 | 18 | 0.005 | NORMAL |
| 14170 | B_CSM2_h8_nt2 | stage B | 2 | 0 | 0 | 0.0 | NORMAL |
| 14171 | B_CSM2_h8_nt4 | stage B | 2 | 1 | 2 | 0.00056 | NORMAL |
| 14172 | B_CSM2_h8_nt8 | stage B | 2 | 2 | 4 | 0.00111 | NORMAL |
| 14173 | B_CSM2_h8_nt16 | stage B | 2 | 5 | 10 | 0.00278 | NORMAL |
| 14174 | B_CSM1_h20_nt2 | stage B | 2 | 2 | 4 | 0.00111 | NORMAL |
| 14175 | B_CSM1_h20_nt4 | stage B | 2 | 8 | 16 | 0.00444 | NORMAL |
| 14176 | B_CSM2_h20_nt2 | stage B | 2 | 2 | 4 | 0.00111 | NORMAL |
| 14177 | B_CSM2_h20_nt4 | stage B | 2 | 5 | 10 | 0.00278 | NORMAL |
| 14178 | B_CSM1_h8_nt16_root025 | stage B | 2 | 9 | 18 | 0.005 | NORMAL |
| 14179 | B_CSM2_h8_nt16_root025 | stage B | 2 | 5 | 10 | 0.00278 | NORMAL |
| 14180 | B_LIN_h8_nt8 | stage B | 2 | 1 | 2 | 0.00056 | NORMAL |
| 14181 | B_LIN_h8_nt16 | stage B | 2 | 4 | 8 | 0.00222 | NORMAL |
| 14182 | C1u_L1_f1.25_K1e6 | stage C diagnostic (unpinned) | 4 | 221 | 884 | 0.24556 | NORMAL |
| 14183 | C1p_L1_f1.25_K1e6 | stage C (pinned) | 4 | 216 | 864 | 0.24 | NORMAL |

The ledger holds 19 rows for this tranche (benchmark `HronTurekFSI1`); the diagnostic job and the later job B_CSM1_h8_nt4 share the case id, distinguished by input hash and run directory.

## 2. Failed and diagnostic attempts (all kept)
1. **Job 14165** (generator v1, `inputs/stageB/DIAG_v1_signtest_B_CSM1_h8_nt4.feb`): normal termination, but u_y = +61.93 mm (upward) and reaction forces not logged properly: the const body load with y = -2 acts opposite to the stated acceleration, and the plot variable `reaction forces` was missing. Diagnostic only; fixed in generator v2 (sign flipped once, variable added); all stage-B inputs regenerated. See `STAGE_B_STRUCTURAL.md`.
2. **Reaction forces zero in every stage-B run** (after the fix as well): explained from source (fixed DOFs have no equation number); the planned equilibrium check was replaced by a stress-based free-body check made from the stored xplt (done after the runs, for this package).
3. **Job 14182, C1u** (unpinned): normal termination, not steady, junction-node pressures 471 / 588 Pa outliers; pre-declared pin rule fired; kept as a diagnostic.
4. **Job 14183, C1p**: normal termination, steady, gates C3-C5 failed; stage C stopped by the pre-registered rule (level 2, K x 10, pin sensitivities and the one-way prediction were therefore never run: NOT TESTED).
5. Analysis-side corrections made before results were interpreted: the control-volume formula in the written pre-registration had the sign reversed (the code initially followed it and gave an absurd lift of 47 N/m; the sign was then corrected in the code; the value of F3 with nodal pressure is nevertheless contaminated by the parity mode, not by the sign).

## 3. Corrections and refinements relative to the archived working reports (`docs/archive/`)
- The archived results report states that the parity mode is "30-80 times the 0.5 Pa cross-flag signal". Precise statement: the uniform load equivalent to the FSI1 deflection is 0.48 Pa (analytical estimate); the parity-cleaned cross-flag pressure difference on the rigid flag is -2.15 Pa; the raw nodal difference is -10.6 Pa; the alternation amplitude is 2.6-4.1 Pa along the flag walls and 15.8 Pa along x = 1.0 m. The ratio to the signal therefore ranges from about 1.2-2 (wall chains) to 7-33 (vertical line), not uniformly 30-80.
- The archived "P2" wall-row-offset values (-81 Pa top, +3.8 Pa bottom) come from a quadratic extrapolation that amplifies an alternating pattern about 7 times; they are not pressure offsets (see `STAGE_C_CFD1.md`). The gate failure is established by the along-wall amplitude a_t.
- The archived reports describe the F2/F3b spread as 4.6 % (of F3b); it is 4.8 % of F2, 4.6 % of F3b; its decomposition is in `STAGE_C_CFD1.md` section 5.
- The archived reports call the flag root "inferred from the stage-B agreement"; the independent geometric justification (Featflow definition page) is in `STAGE_B_STRUCTURAL.md`.
- Stage-B equilibrium, F4 (weak-form force) and F-actual (traction-formula reconstruction) are new in this package.
- The archived pilot plan (stages A-D, budget 6 CPU-h) and roadmap were planning documents; only stages B and C were authorised and run.
