# Pre-registration and acceptance criteria, as run

The pre-registration is `archive/RTR_PILOT_PLAN.md`. Its sha256 was stored before the first submission (`rtr_agent/PLAN_SHA256.txt` of the campaign): **290bf60efc0e9de1153e1b50003424fe1c9bcb70b1efdd80e52fa89b8f8e7834**. Re-verified here: the first 15,916 bytes of the archived file (everything before the section "APPENDED AFTER EXECUTION") hash to exactly that value; the archived file as a whole (with the appended execution notes and the later stall-diagnostic pre-registration) is `c2796334a0bfd1805fd0c53af18512b8c88251696e74870169a0a4a5c0023a90`. The appended text states that it alters no gate. The campaign's own record of generator and reference hashes at plan time is in `../provenance/INCONSISTENCIES.md` (the generator was edited after plan time; the gated runs still regenerate bit-identically).

Fixed settings: stock solver, K_f = 1e7 Pa, rho_inf = 1, dt = 0.01 s (200 steps to 2 s), E_mesh = 1 Pa, vtol = ftol = 1e-6. Cap 3.0 allocated CPU-h; actual 0.1350.

## Gates as pre-registered and outcome
| Stage | Gate (fixed before the runs) | Outcome |
|---|---|---|
| RT-A law1 | all six Cauchy components within 1e-8 of max abs stress at 5 steps | 1.2e-12, PASS (job 14189; the first two designs 14186/14188 had free faces and are not evidence) |
| RT-A rotcheck | node displacement within 1e-10 m at 40 times | 7.0e-13 m, PASS (14187) |
| RT-A twist | E_max(T_A) at M3 at most 1e-3 and order (M2, M3) at least 1.8 | 2.2e-4, order 1.99, PASS (14193 M1, 14195 M2, 14196 M3) |
| smoke | normal termination, 20/20 steps | fluid module PASS; coupled with MR law FAIL (Newton stall at step 5): see `COUPLED_STALL.md` |
| RT-B B1 | velocity E_L2 at most 2 % at M1; ratio M1/M2 at least 3 | 1.747 %, ratio 1.02: ratio FAIL |
| RT-B B2 | E_max(T2) at most 10 % (M1), 5 % (M2), order at least 0.8 | 12.87 %, 7.47 %, order 0.78: FAIL (expected to fail in the plan, predicted 11.7 / 6.0 %) |
| RT-B B3 | pressure parity (theta scatter and radial high-frequency part) at most 5 % of max exact pressure | 0.09 % / 0.06 % (radial), theta scatter is round-off: PASS |
| RT-B B4 | net outer flux at most 1e-4 of wall flux scale | 1.7e-6 / 1.8e-6: PASS |
| RT-B B5 | all steps converge, mean Newton iterations at most 12 | 200/200, mean 1.86: PASS |
| RT-C | only if RT-B passes all gates | RT-B failed; a later coordinator authorisation (pre-registered in the plan appendix before the diagnostic run) allowed one stall diagnostic and RT-C M1 then M2 if the stall resolved |

Declared classification rules (plan, final section): VERIFIED AGREEMENT needs mesh and time convergence and was declared not reachable in this pilot "by construction" (no dt study, no M3 for RT-B). The static RT-A test has no time dependence and has M1-M3; its VERIFIED AGREEMENT label is applied to the static solid only (see the inconsistency note). DIAGNOSTIC DISCREPANCY requires a failed gate and a quantitative mechanism prediction confirmed within +-25 %.

## Deviations, amendments and post-hoc steps (all labelled)
- A1 (before any scientific RT-B/RT-C run): `min_residual` 1e-18 instead of 1e-30 and the unrecognised `dtol` tag removed from the fluid-module solver block, because the C-infinity ramp starts at round-off and FEBio's relative norms cannot be met at amplitudes of 1e-30.
- A2: diagnostic-only generator options (`--solid svk`, `--three-field 0`, `--clamp-outer-mesh`, `--no-solid`, `--mu-f`, `--bc-time-shift`), never used for a gated run.
- Generator defects found and fixed (the failed runs are preserved): duplicate node-set name (jobs 14190, 14191), unrecognised tag (14197), fluid-FSI material id (14206), law1 test design (14186, 14188).
- **Post hoc, not a gate:** the dt/2 time-label offset was found after RT-B (no pre-registered check existed); the aligned numbers are a post-hoc reanalysis, confirmed by the single permitted diagnostic run 14211 (input table advanced by dt/2).
- A3 and A3b (coupled runs): rho_inf = 0.9 and a reformation cap of 25, adopted after the stall diagnostic although its sub-criterion (mean at most 6 iterations) was missed (mean 6.95): the coordinator's rule keyed on the stall, and the RT-C gate itself is at most 12 iterations.
- Differences between the earlier reassessment document (experiment 1 definition) and the executed plan: dt 0.01 s with 200 steps (reassessment: 0.005 s, 400 steps); RT-A twist amplitude 0.5 rad, five times the member's amplitude, to make the finite-strain and bulk terms visible (reassessment: 0.1 rad); the plan's analytical traction prediction (11.7 / 6.0 / 3.1 %) replaced the reassessment's cruder 7 / 3.5 / 1.7 % while the reassessment's gates were kept unchanged.
- RT-C M1 (job 14213) ran with the heavy-lock rule; no M2.
