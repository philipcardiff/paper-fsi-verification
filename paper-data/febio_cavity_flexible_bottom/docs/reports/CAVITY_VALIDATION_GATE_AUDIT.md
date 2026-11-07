> **ARCHIVED WORKING REPORT (verbatim except path/host sanitisation).** Written during the investigation; file paths such as `runs/`, `tools/`, `cavity_agent/` refer to the original working tree, not to this package (see `docs/INDEX.md` for the package equivalents). Later reports supersede earlier ones where they disagree; known numerical/wording inconsistencies and their resolution are listed in `provenance/INCONSISTENCIES.md`. Nothing here claims a verified FEBio result.

# CAVITY VALIDATION-GATE AUDIT (coordinator, 2026-10-30)
Read from the CAVITY agent's own documents (`CAVITY_FLUID_ONLY_RESULTS.md`, `CAVITY_FSI_PILOT_RESULTS.md`, `CAVITY_NUMERICAL_SENSITIVITY.md`) and the ledger. No check is re-marked.

## 1. Gate record
| Gate / check | Pre-registered criterion | Outcome (as recorded by the agent) | Status I record |
|---|---|---|---|
| C1 (run 14012) A3 steady Fy_raw (<1e-4) | steady over last 10 % | 8.2 % drift from undetermined pressure at the 4 floor-end corner nodes | **FAIL** (preserved) |
| C1 A4 upstream dp/dx within 2 % | raw fit, 8 columns | -0.1957 vs -0.1212 (+61 %) | **FAIL** |
| C1b (14015) A3, A5-A8 | as above | pass, after pinning corner dilatation to 1.15 Pa | pass |
| C1b A4 | as above | raw -0.1945 (+60 %); same data pair-averaged -0.3 %..+0.0 %; agent diagnosed an aliasing defect in the criterion | **FAIL as written** (criterion defect diagnosed, not a pass) |
| C1b raw interface force | not a pass/fail item | Fy_raw -4.6637 N/m, odd-even bracket +/-2.8 % | uncertainty carried |
| C2 pilot P1 normal termination | 0 errors | exit 1 at t = 975 s (`Max nr of reformations reached`, dilatation-increment floor) | **FAIL** |
| C2 P2 steady (spread < 1e-4 over [1125,1250]) | window never reached | spread 7.6e-5 over its last 10 % | window not met; conditional steady estimate only |
| C3 M2 P1 | normal termination | exit 1 at t = 872.5 s, same mechanism | **FAIL** |
| C3 M2 steady spread | < 1e-4 | 1.79e-4 relative on u_y(A) over the last 10 % | **not met** (above the criterion) |
| P3 flux / mass | 1e-3 | 1.2e-4 (M1), 1.9e-5 (M2) | pass |
| S2 pressure-transfer check | within 2 % | 0.025 % | pass (see section 3) |

## 2. Why the coupled runs proceeded
The CAVITY agent states (CAVITY_FLUID_ONLY_RESULTS.md, "Gate decision") that its own pre-registered rule (A3-A8 must pass) was **not** satisfied, that A4 failed, and that it proceeded to one coupled pilot "as a diagnostic ... not counted as a verification result". So: the gate was not retroactively marked as passed; the failure and the decision to proceed anyway are both recorded. Two points the record does not support:
1. The earlier explicit gate from your previous brief ("no coupled pilot until fluid-only passes") was replaced by the new brief, whose C2 follows C1 without requiring A4; the agent's own stricter rule was relaxed by the agent, with a stated reason (the A4 criterion was ill-conditioned for a parity-unbalanced window) but **no independent confirmation that the pressure field is accurate**, because the raw-force uncertainty (+/-2.8 %) could not be reduced against a reference.
2. The A4 explanation (aliasing of the odd-even mode into a short-window slope) is plausible and numerically consistent (pair-averaged gradient agrees to 0.3 % or better) but was derived from the same data; it was not tested by a run with a parity-balanced window prescribed in advance. I therefore keep A4 as FAIL and the force comparison as unverified, and the C2/C3 results as "diagnostic".
3. C2 and C3-M2 **did not terminate normally** (P1 failed) and M2 did not meet P2. They were nevertheless analysed and reported, correctly labelled PARTIAL; they must be called conditional steady estimates, not clean completions.

## 3. S2: what it actually tested
Source: CAVITY_FSI_PILOT_RESULTS.md section P4/S2 and `S2_plate_hex8_fsiM2pressure.feb`. The plate-only hex8 80x8 StVK model was loaded **only by the face-mean nodal interface pressure** of C3_fsi_M2 at the last state (follower pressure). Result: u_y(A) -0.226176 vs coupled -0.226233 (0.025 %). Therefore:
- It tested **pressure traction alone**, not pressure plus viscous traction. The viscous load (Fx_visc about 0.02 N/m, at most 0.4 % of F_y) is absent; u_x(A) already differs (plate-only 3.92e-4 vs coupled 3.35e-4 m, 17 %).
- It is a **self-consistency test**: it feeds FEBio's own fluid pressure to the same structure. It shows that the deflection is explained by that pressure field, i.e. that the `fluid-FSI traction` load on the plate is not wrong by more than 0.03 % *in the vertical, pressure part*. It **cannot** detect an error in the pressure itself (level, boundary condition, odd-even mode), nor a wrong total force: a pressure field 5 % too low would pass identically.
- Face-averaging of the nodal pressure cancels most of the odd-even mode, so S2 would not reveal an alias in the raw nodal force either.
- It does not test force normalisation or the plate-equilibrium reaction (the P4 plate-equilibrium test was pre-registered and is recorded by the agent as **NOT verified** because reaction logging gave zeros in the solid module and no xplt reader existed; the xplt audit requested by the coordinator addresses this).

## 4. Consequence
The force deficit (-5.6 % at M2, -7.9 % at M1, vs solids4foam at the same dx) is currently **unaccounted for**: S2 excludes the vertical pressure-load transfer as the cause, which leaves the pressure field itself (boundary conditions, discretisation, odd-even mode) and the physical-model difference. See `CAVITY_PRESSURE_FORCE_AUDIT.md`.
