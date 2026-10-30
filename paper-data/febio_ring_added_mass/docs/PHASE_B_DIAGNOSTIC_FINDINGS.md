# PHASE B DIAGNOSTIC FINDINGS

Date: 2026-10-29. Pristine FEBio 4.13.0.32ae206ff (sha256 prefix `f5ae1ed6f347`) for every run below; `rhoi = 1`; 1 core per run; main/dev partitions; per-node fluid-IC workaround. Full per-run data: `results/diagnostic_matrix.csv` (built by `tools/build_matrix.py` from each run's `report.json`; nothing hand-entered). Audit: `PHASE_B_COUPLING_AUDIT.md`. Acceptance thresholds are those of `PHASE0_AUDIT.md` §9 and were not changed.

## 1. Outcome in one paragraph

A controlled comparison demonstrates that **the strong nodal no-slip constraint on the rigid outer wall is the dominant cause of the −6 % wet/dry-ratio error and of the non-physical pressure** in the original model. Replacing it with FEBio's natural (frictionless) wall condition, and changing nothing else except the matching initial-velocity coefficient, moves the ratio error from **−6.01 % to −1.11 %**, the effective added-mass ratio from 1.270 to 1.050 (exact 1.008), removes a spurious uniform pressure (3.9 Pa → 0.000 Pa at the fundamental frequency), and brings the integrated kinetic energy to 0.9986–0.9995 of its initial value at integer periods. **The 0.3 % screening tolerance is not met by any run** (best −1.11 %), so stopping condition 2 is not satisfied; condition 1 (root cause demonstrated) is satisfied for the main effect, and condition 3 (five-run budget) is exhausted. What remains of the −1.1 % is not explained.

## 2. Runs (all 400/400 steps, normal termination, no Jacobian/NaN/inversion messages, 1 core)

Counting against the five-run budget: R1 (baseline rerun after the audit defect), D1, D2, D4, D5. D3 (hex20 dry ring) was **replaced**: D1 already shows the dry error converges at observed order 2.04, so a hex20 test is lower value, and equal-order hex20 fluid support in fluid-FSI is unverified in the tagged source. Earlier runs B1–B4 are prior evidence (v1 IC), shown for the comparisons.

| run | what differs | mesh (ring / fluid, hex8) | wall time | ω fit (rad/s) | error vs ref | wet/dry ratio | ratio error vs 0.706256 | μ_eff (exact 1.008) |
|---|---|---|---|---|---|---|---|---|
| B1 dry | – | 4×48 | 12 s | 2.62106 | **+2.584 %** | – | – | – |
| B2 wet | v1 IC, no-slip wall, E_mesh 1, K 1e8 | 4×48 / 12×48 | 207 s | 1.72889 | −4.19 % | 0.65961 | **−6.60 %** | 1.298 |
| B3 wet | K 1e7 | same | 210 s | 1.72926 | −4.17 % | 0.65976 | −6.58 % | 1.297 |
| B4 wet | K 1e9 | same | 244 s | 1.72701 | −4.30 % | 0.65890 | −6.71 % | 1.303 |
| **R1 wet** | corrected IC (v2); otherwise = B2 | same | 260 s | 1.73982 | −3.58 % | 0.66378 | **−6.01 %** | 1.270 |
| **D1 dry** | mesh factor 2 | 8×96 | 48 s | 2.57108 | **+0.628 %** | – | – | – |
| **D2 wet** | mesh factor 2; v2 IC; no-slip wall | 8×96 / 24×96 | 1325 s | 1.76738 | −2.06 % | 0.68741 (vs D1) | **−2.67 %** | 1.116 |
| **D4 wet** | E_mesh 1e-3 (from 1 Pa); else = R1 | 4×48 / 12×48 | 300 s | 1.73703 | −3.74 % | 0.66272 | −6.16 % | 1.277 |
| **D5 wet** | **outer wall free-slip**; v2 IC adjusted for it; else = R1 | same | 239 s | **1.83067** | +1.45 % | **0.69845** | **−1.11 %** | **1.050** |

Wet/dry ratios use the dry run with the same ring mesh, solver, integrator, tolerances, IC and analysis (B1 for the 4×48 wet runs, D1 for D2). Frequencies are damped-harmonic least-squares fits over 0–4T; the zero-crossing diagnostic agrees to ≤ 1e-4 relative in every run and the window spread (0–4T, T–4T, 0–3T, 0–2T) is ≤ 0.026 %. Input hashes, job ids and all other fields are in the CSV.

## 3. Separation of effects (model form, discretisation, solver)

Ratio-error changes isolate one change at a time (percentage points of the wet/dry ratio):

| effect | comparison | change in ratio error | class |
|---|---|---|---|
| initial-velocity incompatibility (audit F1) | B2 → R1 | −6.60 → −6.01 (**+0.59**) | input/IC defect |
| fluid-mesh regulariser stiffness 1 Pa → 1e-3 Pa | R1 → D4 | −6.01 → −6.16 (**−0.15**) | model-form parameter; ω shifts 0.16 % |
| compressibility, K 1e7 / 1e8 / 1e9 (v1 IC) | B3 / B2 / B4 | −6.58 / −6.60 / −6.71 (spread 0.13) | model form; **not** the cause |
| **rigid outer wall: no-slip (nodal w = 0) → free-slip** | **R1 → D5** | **−6.01 → −1.11 (+4.90)** | **model form / discretisation of the BC** |
| spatial refinement ×2, no-slip wall | R1 → D2 | −6.01 → −2.67 (**+3.35**), two levels: apparent order 1.17 | spatial (no-slip artefact shrinks with h) |
| spatial refinement ×2, dry ring | B1 → D1 | +2.584 → +0.628 %, **order 2.04** | spatial (structural) |
| time step | not varied | – | time error estimated, not measured: midpoint-rule phase error ≈ −(ωΔt)²/12 ≈ −3e-4 (4e-5 in the ratio) [INF] |
| nonlinear/linear solver tolerance | not varied | – | all runs converged every step to the same tolerances; effect **not measured** |

The wall effect and the refinement effect point the same way: a strong nodal no-slip constraint over a cell that is ~100× thicker than the physical Stokes layer (δ = √(2ν/ω) ≈ 3.3e-4 m against a 3.75e-2 m cell) is a discretisation artefact that shrinks with cell size (D2) and disappears when the constraint is removed (D5).

## 4. Pressure, interface force, energy, damping, volume

Pressure decomposed into angular components in the fluid layer next to the ring (`tools/pressure_modes.py`; n = 2 is the physical signal, analytic amplitude m_a·ω²·D):

| run | n = 2 amplitude at ω (Pa) | n = 2 residual rms (Pa) | **n = 0 amplitude at ω (Pa)** | n = 0 at 2ω (Pa) |
|---|---|---|---|---|
| B2 (v1 IC) | 0.290 | 0.875 (period-2 sawtooth) | 1.22 | 0.87 |
| R1 | 0.325 | 0.185 | **3.91** | 1.04 |
| D4 | 0.325 | 0.185 | 3.92 | 1.05 |
| D2 (f2) | 0.324 | 0.046 | 2.19 | 0.55 |
| **D5** | **0.336** (predicted 0.328) | **0.011** | **0.000** | **0.026** |

- n = 2 pressure at the added-mass relation: effective added mass / analytic, layer 1: 1.004 (B2), 1.007 (R1), 0.999 (D2), 0.982 (D5); layer 0 is 2 % to 5 % higher (cell centre closer to the interface, ≈2 % by the finite-difference estimate); phase 178.3–178.4° (anti-phase, correct sign).
- Integrated n = 2 radial loading from layer 0 against the analytic added-mass loading: R1 2.68e-2 vs 2.55e-2 N (+5 %); D5 2.77e-2 vs 2.71e-2 N (+2.3 %).
- The spurious uniform (n = 0) pressure at the fundamental frequency is present with the no-slip wall (1.2–3.9 Pa, 12× the physical signal in R1) and **vanishes with the free-slip wall**. Linear theory decouples it from the n = 2 mode; its mechanism (a net discrete volume change at the over-constrained wall, K-amplified) is an inference. D5's only n = 0 content, 0.026 Pa at 2ω, is the physically expected second-order effect.
- Fluid volume: max |mean(J−1)| ≤ 8e-8 in every run (K·that = the n = 0 pressures above); max |J−1| ≤ 4e-8 except B3 (3.6e-6).
- Amplitude loss per period (rhoi = 1, fitted decay): ≤ 2.1e-4 in all wet runs, 1e-9 dry; B4 grows by 2.7e-4 (far below s4f's growth tolerance 2e-3). **No numerical damping at this level.**
- Energy: total kinetic energy at displacement zero-crossings of the same sign (integer periods): D5 0.9986, 0.9995, 0.9987, 0.9990 of initial (conserved to 0.14 %); R1/D4 drift to 0.96–0.99; D2 0.991–0.997; the earlier defective-IC runs are worse: B2 0.82–0.86, B4 0.81–0.93 and **B3 (K = 1e7) grows to 2.68× its initial kinetic energy by 3T**, i.e. that configuration is not merely inaccurate but energetically unstable. The alternating +6–8 % at half-periods is an artefact of the corrected IC's stationary field (its cross term with the reversing potential flow changes sign each half period), present in every v2 run, and is not energy gain.
- Ring velocity consistency (rms of v − du/dt as % of peak): dry 0.09 %, D5 0.37 %, D2 0.27 %, R1/D4 1.1 %, B2 6.7 %.
- Solver evidence: every accepted step converged; average residual evaluations per step 3.8 (dry) to 4.0 (wet) with full Newton; no inverted-element or Jacobian messages; no "Max nr of iterations reached" events. Tolerances unchanged (dtol=vtol=ftol 1e-6, etol 1e-8, rtol 1e-10). Convergence tolerances were not swept.

## 5. S1–S6 screening status (thresholds unchanged)

| gate | result |
|---|---|
| S1 dry ω within 0.5 % | **FAIL** at hex8 factor 1 (+2.58 %) and factor 2 (+0.63 %); second-order convergence observed (2.04); factor 4 would be expected near 0.15 % [INF], not run |
| S2 wet ω within 0.5 %, ratio within 0.3 % | **FAIL** in every run (best D5: ω +1.45 % because the dry ring is 2.58 % stiff; ratio −1.11 %) |
| S3 K convergence | **INDETERMINATE/untested under the corrected model**: K was varied only with the defective v1 IC and no-slip wall |
| S4 amplitude loss ≤ 0.5 % per period | PASS in all runs |
| S5 no acoustic/pressure ringing | FAIL for all no-slip-wall runs; **PASS for D5** (n = 0 at ω absent, n = 2 residual 3 % of signal) |
| S6 archived | PASS |

## 6. What is established, what is not

Established by controlled comparison: the IC defect (small effect, +0.6 pt); the regulariser stiffness is not responsible (−0.15 pt); compressibility over K = 1e7–1e9 is not responsible; the strong no-slip on the rigid outer wall causes most of the error and the spurious pressure; the dry hex8 ring converges at second order.

**Not established:**
- The origin of the remaining −1.1 % at D5. Candidates, none tested: the no-slip condition at the **interface** (kept per your instruction), spatial error of the wet mesh (D2-type refinement with the wall free-slip was not run), K, dt, solver tolerance, the ring's own residual.
- Mesh convergence of any wet quantity: only two levels exist, and not for the corrected (free-slip) model.
- Temporal and solver-tolerance errors.
- Whether the free-slip rigid wall is the faithful analogue of the reference. The reference's own statement is that the unresolved viscous layer makes the computed flow effectively inviscid and the inviscid values are the references (viscous frequency shift −0.026 % at this density), which supports free-slip as the consistent inviscid limit, but FEBio imposes it weakly (natural condition w_n = 0 on a curved wall); this is a model-form choice that must be disclosed.
- Anything about other densities, other modes, or other benchmarks.

**Statement of scope.** These results do **not** show that FEBio is incapable of modelling this benchmark. They show that, with equal-order hex8 fluid-FSI elements on the 12×48 and 24×96 meshes, a strongly imposed no-slip wall gives a −6.0 % (coarse) and −2.7 % (refined) wet/dry-ratio error, and that a free-slip wall gives −1.1 % on the coarse mesh; the configuration that would settle the remainder has not been run.

## 7. Items that are limitations of this study

The corrected IC is partly non-modal (adds KE, 1.33× analytic) and is a departure from both the solids4foam impulsive start and the v1 field. The wall free-slip run also changed the wall coefficient of that IC (by construction). Two runs at factor 2 only (D1, D2); D2 used 1 core for 22 minutes, the 1-core choice was kept for comparability across runs (Phase A showed thread count changes the Newton path). The Phase 0 hypothesis that compressibility would limit accuracy was wrong for this case; the earlier hypothesis that hex8 locking limits the dry ring was also not supported (convergence is second order).

## 8. Reproduction

`tools/make_ring_feb.py` (flags: `--ic-field compatible`, `--wall slip`, `--e-mesh`, `--mesh-factor`, `--k`), `tools/run_febio_job.sh`, `tools/run_report.py`, `tools/pressure_modes.py`, `tools/audit_ic.py`, `tools/audit_pressure_force.py`, `tools/build_matrix.py`. Run directories under `runs/phaseB/` hold inputs, logs, probes, `report.json` and `run_meta_*.json`. Originals and failed cases are untouched.
