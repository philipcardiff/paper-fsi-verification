> **ARCHIVED WORKING REPORT (verbatim except path/host sanitisation).** Written during the investigation; file paths such as `runs/`, `tools/`, `cavity_agent/` refer to the original working tree, not to this package (see `docs/INDEX.md` for the package equivalents). Later reports supersede earlier ones where they disagree; known numerical/wording inconsistencies and their resolution are listed in `provenance/INCONSISTENCIES.md`. Nothing here claims a verified FEBio result.

# FEBIO FLUX-PRESERVING PROJECTION (round 2), 2026-10-30
Derivation: `FEBIO_PRESSURE_PROJECTION_DERIVATION.md` section 13. The alpha = 3/4 form of round 1 FAILED mass conservation and stays rejected (Q_out/Q_in - 1 = -1.55 % on C1b).
## A. Pre-registration (before any round-2 run)
Formulation: redistribution (a) = (c) of section 13.3; alpha = 1.0e-3 (section 13.5); fixed-mesh fluid domain only. Budget for the round: 0.75 CPU-h (running table in `FEBIO_EXPERIMENTAL_BUILD.md`, section "flux-preserving round").
Element/assembled tests (must pass before any run): sum of free rows of the stabilisation residual = 0 to machine precision on an assembled random mesh with constrained outlet/rim/pins (old form: nonzero = -sum constrained rows); constrained rows contribute 0; dS/de and dS/dx vs finite differences <= 1e-8; constant p -> 0; alpha = 0 bitwise no-op; parity energy as before.
Tiny exact Poiseuille channel (the stored T3 geometry/settings, 1400 hex8, K 1e6, mu 0.01, dt 0.5, t = 100 s, 2 cores), stock baseline = stored `runs/cavityExp/T3_stock` (no duplicate). Criteria (STOP at the first failure; no second alpha):
| # | criterion | threshold |
|---|---|---|
| P0 | rim-once with alpha = 0 reproduces the stock steady state (separate control run) | fields at t = 100 within 1e-8 relative of T3_stock |
| P1 | mass flux: Q(x)/Q_in - 1 at inlet, every Mon section (x = 2, 7, 10, 12, 13) and outlet | max abs < 1e-4 (prediction: 1e-5 interior, ~0 outlet) |
| P2 | pressure gradient (fit x = 2..13) vs the discrete Poiseuille value -0.12/(1 - h^2) = -0.121212 | within 0.2 % |
| P3 | velocity profile at x = 7, 10, 12 vs the nodal inlet profile (the exact nodal parabola of the actual inlet flux) | max abs difference <= 0.5 % of the centreline value |
| P4 | y-alternation of p (rms of (p_{j-1} - 2 p_j + p_{j+1})/4 on the columns; stock 6.3e-4 Pa) | <= 6.3e-5 Pa (10-fold reduction) |
| P5 | independent re-assembly (python, straight channel): momentum residual and continuity residual including the PPP~ term | relative <= 1e-5 |
| P6 | nonlinear convergence | normal termination, no reformation-limit event, iterations/step <= 2 x stock |

## B. Tests (all on the final round-2 library, sha256 6627dadd...; binary febio4 85c344a7... unchanged)
| test | result |
|---|---|
| assembled random mesh 7x4x1 (perturbed nodes, random e_f, 16 constrained nodes: outlet column, rim corners, 2 pins), alpha = 1e-3, library routines | sum over ALL rows old -4.4e-20, new -2.9e-20; **sum over FREE rows: old -9.7e-6 (= -sum of constrained rows, 80 % of max|S|), new -2.9e-20 (2e-15 of max|S|)**; constrained rows of the new form 0; PROVEN (machine precision) |
| 40 random elements/masks incl. all-free and all-constrained | free-row sum 1e-14 of max; dS/de vs FD 1.2e-12; dS/dx vs FD 3.2e-9; constant p 1e-16; alpha = 0 with mask bitwise unchanged |
| regression of round-1 element test | unchanged (parity energy exact 2e-16, FD 2.5e-13 / 3.5e-10) |
| P0 rim-once with alpha = 0 (T5_rimonce_a0) vs stored stock T3 at t = 100 | **PASS: 5.9e-12 relative on all Mon files** (same fixed point); 274 vs 492 iterations, 39 s vs 61 s. So the rim change alone (separated from the projection) does not change the converged state |
## C. Tiny exact Poiseuille channel, alpha = 1.0e-3 (job 14119, `runs/cavityExp/T5_ppp_a1e-3`), baseline stored stock T3 (job 14105), round-1 alpha 3/4 (job 14106) for comparison
| criterion | stock | alpha = 3/4, round 1 | alpha = 1e-3, round 2 | verdict |
|---|---|---|---|---|
| P1 max abs (Q/Q_in - 1) at inlet, x = 2,7,10,12,13, outlet (nodal column flux, the stock measure) | 1.7e-6 | 7.4e-3 | **4.86e-4** (x2 -4e-6, x7 1.5e-5, x10 7.5e-5, x12 1.9e-4, x13 2.9e-4, outlet 4.9e-4) | **FAIL (> 1e-4)** |
| P2 dp/dx fit vs -0.121212 | -0.04 % | -0.79 % | -0.29 % | FAIL (> 0.2 %) |
| P3 profile vs inlet nodal profile (max, of centreline) | 6.2e-5 | 7.4e-3 | 2.0e-4 | pass |
| P4 y-alternation rms max over sections (Pa; stock 6.3e-4) | 6.3e-4 | 1.2e-5 | 3.0e-4 at x12 (6.0e-5 at x7, 4e-6 at x2) | FAIL (> 6.3e-5 at x >= 10: only 2x reduction there) |
| P5 independent re-assembly: momentum / continuity incl. PPP~ | 1.6e-9 / 1.4e-7 | n/a | 1.5e-9 / 1.43e-7 (without the term 7.6e-5); sum S~ over free rows 1e-18 | pass |
| P6 convergence | 492 it, normal | 346 it | 274 it (1.37/step), normal, no reformation-limit | pass |
**Decision by the pre-registered rule: STOPPED at the first failed criterion (P1). No second alpha, no tuning.**
### Diagnosis of the P1 failure (login-node post-processing of the same run, no extra computation)
* The stabilisation behaves as derived: cumulative sum of the redistributed term S~ upstream of x = 2..13 is -1.0e-5 ... -1.8e-5 (prediction alpha (h/H)^2 = 1.0e-5), the sum over free rows is 1e-18, and the solved rows satisfy C + b + S~ = 0 (1e-7). The row-balance quantity is therefore conserved to ~2e-5, as predicted: the element-pair-mean flux error, the quantity the continuity rows actually constrain, is **2.0e-5 max** (< 1e-4).
* The nodal column flux Q_j is NOT constrained individually: the rows fix only Q_j + Q_{j+1} (the cut element sees the mean flux of its two node columns). The remaining freedom is an **x-sawtooth in the nodal flux** (alternating sign Q_j - 1: +40e-5, -44e-5, +45e-5, -48e-5, +49e-5 for j = 136..140), growing from 2e-5 at the inlet to 4.9e-4 at the outlet; the stock solution has it at 2e-6. So the failed measure is an x-parity mode of the velocity/flux field of the unstabilised x-parity family, amplified ~300x relative to stock when the pressure parity is damped at this alpha; the cause of the amplification is UNRESOLVED (hypothesis: the stock x-sawtooth is held small by the coupling to the pressure x-sawtooth that the term now damps; not tested). The y-alternation also persists near the outlet (3e-4 Pa).
* Therefore the pre-registered criteria are failed; the failure is not an assembly or conservation-bookkeeping error. Conservation in the global/row sense: PROVEN to machine precision (formulation (a)). Nodal-column flux: FAILED.
### CPU-h of the round
Rebuild 14117 (2 x 6 s), T5 rim-once control 14118 (2 x 40 s), T5 alpha = 1e-3 14119 (2 x 36 s): 164 CPU-s = **0.046 CPU-h of 0.75** (from sacct). The element/assembled tests ran on the login node (seconds, not counted).
## D. Labels and answer on a later fluid-only cavity test
* Global mass conservation via redistribution (a)=(c): PROVEN (machine precision). Poiseuille cut defect alpha (h/H)^2 cannot vanish for local P0 PPP: PROVEN (derivation) and observed in magnitude (1e-5..2e-5 pair-mean). Pre-registered nodal flux/gradient/parity criteria: FAILED. Rim-once separated: PROVEN harmless at alpha = 0 on this case. FSI domain: compiled/element-tested only, NOT TESTED in a run.
* A fluid-only cavity test is **not justified now**: the straight channel, the cheapest exact case, already fails three pre-registered criteria for this formulation at its stability-derived alpha, and the failure involves an x-parity flux mode that a cavity run (about 0.05 CPU-h per C1b-size run at M1 with this solver, 0.2 CPU-h for the stock comparator) would only reproduce with less diagnostic clarity. Prerequisite: understand/suppress the x-sawtooth (for example a derived stabilisation acting on the x-parity velocity/flux family, or a jump-based/macro-element projection that vanishes on linear p), verified first on the tiny channel (about 0.02-0.05 CPU-h per run).
