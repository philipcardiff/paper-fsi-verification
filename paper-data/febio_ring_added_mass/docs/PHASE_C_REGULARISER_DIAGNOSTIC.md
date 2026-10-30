# PHASE C REGULARISER DIAGNOSTIC (one run: E_mesh 1 Pa -> 0.001 Pa)

Status: ONE additional simulation, `runs/phaseC/C2_Emesh1e-3_wet_f1` (Slurm 13712, 1 core, partition main, pristine v4.13, 273 s, 120 MB). Total project-Phase-C runs: 9 (8 + this one). Does not establish Category A.

## 1. Controlled change
- Baseline: C0 (`runs/phaseC/C0_wet_f1_wallslip_ifslip`), M1, free-slip wall + free-slip interface (inviscid-limit formulation), K = 1e8, rho_f = 140, rhoi = 1, 100 steps/period, compatible IC v2 (per-node), tolerances 1e-6/1e-6/1e-6/1e-8/1e-10.
- The generator was first re-run with `--e-mesh 1` and reproduced the stored C0 `.feb` **byte-for-byte** (cmp identical), then re-run with `--e-mesh 0.001`. `diff` of the two files is **exactly one line**: `<E>1</E>` -> `<E>0.001</E>` in the fluid-FSI regulariser solid. No other parameter, mesh, IC, solver or tolerance changed.
- Dry reuse: `runs/phaseB/B1_dry_hex8_f1` has only the "isotropic elastic" ring material (no fluid-FSI domain, so no regulariser), same mesh (4x48), same steps/period (100), periods (4), rhoi, tolerances and max_refs as C0 (`model_meta.json` compared). Dry omega = 2.62106 is unchanged, so only the wet run needed repeating.

## 2. Pre-registered interpretation (written before post-processing)
(a) offset removed: |ratio error| < ~0.05 %; (b) changed but not eliminated; (c) negligible (< ~0.01 pp change); (d) numerical problem (non-convergence, J < 0, phase/amplitude growth). Hypothesis: the E_mesh = 1 Pa stiffening causes ~ -0.16 % in wet omega (Phase B D4 measured -0.160 % in the no-slip family).

## 3. Results (`results/phaseC_regulariser_diagnostic.csv`, `results/phaseC_time_iter_sensitivity.csv`)

| quantity | C0 (E_mesh = 1) | C2 (E_mesh = 1e-3) | change |
|---|---|---|---|
| wet omega (fit) | 1.854296 | **1.851254** | **-0.1640 %** (zero-crossing fit 1.851254, agrees) |
| wet omega error vs 1.80451344 | +2.759 % | +2.591 % | mesh-dominated at M1 |
| wet/dry ratio (dry 2.621064) | 0.707459 | **0.706299** | |
| ratio error vs 0.70625572 | +0.1704 % | **+0.0061 %** | **-0.1643 pp** |
| mu_eff (effective added-mass ratio) | 0.9980 | 1.0046 | +0.66 pp |
| n=2 pressure amplitude | 0.32543 Pa | 0.32489 Pa | -0.16 % |
| n=2 m_a,eff / analytic, layer 0 / 1 | 0.9783 / 0.9407 | 0.9783 / 0.9407 | none visible |
| n=2 pressure phase | 178.1 deg | 178.2 deg | ~0 |
| integrated n=2 ring load / analytic at measured omega | 2.684e-2 / 2.742e-2 N (-2.12 %) | 2.679e-2 / 2.737e-2 N (-2.12 %) | none |
| interface leakage rms / max (relative) | 1.045 % / 1.253 % | 1.045 % / 1.253 % | none |
| net leaked flux (peak / peak ring flux); peak leaked volume | 0.432 % ; 1.650e-6 | 0.433 % ; 1.649e-6 | none |
| mean(J-1) over fluid, max over time | 5.07e-10 | 5.09e-10 | none |
| max element |J-1|; min/max J-1 | -2.9e-9 / +3.9e-9 | -2.9e-9 / +3.9e-9 | none |
| max nodal displacement; max relative velocity | 8.06e-4 m ; 6.151e-3 | 8.08e-4 m ; 6.151e-3 | +0.16 % |
| amplitude loss per period | 3.1e-6 | 3.1e-6 | none |
| Newton iterations/step | 4-5 (median 4) | 4 (min/med/max 4/4/4) | |
| dilatation margin worst (achieved/tol) | 0.48 | 0.40 | better |
| residual / energy / displacement / velocity worst margins | all <= 1e-4 of tol (per C0) | residual 1.0e-5, energy 1.4e-10, disp 1.5e-4, vel 1.4e-4 | similar |
| Jacobian / NaN / max-refs events | 0 | 0 | none |
| wall time / RSS | 256 s / 114 MB | 273 s / 120 MB | load-dependent |

## 4. Energy
- Total KE at the 9 zero crossings over initial: C2 `[1.0, 1.0036, 0.9992, 1.0037, 0.9999, 1.0035, 0.9993, 1.0037, 0.9999]` (C0 `[1.0, 1.0002, 0.9996, 0.9995, 0.9995, 1.0, 0.9999, 1.0001, 0.9997]`).
- Integer-period (same-sign) points: 0.9992-0.9999 in C2 versus 0.9995-1.0000 in C0: conserved to <0.1 %, no drift.
- **New observation**: with E_mesh = 1e-3 the half-period KE is consistently 0.36 % high (alternating pattern). It is absent at E_mesh = 1. It has no visible effect on frequency (window spread 2.2e-5 %), pressure or force. I did not investigate its origin: candidates are a different share of mesh-motion kinetic content in the lumped-mass measure or the non-modal IC cross-term being less damped. Total (kinetic + strain) energy was not independently computed; only KE is logged. Treated as an open, bounded (0.4 %) item.

## 5. Classification of outcome
- **(a) removes the systematic offset: yes at M1**. The shift (-0.164 pp in the ratio) equals the previously measured regulariser effect in the no-slip family (D4: -0.160 % in omega), and the remaining ratio error is +0.006 %. The mechanism is confirmed as a stiffening of the wet ring by the mesh-regulariser solid, as hypothesised.
- (b)-(c) do not apply. (d): no numerical problem: converged in 4 iterations/step everywhere, Jacobians sane, no drift, margins as good as or better than baseline.
- Everything else in the n=2 pressure/force/leakage diagnostics is **unchanged**, so the regulariser is not the cause of the -2.1 % force deviation or the 1 % leakage at M1.

## 6. Caveats (do not skip)
1. **M1 ratio near zero rests on cancellation.** Both dry and wet omega are still +2.6 % high at M1; the ratio agrees only because spatial errors cancel (they cancel to 0.0075 pp between M1 and M2, per Phase C). The result that matters is the ratio at M2/M3 with E_mesh reduced, not tested.
2. A single regulariser value was tried. Linear scaling (effect ~ E_mesh) predicts a residual of ~0.0002 % at 1e-3 but was not verified; the +0.006 % residual is within the roughly 0.01 pp resolution of the existing time/K/mesh effects and should not be read as a measured remnant.
3. The n=2 force still deviates -2.1 % from m_a*u_r at M1 (halving per refinement, and equal to the layer-0 exact-profile sampling offset: -2.107 %, -1.074 %, -0.542 % at M1-M3 vs measured -2.121, -1.065, -0.528 %), i.e. it appears to be a measurement/sampling effect of the pressure used (layer-0 element centre), not a force error. This inference is from existing M1-M3 data, not a new test.
4. Interface leakage (1.0 % rms at M1, 0.32 % at M3) is a weakly enforced condition and remains bounded rather than eliminated.
5. Iterative error not quantified (tight-tolerance attempt failed on a dilatation floor); viscosity sensitivity untested; this is the inviscid-limit formulation, not the finite-viscosity benchmark.
6. No Category A claim. E_mesh = 1e-3 Pa is a numerical regularisation parameter, so lowering it moves the model toward the unregularised physical problem rather than adjusting physics; this is the justification for selecting it, not that it improved agreement.

## 7. Figures reviewed
All six Phase C figures were opened and corrected (presentation only; CSV numbers unchanged: `results/phaseC_mesh_convergence.csv` re-compared with the pre-edit printout). The list of fixes is in `PHASE_C_VERIFICATION_RESULTS.md` section 7.

## 8. M2 confirmation (one run, `runs/phaseC/C2_Emesh1e-3_wet_f2`, Slurm 13713, 1 core)
Same controlled change at M2 (regenerated baseline `.feb` identical to `C1_M2_wet_f2`; diff = the one `<E>` line; dry D1 reused, no regulariser). Data: `results/phaseC_regulariser_diagnostic.csv`.

| | M2, E_mesh = 1 | M2, E_mesh = 1e-3 | change |
|---|---|---|---|
| wet omega | 1.819071 | 1.815972 | -0.170 % |
| wet omega error | +0.807 % | +0.635 % | |
| ratio (dry 2.571081) | 0.707512 | 0.706307 | |
| ratio error | +0.178 % | **+0.0072 %** | -0.171 pp |
| mu_eff | 0.9977 | 1.0045 | |
| n=2 m_a,eff layer 0 / 1 | 0.9889 / 0.9686 | 0.9889 / 0.9686 | none |
| leakage rms / max | 0.40 % / 0.49 % | 0.40 % / 0.49 % | none |
| mean(J-1) max | 5.3e-10 | 5.3e-10 | none |
| KE at integer periods | 0.9994-0.9998 | 0.9991-0.9999 | |
| KE half-period excursion | none | up to +0.09 % (M1: +0.36 %) | shrinks with refinement |
| iterations/step | 4-5 | 4-5 | |
| dilatation margin worst | 0.44 | **0.76** | smaller margin; still met |
| wall / RSS | 1493 s / 342 MB | 1597 s / 366 MB | |

- The regulariser effect is mesh-independent: -0.164 pp at M1, -0.171 pp at M2. The ratio error with E_mesh = 1e-3 is +0.0061 % (M1) and +0.0072 % (M2): the offset is removed within the resolution of the other effects (~0.01 pp), and the ratio is flat over M1 -> M2, so M1's agreement was not a cancellation artefact.
- The wet absolute frequency error (+0.635 % at M2) is spatial (dry +0.628 %) and cancels in the ratio, as before.
- The half-period KE excursion shrinks with refinement (0.36 -> 0.09 %), so it is a coarse-mesh effect, not a defect.
- Dilatation tolerance margin is thinner (0.76 of tolerance, vs 0.44); iterative error remains unquantified, and this is one more reason to treat that as open.
- E_mesh was tested at one reduced value only; no E_mesh sweep or M3 repeat (not needed for the decision).
