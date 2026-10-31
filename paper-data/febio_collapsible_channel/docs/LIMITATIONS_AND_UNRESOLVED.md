# Limitations and unresolved questions

## What this evidence does not establish
- **The coupled FEBio-vs-oomph-lib comparison is NOT ESTABLISHED.** No full 140-step coupled run was made (the pre-registered stop rule after the F1 gate and the refuted hypothesis H_D). The 20-step diagnostic does not bracket the first trough, so the first-trough value and time, the trough-to-trough period, the rebound and the second trough are NOT ESTABLISHED for FEBio.
- No claim about solids4foam, the manuscript, or the oomph-lib reference beyond the reference values used. No FEBio source change is proposed; no manuscript change is requested.
- One mesh level (fluid 7,040 cells, wall 320x8), one time step (0.025 s), one K (1e7), one rhoi (0.9), one E_mesh (1). No mesh or time-step refinement, no K sweep. VERIFIED AGREEMENT is not available for any quantity.
- The fluid mesh is not the tutorial mesh (7,040 vs 3,200 cells; the only hex8 wall that met its static gate is conformal with 320 columns). Tutorial 80x8 hex8 wall: 11.6 % too stiff (PROVEN).

## Limitations of the evidence that exists
- F1 vs F2 is a reconstruction comparison; which of the two is the physical load is UNRESOLVED (neither family is exact at the nodes). The free-body identification supports that the load the wall receives equals F1, using the S1 static curve (320x8 hex8, 0.76 % stiff) and the static shape psi; it is valid only for |u_mid| <= 0.1425 m (t <= 0.35 s) and its static-shape assumption is poor at t <= 0.15 s.
- Steady F1 gate: the relative spread is measured against a small steady force (design weakness stated in advance; the dynamic test showed the effect is much larger in the dynamic regime).
- Whether the y-parity amplitude falls under refinement (cavity: h^2-like) was not measured here (no second fluid mesh; the cheapest informative extra test, see the recommendation).
- The wall lag matches the load-balance law in sign and size but is not proof of a sole cause: hex8 wall 0.76 % too stiff (period effect about 0.4 %, analytical), rhoi 0.9 and dt 0.025 (oomph itself differs by 6.3 % of the peak at dt 0.025 against its converged value, its dt = 0.025 trough -0.2129 m vs -0.2160 m) and K = 1e7 compressibility ((omega L/c)^2 = 0.4 %) were not separated.
- The load-balance sensitivities are analytical estimates (shape factor, exponents) not re-derived beyond arithmetic.
- ALE robustness for a 22 % collapse was only exercised up to 0.5 s (no negative Jacobian; 20 steps): later collapse/rebound not tested.
- Junction pin rule (b) literal trigger (7e-4 / 5e-4 Pa) was judged a non-event without a run; not tested that pins leave F1 - F2 unchanged for this geometry (FSI1 showed they did).
- The pre-registration text was overwritten in place, no earlier version exists, and its written clock times disagree with the job times (INC-1): ordering of A1 relative to the F1 run cannot be proven from the files.
- Reference: the oomph-lib CSVs were copied from solids4foam; the oomph-lib driver was not re-run; the static value was re-derived only with an independent beam code of the same model class (about 0.05 %); the reference README's "period about 0.8 s" does not match the CSV (INC-2).
- Not independently re-read: the FEBio source statements (`FEFluidFSITraction` nodal pressure, v4.13 has no stabilised pressure operator) are carried from the earlier campaign documents (cavity and FSI1 packages).

## Unresolved questions
1. Does the y-parity amplitude scale with the smooth load (added-mass pressure) or with viscous/convective forcing? D1 (amplitude 58 Pa at t = 0.25 s, growing with the load: D mean 51 -> 203 Pa from t = 0.10 to 0.30 s) suggests it tracks the load, F1 (steady, 0.9 Pa) does not contradict; not separated.
2. What is the FEBio full-run trough, period and decay as the parity-biased answer (not run)?
3. Does a 2x finer-in-y fluid mesh reduce the y-parity?
4. What is the physically correct wall pressure at the interface nodes (no reference available)?
