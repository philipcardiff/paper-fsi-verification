# Load-balance analysis (analytical estimate [A], from BENCHMARK_REASSESSMENT_AFTER_FSI1.md section 4; derived before any run of this pilot)

**Withdrawn assumption.** The earlier roadmap said the collapsibleChannel displacement is "driven by the external load and so insensitive to the force floor / pressure parity". The reassessment WITHDREW that statement as an unfounded assumption: at the first trough the fluid added-mass pressure is the larger part of the net wall load. The pilot was designed around that correction (parity measured before the coupled run). The arithmetic below was re-evaluated here (`scripts/reverify.py`, last rows); it is an estimate from the stored reference peak, not a simulation result.

- Wall: clamped plane-strain beam, L = 10 m, h = 0.05 m, E_eff = 5e8 Pa, EI = 5,208 N m, rho_s h = 0.05 kg/m^2 (negligible inertia: about 0.6 Pa at omega = 7.85 rad/s and 0.2 m).
- Midpoint pressure balance p_ext - p_f = p_el(u0): bending part 384 EI u0/L^4 = 200 Pa/m x u0; membrane part (parabolic shape, tension T = E_eff h (8/3)(u0/L)^2) 8 T u0/L^2 = 5.33e4 Pa/m^3 x u0^3. Check on the static deflection 0.14363 m: 28.7 + 158 = 186.7 Pa against 200 Pa applied (7 % low from the shape assumption; recomputed 186.75 Pa).
- First trough u0 = 0.2160 m: membrane 537 Pa + bending 43 Pa = 580 Pa (recomputed 537.5 + 43.2), times the shape scaling 1.071 = 621 Pa, so p_f = 200 - 621 = about -420 Pa (range -380 to -420): the fluid pressure is about twice p_ext and 65 % of the net wall load. At late times the fluid pressure is the Poiseuille drop, about 3 % of p_ext.
- Sensitivity to a fractional bias eps of the applied fluid load: the restoring pressure scales as u^2.8 at the trough (membrane dominated), so d ln u / d eps = 0.65/2.8 = 0.23 (trough); added-mass-dominated period: d ln T / d eps = 1/2.

| bias eps of the nodal fluid load | trough error 0.23 eps | period error eps/2 | source |
|---|---|---|---|
| 0.5 % | 0.12 % | 0.25 % | smooth channel (hypothesis) |
| 5 % | 1.2 % | 2.5 % | lower edge of the cavity F1 deficit |
| 13 % | 3.0 % | 6.5 % | upper edge of the cavity F1 deficit |
| 30 % | 6.9 % | 15 % | FSI1-like one-signed worst case |

**Use in this pilot.** The measured steady F1 vs F2 (-6.99 % F_y, -8.63 % G) and, much larger, the D1 spread (+19 to +65 % in loading) fall in or above the last rows of the table; D1's wall lag (12.8 % of the peak at t = 0.425 s) has the sign and size that the law gives for such an amplification. That is a consistency statement. Not shown: that the parity bias is the only cause of the lag (see `LIMITATIONS_AND_UNRESOLVED.md`); that F2 is the physically correct load; the law itself was not re-derived beyond the arithmetic above, and its shape factor 1.071 and sensitivity exponents carry an uncertainty of order 10 % that was not quantified.
