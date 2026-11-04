# Overnight-2 preregistration B: Hron-Turek FSI1, CFD1 rigid-flag refinement (decisive gate for the nodal-pressure wall force)

Hashed before any run of this stage. **Selected benchmark: Hron-Turek FSI1** (selection record in `B_SELECTION.md`). This is a NEW exploratory study; it follows the earlier stage-C gate failures (level 1 only: F1 drag +15.73 %, lift +157 %, odd-even a_t 4.1 and 2.6 Pa against the 0.048 Pa gate) and does not reinterpret them.

## Question
At the stage-C resolution (8,330 hex8 cells, factor 1.25) the nodal-pressure wall-integral force F1, which is what `fluid-FSI traction` applies (nodal dilatation interpolated on the face), is +15.7 % in drag against the Featflow converged reference 14.2929 N/m and +157 % in lift against 1.11905 N/m. Is that a resolution effect, as the ring n = 2 nodal force was (-2.1 %, -1.06 %, -0.53 % over three levels), or does the error persist under refinement?

## Configuration (fixed; nothing physical changes)
Stage-C model: FEBio `fluid` module, hex8 slab TZ = 0.015 m (w_z = 0), rho 1000, mu 1 Pa s, K 1e6, rhoi 0, dt 0.5 s, 80 steps to t = 40 s, inlet ramp 2 s, rigid cylinder + flag (CFD1, Re = 20), generator `tools/make_hronturek_cfd_feb.py` unchanged, pristine v4.13. Levels by the generator factor: L1 = 1.25 (existing run C1p, 8,330 cells), **L1.5 = 1.875, L2 = 2.5** (new). Root-junction pins by the same declared rule as C1p: an unpinned run first, pin value = mean of the neighbouring obstacle-wall nodal pressures at the last step (analyser `hronturek_cfd_analyze.py`), then the pinned run. Forces per unit depth (N/m), reconstructions exactly as stage C (F1 nodal-pressure integral plus viscous part; F2 element-mean; F3 momentum balance; definitions unchanged). Reference: Featflow level 6+0, drag 14.2929, lift 1.11905 (fetched 2026-10-31; equality with the Turek-Hron 2006 tables is unconfirmed, as before).

## Predictions and criteria
- **B1 (monotone):** |F1 drag error| decreases L1 > L1.5 > L2 and |F1 lift error| decreases likewise.
- **B2 (rate):** |F1 drag error| at L2 is at most 0.7 times that at L1 (first-order would give 0.5).
- **B3 (parity):** the flag-face odd-even amplitude a_t (definition of stage C) falls by a factor of at least 1.5 from L1 to L2.
- **B4 (soundness each level):** 80/80 steps, no Jacobian/NaN messages, no "Max nr of reformations", |Q_out/Q_in - 1| at most 1e-3, drag spread over the last 20 % at most 1e-3.
Stop rule: L1.5 (pinned) is evaluated first; if |F1 drag error| and a_t both fail to decrease relative to L1, L2 is cancelled and the result is a documented negative. Observed order is stated only with three levels in a consistent asymptotic range (here the level ratios are 1.5 and 1.33, so differences-based order is reported with that caveat). Exploratory; the result alone cannot establish FSI1 verification, and no coupled run is attempted unless B1 to B4 hold and a coupled test is separately registered.
Budget: at most 12 CPU-h, 4 cores per run.
