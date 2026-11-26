# Ring added-mass problem with the full interface pressure condition

Source of Table `tab:ring`, the bullets of Section `sec:ring`, the ring rows of
Table `tab:refvals_ring` and `../reference_values/reference_values.csv`, panel (b)
of Figure `fig:ring_febio_ratio`, and the IQN-ILS column of the ring tolerance
sweep (`../tolerance_sweeps/ring_strong_m2_n100_strict_iqn_mwp.csv`).

## Provenance

- **Build:** solids4foam `development` 7a2f0a63e with the interface-pressure
  changes of branch `fix/movingWallPressure-721`; OpenFOAM v2512, GCC, UCD Sonic
  (Linux).
- **Study:** the tutorial's full `./Allverify --cores 16` study (35 serial runs),
  unchanged except for the interface pressure of the IQN-ILS set-up:
  `0/fluid/p.dirichletNeumann` uses `movingWallPressure` (wall acceleration,
  viscous and body-force terms, `viscousTerm yes`; its mean is removed in this
  enclosed fluid) instead of `zeroGradient`. The Robin-Neumann runs use
  `elasticWallPressure` with the viscous term (`viscousTerm yes`, the default).
- **Earlier values:** the zero-gradient results (solids4foam
  `ringAddedMass/verification/README.md`, v2412, xenosim) are listed next to the
  new ones in `tables.md`.

## Files

| File | Content |
|---|---|
| `ringAddedMass_results.csv` | one row per run, written by the tutorial's verification script |
| `verification_summary.md` | the tutorial's summary table and acceptance checks (all pass) |
| `tables.md` | the paper tables, new against old values, orders, coupling comparison, amplitude loss and decay rates, and the tolerance sweep |

## Main differences from the zero-gradient interface

- Strong level, mesh 1: amplitude loss 2.10% to -0.001% per period; frequency
  -0.13%; ratio error +0.068% to -0.061%.
- Ratio orders (mesh 1, 2, 4): 1.90, 1.84, 1.90 (weak, moderate, strong); were
  1.77, 1.73, 1.68. Time-step orders 1.89-1.90.
- IQN-ILS and Robin-Neumann agree to 4e-4% in frequency (mesh 2).
