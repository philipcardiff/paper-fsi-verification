# cavityFlexibleBottom quasi-monolithic, solids4foam #618 build

- **Runs:** solids4foam `1252bd2e5` (with #618, the transposed viscous term
  that supersedes #611). Run by the solver-paper campaign (`rf-cav-qm-L1`–`L4`
  in `MeluXina:/project/home/p201136/philip-ma/runs-rfix/`). Steady state at
  t = 300 s. Source: paper-monolithic-fsi `evidence/m-runs-a-rfix`,
  `paper-data/campaign/m-runs-a-rfix/RFIX_findings.md`.
- **Sampling:** those runs monitor the plate underside. The top-surface values
  (point A, x = 4, plate top) were read from `300/solid/pointD` at the nearest
  solid mesh point, on 27 Oct 2026. The same read at the underside reproduces
  the runs' monitor output to all six digits.
- **Forces:** total `F_y` from `postProcessing/fluid/forces/0/force.dat` at
  t = 300 s.
- **Previous values:** the `_597` columns are the #597-build values in
  solids4foam `f356c3242`,
  `cavityFlexibleBottom/verification/results/four_level_partitioned_monolithic.csv`.
  The largest change is 1.3e-4 relative in u_y and 1.8e-4 in `F_y`.
