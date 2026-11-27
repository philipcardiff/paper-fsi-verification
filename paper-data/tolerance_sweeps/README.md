# Coupling-tolerance sweeps with the strict test

Iterative-error evidence for Section `sec:iterative_errors` (Table `tab:tol_sweeps`).
Strict test (`requireAllResidualMeasures yes`), IQN-ILS `couplingReuse 0`.
solids4foam `a421f063d`, OpenFOAM-v2512, MeluXina (jobs 5410800 ring, 5410835 collapsible), 2026-11-01.

- `ring_strong_m2_n100_strict.csv`: ringAddedMass, strong level, mesh 2, 100 steps
  per period, standard solid; IQN-ILS and Robin–Neumann at tolerances 1e-4 to 1e-8.
  The ratio error is the wet/dry frequency ratio against the exact ratio, with the
  dry ring at the same discretisation (`ts-ram-dry-m2-n100`). Driver: `ring_tol.py`
  (wraps the campaign driver `ma_run.py`; the Robin–Neumann pressure and flux
  tolerances are 10x and 5000x the displacement tolerance, the tutorial ratios).
- `ring_strong_m2_n100_strict_iqn_mwp.csv`: the IQN-ILS rows of Table `tab:tol_sweeps`,
  rerun with the full interface pressure condition (`movingWallPressure` with
  the viscous term and its mean removed, instead of `zeroGradient`); solids4foam
  7a2f0a63e with branch `fix/movingWallPressure-721`, OpenFOAM v2512, UCD Sonic.
  Dry ring at the same discretisation: ω = 2.7125759973 (Sonic run). Written by
  `ring_tol_mwp.py` from the run directories. The IQN-ILS rows of
  `ring_strong_m2_n100_strict.csv` (zero-gradient interface) are superseded; its
  Robin-Neumann rows (acceleration-only Robin condition, MeluXina dry ring) are
  superseded by `ring_strong_m2_n100_strict_rn_full.csv`.
- `ring_strong_m2_n100_strict_rn_full.csv`: the Robin-Neumann rows of Table
  `tab:tol_sweeps` (27 Nov 2026), rerun with the full interface pressure condition
  (`elasticWallPressure`, `viscousTerm` default `yes`); solids4foam PR #777 head
  b9b5e2b9e, OpenFOAM v2512, GCC, UCD Sonic (job 1192806). Same settings as above
  (Robin pressure and flux tolerances 10x and 5000x the displacement tolerance);
  ratio error against the same Sonic dry ring as the IQN-ILS column
  (ω = 2.7125759973). Written by `ring_tol_rn_full.py` (`ring_tol_mwp.py` with RN
  labels); job scripts `ringtol_rn.slurm`. `ab_viscousTerm_no_tol1e-6.txt`: A/B at
  tolerance 1e-6 with `viscousTerm no` (`ringtol_rn_novisc.slurm`), ratio error
  -0.008497% against -0.008491% with the full condition. Run directories on Sonic
  `/scratch/pcardiff/mwp721/ring/cases/tol_rn_*` (not copied).
- `collapsible_F3_dt160_strict.txt`: collapsibleChannel, fluid level 3, solid
  160x8 high order, dt = 1/160 s; Robin–Neumann and IQN-ILS at 1e-3 to 1e-6
  (driver summary lines). Driver: `coll_tol.py` (wraps the tutorial's
  `collapsible_channel_verification.py`; Robin–Neumann pressure tolerance equals
  the displacement tolerance, as in the tutorial). IQN-ILS failed at every
  tolerance (1e-5, 1e-6: stalled at nOuterCorr = 200 early in the run; 1e-3, 1e-4:
  coupling divergence after the external-pressure ramp, ending in the
  interface-to-interface face-mismatch check); Robin–Neumann 1e-6 stalled at
  nOuterCorr in the first step.
