# IQN-ILS coupling-reuse fix: results with `couplingReuse > 0` may need regenerating

Note for the agents working on the manuscript (2026-10-19). No `.tex` file is
touched.

## Summary

IQN-ILS with `couplingReuse > 0` (reuse of secant modes from previous time
steps) was not robust in solids4foam: the first update of a time step could
move the interface by orders of magnitude too much, and the solid then
diverged (solids4foam issue
[#587](https://github.com/solids4foam/solids4foam/issues/587)). The cause is
found, and the fix is in solids4foam PR
[#598](https://github.com/solids4foam/solids4foam/pull/598) (branch
`fix/issue-587-iqnils-reuse`, open at the time of writing).

**Results computed with `couplingReuse 0` and none of the optional IQN-ILS
switches listed below are unaffected**: with the fix, they are identical to
`development` on every case tested (`TurekHron2006` FSI3, `membraneRoof`,
`beamInCrossFlow`, `collapsibleChannel`).

**Results computed with `couplingReuse > 0` change, and should be regenerated
with the fix** before they are used (see
[Paper results that may need regenerating](#paper-results-that-may-need-regenerating)).

## Cause

Each Dirichlet-Neumann time step starts with a solid predictor: a solid solve
with the previous step's fluid traction. With `predictSolid yes` (the
default), solids4foam applied the reused secant modes in iteration 1 to the
residual of that predictor. That residual does not come from the fluid-solid
coupling map the modes approximate, because the fluid has not yet been solved
at the new time. The least-squares model can then return a correction many
times larger than the true one.

preCICE and CoCoNuT never do this: their first quasi-Newton update of a step
is always applied to the residual of a complete fluid-solid evaluation.
CoCoNuT's `iqni` runs predictor, fluid, solid, then update, and preCICE's
accelerator only receives fresh coupling output.

## The fix (solids4foam PR #598)

- **Main change.** With `skipModesOnSolidPredictor yes` (the new default),
  iteration 1 of a step uses the predictor/relaxation update, and the modes
  are applied from iteration 2, after a complete fluid-solid iteration.
  `reusePreviousStepModes` now acts at that iteration too.
- **New defaults with `couplingReuse > 0`.** `relMinSignificant` defaults to
  `1e-3`, a preCICE-QR2-style relative filter. The modes of a step accepted
  unconverged (`allowUnconvergedCoupling yes`) are discarded.
- **Fixes to optional paths:**
  - the `preciceStyleCouplingQR` solve had the wrong sign;
  - `residualSumPreconditioning` did not scale the residual;
  - `maxReuseUpdateNormRatio` had no effect with `couplingReuse > 0`;
  - the `qrFilterTolerance` test was inconsistent with
    `normalizeCouplingColumns`.

`skipModesOnSolidPredictor no` with `relMinSignificant -1` restores the old
behaviour of the default path.

## Evidence

OpenFOAM-v2412, MeluXina, serial. Iterations are the mean coupling iterations
per time step. "Old" is solids4foam `development` at 76ae126, "new" is
PR #598.

**`TurekHron2006` FSI3** (2x tutorial mesh, Δt = 1e-3 s, coupled from 2 s, run
to 6 s, `nOuterCorr 60`, `allowUnconvergedCoupling yes`; iterations for
t ≥ 2 s):

| Setting | Old | New |
| --- | --- | --- |
| `couplingReuse 0` | 6.9 | 6.9 (identical) |
| `couplingReuse 5` | diverged at 2.91 s | 3.4 |
| `couplingReuse 15` | diverged at 2.01 s | 2.9 |
| `couplingReuse 5`, `relMinSignificant 1.2e-3` | diverged at 4.88 s | – |
| `couplingReuse 15`, `relMinSignificant 1.2e-3` | diverged at 3.18 s | – |
| `couplingReuse 5`, strict test | diverged at 3.23 s | 5.0 |
| `couplingReuse 5` or `15`, `predictSolid no` | completes, 3.7–4.5 | – |

The strict test is `requireAllResidualMeasures yes`; the strict run with
`couplingReuse 0` needs 9.6 iterations per step.

- Every old failure was in **iteration 1** of a step, straight after a
  normally converged step, with the interface update at 0.3–17 m against an
  amplitude of ±30 mm. A relative filter alone delayed but did not prevent it.
- All completed runs agree on the flap-tip displacement over 5–6 s within
  0.4 mm.

**Other cases**, `couplingReuse 5` (reuse 0 results are identical old and
new):

| Case | Reuse 0 | Old, reuse 5 | New, reuse 5 |
| --- | --- | --- | --- |
| `membraneRoof` (to 12 s) | 6.0 | diverged at 1.32 s | 2.8 |
| `beamInCrossFlow` (to 4 s) | 11.5 | diverged at 1.8 s | 10.3 |
| `collapsibleChannel` (to 3.5 s)\* | 14.4 | 12.3 | 12.3 (identical) |

\* `collapsibleChannel` uses `predictSolid no`, so the main change does not
apply to it.

`3dTube/regressionTest.sh` (`couplingReuse 2`) passes with the fix.
`mokLidDrivenCavity` is covered in
[Tutorial defaults changed](#tutorial-defaults-changed).

**preCICE comparison** (`turek-hron-fsi3`, OpenFOAM `pimpleFoam` + the
solids4foam solid of precice/tutorials#915, preCICE 3.4.1, IQN-ILS):

- **Tutorial tolerance (1e-4).** No run diverged: parallel- and
  serial-implicit coupling, `time-windows-reused` 15, 5 and 0, with the QR2
  filter at 1.2e-3 or effectively off.
- **Strict tolerance (1e-6, below the noise floor).** Parallel coupling with
  reuse 15 and **no filter** also diverged, after three windows capped at 100
  iterations.

preCICE's robustness with reuse therefore comes mainly from the consistent
first update. Its relative filter matters once the tolerance is at the noise
floor.

## Paper results that may need regenerating

Only runs with `couplingReuse > 0`, or with one of the optional switches
`reusePreviousStepModes`, `preciceStyleCouplingQR`,
`residualSumPreconditioning`, `maxReuseUpdateNormRatio`,
`normalizeCouplingColumns`, `qrFilterTolerance` or `combinedCouplingSystem`,
are affected. Check the `fsiProperties` of each run. Among the tutorials,
only `3dTube` (`couplingReuse 2`) and `hessenthalerFsi` (`couplingReuse 5`)
set `couplingReuse > 0`.

- **Hessenthaler** (`sec:hessenthaler`; the drift paragraph in
  `sectionCrossBenchmark.tex`). The campaign used the tutorial's
  `couplingReuse 5`. The default-test drift study (40 000 steps to 80 s, about
  three iterations per step) and the strict-test restarts (10–11 iterations
  per step) were computed with the old reuse algorithm. Their iteration
  counts, and possibly the accumulated-drift values, change with the fix.
  The strict-test runs from t = 0 that diverged with reuse 5 (F1S2 at 3.57 s
  and 5.45 s, F2S2 at 1.08 s; issue #587) show the defect fixed here, not a
  property of the case. Strict runs with `couplingReuse 0` are unaffected.
  If reuse is kept, regenerate with the fix; the campaign can then use reuse
  with the strict test, which previously diverged.
- **3dTube** (`sec:3dtube`; `paper-data/3dtube_level3_findings.md`). If the
  IQN-ILS runs used the tutorial's `couplingReuse 2`, the IQN-ILS solutions
  (levels 1–3) and iteration counts (about 15.5 per step), and the Robin vs
  IQN-ILS comparison, change at the coupling-tolerance level. The
  Robin-Neumann results are unaffected.
- **beamInCrossFlow IQN-ILS option sweep** (solids4foam PR #518, quoted in
  `benchmark_inventory.md`). If the sweep used `reusePreviousStepModes`,
  `maxReuseUpdateNormRatio`, `preciceStyleCouplingQR`,
  `residualSumPreconditioning` or `normalizeCouplingColumns`, those results
  were obtained with the defects listed above (the `preciceStyleCouplingQR`
  solve had the wrong sign).
- **Any IQN-ILS iteration count or robustness statement** for a case run with
  `couplingReuse > 0`. Once the fix is merged, `couplingReuse 0` is no longer
  needed for robustness, if the text recommends it anywhere.

Unaffected, because they were run with `couplingReuse 0`: ringAddedMass,
womersleyTube, collapsibleChannel, Hron–Turek FSI1/FSI2/FSI3 (the FSI3 4x
study used "no reuse": `hronturek_fsi3_4x_findings.md`), Mok, and the
Robin-Neumann results throughout. **However, several of these tutorials now
default to `couplingReuse 5`**, so a rerun from the current tutorial settings
no longer reproduces them; see the next section.

## Tutorial defaults changed

PR #598 also sets `couplingReuse 5` in the IQN-ILS set-ups of six tutorials,
four of which are paper cases:

- `ringAddedMass`;
- `womersleyTube`;
- `TurekHron2006` (Hron–Turek), which covers FSI1, FSI2 and FSI3, because
  `Allrun fsi1|fsi2` edit the same `fsiProperties.iqnils`;
- `mokLidDrivenCavity`;
- `membraneRoof` and `blobInTreacle`, which are not paper cases.

Full tutorial runs (OpenFOAM-v2512, Linux, serial), mean coupling iterations
per time step:

| Tutorial | `couplingReuse 0` | `couplingReuse 5` |
| --- | --- | --- |
| `TurekHron2006` (FSI3, 6 s, incl. the uncoupled 0–2 s) | 5.0 | 2.6 |
| `TurekHron2006 fsi1` (30 s; t ≥ 2 s) | 3.7 | 2.4 |
| `TurekHron2006 fsi2` (10.5 s; t ≥ 2 s) | 4.6 | 2.6 |
| `ringAddedMass` | 4.0 | 2.4 |
| `womersleyTube` | 14.6 | 8.4 |
| `mokLidDrivenCavity` | 4.3 | 2.6 |
| `membraneRoof` | 5.9 | 2.7 |
| `blobInTreacle` | 6.2 | 4.2 |

All six regression tests pass unchanged with reuse, so the QoIs move only at
the coupling-tolerance level: for example, the final FSI1 tip displacement
changes by 0.01%, and the FSI2 `u_y` range over 9.5–10.5 s by about 0.1%. But:

- **Results are no longer reproduced from the tutorial settings.** A paper
  result from these tutorials, computed with `couplingReuse 0`, is reproduced
  only with `couplingReuse 0`. Either state the coupling settings with each
  result, or regenerate with the new defaults.
- **Iteration counts change.** Every IQN-ILS iteration count in the paper for
  these cases (for example ringAddedMass 2.99/3.95/4.91, womersleyTube 14.6 or
  14.8, FSI3 about 7–9, and the IQN-ILS vs Robin-Neumann cost comparisons)
  corresponds to `couplingReuse 0`. With reuse, IQN-ILS needs 35–50% fewer
  iterations, which changes the Robin-Neumann/IQN-ILS cost comparisons. For
  example, in womersleyTube IQN-ILS with reuse took 8.4 iterations per step
  in this run, against 9.3 quoted in the paper for Robin-Neumann and 14.6 for
  IQN-ILS without reuse.
- **The `verification/` READMEs and `Allverify` studies** of these tutorials
  were computed without reuse and have not been regenerated. `Allverify`
  `--reuse` treats changed tutorial files as new settings, so a rerun
  recomputes them with reuse.
- **ringAddedMass amplitude loss.** The unexplained mesh-dependent IQN-ILS
  amplitude loss (`benchmark_inventory.md`) was observed without reuse, and
  has not been re-examined with it.

`beamInCrossFlow`, `collapsibleChannel` (paper case), `fillingElasticContainer`
and `perpendicularFlap` keep `couplingReuse 0`: reuse gave about 9% fewer
iterations or none, and `collapsibleChannel`'s worst step went from 36 to 59
iterations.

## Manuscript text to revisit

- `sectionNumericalModel.tex`, Dirichlet–Neumann paragraph and
  Algorithm `alg:partitioned_fsi`. With the fix, when secant modes from
  previous steps are reused, the predictor iteration is not accelerated, and
  the first IQN-ILS update of a step is applied after a complete fluid–solid
  iteration. The filter sentence could state that the relative QR2-type
  filter is on by default with reuse.
- The Hessenthaler and 3dTube items above, once regenerated, with the
  solids4foam commit that contains PR #598.

## Not included in the fix

Keeping the secant of the final iteration of each step, as preCICE and
CoCoNuT do, reduced the FSI3 iterations by about 17%. However, it made
`collapsibleChannel` with `couplingReuse 5` need two to three times more
iterations per step, and it then diverged at 0.15 s. It was left out.

## Runs

MeluXina, `/project/home/p201287/philip/precice-issue587`:

- `study/`: preCICE variants;
- `native/`: old-code and prototype FSI3 variants;
- `prtests/`, `prtests2/`: old vs new on the cheap cases;
- `status.sh`: summary of the above.
