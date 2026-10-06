# Hron–Turek FSI3 4x verification findings

Evidence report for the Hron–Turek FSI3 benchmark (`tab:hronturek`,
`tab:hronturek_ref` in `sectionBenchmarkSuite.tex`) before any manuscript
change. It closes item 2 of `manuscript_todos.md` (FSI3 4x level) and item
S3 (independent Δt study on 2x), and contributes to `reference-provenance`
in `benchmarks.yaml`. No `.tex` file is touched.

The per-quantity numbers are in `hronturek_fsi3_mesh_study.csv`. That file
is generated from `fsi3_iqnils_refinement_study.json` in the solids4foam
branch below.

## 1. Executive scientific conclusion

The third level does **not** confirm that FSI3 converges onto the Featflow
level-4 reference. All three levels completed, and every coupled step met
the interface tolerance.

- **Displacements.** They change monotonically but are not asymptotic. From
  2x to 4x the `u_y` amplitude changes by 6.0% and the `u_x` mean and
  amplitude by about 10%. The observed order along the scaled refinement
  path is about 0.9–1.1, half the formal order. The 2x agreement with
  Featflow (2–3%) was a crossing: at 4x the displacements are 3.5% (`u_y`)
  and about 7% (`u_x`) *past* Featflow, and still moving away.
- **Frequencies.** They converge monotonically onto Featflow (+0.25% at 4x).
  Their observed order (about 0.4) is not interpretable, because the
  time-step and coupling contributions (up to 0.4%) are as large as the
  2x→4x change.
- **Drag mean.** Within 0.3% of Featflow at 4x. The 2x→4x change (0.5%) is
  inside its late-time variability (0.8%).
- **Lift amplitude.** The previous 13% gap does not close onto Featflow. On
  the benchmark definition (extrema of the last period) it is +15.0% at 4x,
  but at Δt = 0.00025 s that value is inflated by coupling noise (Section
  7). With the noise removed (4 ms moving average) the lift amplitude falls
  228.8 → 172.3 → 165.6 N/m: +11.9% at 2x and +7.6% at 4x. The sequence
  decelerates strongly (ratio 0.12) towards about 164–166 N/m, roughly
  6–8% above Featflow. Its apparent order of 3.1 exceeds the formal order,
  so an asymptotic range is not demonstrated.
- **Drag amplitude.** Still rising: +10.4% at 4x smoothed, +18.1% raw
  (noise-inflated). Not converged.
- **Observed order.** No quantity shows an order consistent with the formal
  second order (taken as 1.5 ≤ p ≤ 2.5). No Richardson extrapolation is
  meaningful.
- **Reference uncertainty.** The "≈1%" uncertainty assigned to the FSI3
  reference is defensible only for the drag mean and the frequencies. It is
  about 2% for the `u_y` amplitude, and 3–5% or more for the `u_x` mean and
  amplitude and the drag and lift amplitudes (Section 3).
- **Category.** FSI3 remains category C. The completed study sharpens the
  wording but does not justify promoting it.

## 2. Provenance

- **solids4foam.** Branch `verification/hronturek-fsi3-4x`, commit
  `d49552f4c`, from `development` at `aee0c8e35`. `src/` and the tutorial
  inputs are unchanged. The commit changes only the verification driver,
  the README, the reference JSON (added tables, primary values unchanged),
  a new analysis script and compact results in
  `tutorials/fluidSolidInteraction/HronTurek/verification/reference/fsi3_refinement/`.
  The histories of every run from t = 4.5 s at 1 ms are committed, so all
  numbers can be re-extracted without the run directories.
- **Platform.** OpenFOAM v2512 on xenosim, a shared 192-core Linux node,
  with a private solids4foam build. A 1x replicate with OpenFOAM v2412 on
  MeluXina agrees within 1.6% for all primary QoIs. The v2512 1x/2x runs
  reproduce the earlier v2412/M1 record (README) to within 0.5% on 2x and
  2.4% on 1x.
- **Set-up, unchanged from the 1x/2x study.** IQN-ILS with `predictor yes`,
  `outerCorrTolerance 1e-5`, `nOuterCorr 30`, coupling from t = 2 s, BDF2
  (`backward`) in fluid and solid, `velocityLaplacian` mesh motion with
  quadratic inverse-distance diffusivity, and fluid solver tolerances 1e-6.
  Each run goes to t = 7 s and is analysed over the closing 1 s. The mean
  and amplitude are `(max ± min)/2` over the last full `u_y` period, and
  the frequency is the average over all full periods in the window.
- **Refinement.** Ratio 2 in-plane (blockMesh divisions ×2 with gradings
  kept), one spanwise cell, and Δt halved with h.

| Level | Fluid cells | Solid cells | Δt (s) | Ranks | Wall time | Mean (max) FSI iterations | Largest final residual | Status |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1x | 5 336 | 630 | 0.001 | 1 | 1.4 h | 6.96 (11) | 9.994e-6 | completed |
| 2x | 21 344 | 2 520 | 0.0005 | 8 | 5.5 h | 7.01 (13) | 9.999e-6 | completed |
| 4x | 85 376 | 10 080 | 0.00025 | 32 | 40.8 h | 6.69 (15) | 9.9996e-6 | completed |

The 4x level used 32 ranks rather than the nominal 16. That gives the same
cells per rank as 2x on 8 ranks. A strong-scaling test on one 128-core
node gave 9.1 s/step on 32 ranks, 11.1 on 64 and 18.9 on 128, so more
ranks do not help.

## 3. Existing benchmark and reference status

All values below come from the Featflow FSI tests page (retrieved
2026-10-03). All nine FSI3 table rows and the 2006 tables are now recorded
in the solids4foam reference JSON.

- **A. Featflow table, level 4+0 (15 872 elements), Δt = 0.00025 s.** This
  is the verification reference:
  `u_x = −2.88 ± 2.72 mm [10.93]`, `u_y = 1.47 ± 34.99 mm [5.46]`,
  `drag = 460.5 ± 27.74 [10.93]`, `lift = 2.50 ± 153.91 [5.46]`. It is the
  discretisation of the published `ref_fsi3.point` history
  (`media/fsi/data/fsi3/0p00025`). The driver reproduces it from that
  history to ≤0.6% (drag amplitude 0.8%) and frequencies ≤0.3%.
- **B. Turek & Hron (2006) summary values.** Historical only. They are level
  4 of the 2006 proceedings table at the smaller of its two time steps
  (Δt = 0.0005 s). On the Featflow page these tables survive only inside an
  HTML comment marked "old values", captioned Δt = 0.001, 0.0005. They
  differ from A by 18% in drag amplitude and give the frequency to two
  digits.
- **C. Tuković et al. (2018).** Historical only:
  `−2.72 ± 2.58 [11.07]`, `1.67 ± 33.84 [5.53]`, `459.18 ± 24.86`,
  `1.59 ± 155.9`.

Uncertainty of A, from the Featflow tables themselves:

| QoI | L2 | L3 | L4 | L3→L4 change | Monotone L2→L4? |
|---|---:|---:|---:|---:|---|
| `u_x` mean (mm) | −3.02 | −2.77 | −2.88 | 3.8% | no |
| `u_x` amplitude (mm) | 2.85 | 2.61 | 2.72 | 4.0% | no |
| `u_y` amplitude (mm) | 35.73 | 34.43 | 34.99 | 1.6% | no |
| frequency (Hz) | 5.36 | 5.46 | 5.46 | 0 | – |
| drag mean (N/m) | 458.7 | 459.1 | 460.5 | 0.3% | yes |
| drag amplitude (N/m) | 28.80 | 26.50 | 27.74 | 4.5% | no |
| lift amplitude (N/m) | 146.00 | 149.91 | 153.91 | 2.6% | yes, but +4 N/m per level: not converging |

At level 4 the three tabulated Δt values differ by up to 0.7% in lift
amplitude and 1.0% in drag amplitude. The reference history is not strictly
periodic either: over its 7.9 periods (t = 5–6.44 s) the per-period lift
amplitude drifts from 157.7 to 153.5 N/m (2.7%).

## 4. Spatial / scaled-refinement results

Errors are signed, `(value − ref)/|ref|`. For the `u_x` mean, which is
negative, a negative error therefore means a larger magnitude. The late-time
variability is the range of the per-period values over the 13 periods from
t = 4.5 s to the end of each run. It is the threshold for calling a 2x→4x
change resolved. After saturation the means, amplitudes and frequency still
wander by up to about 1% on time scales longer than the 1 s window, which
the standard two-against-two periodicity test (passed by all runs, worst
1.9%) does not see.

| QoI | 1x | 2x | 4x | Featflow L4 | Error 2x | Error 4x | Ratio | Order | Verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `u_x` mean (mm) | −2.232 | −2.788 | −3.085 | −2.880 | +3.2% | −7.1% | 0.53 | 0.91 | monotone, not asymptotic |
| `u_x` amplitude (mm) | 2.163 | 2.657 | 2.917 | 2.720 | −2.3% | +7.2% | 0.53 | 0.93 | monotone, not asymptotic |
| `u_y` amplitude (mm) | 29.85 | 34.18 | 36.23 | 34.99 | −2.3% | +3.5% | 0.48 | 1.07 | monotone, not asymptotic |
| `u_y` frequency (Hz) | 5.588 | 5.523 | 5.474 | 5.460 | +1.2% | +0.25% | 0.75 | 0.42 | monotone, not asymptotic; Δt/coupling ≈ change |
| drag mean (N/m) | 456.6 | 459.7 | 462.1 | 460.5 | −0.2% | +0.3% | 0.76 | – | 2x→4x within variability |
| drag amplitude (N/m) | 22.87 | 28.03 | 32.76 | 27.74 | +1.0% | +18.1% | 0.92 | 0.13 | noise-inflated |
| drag amplitude, smoothed | 22.63 | 27.63 | 30.63 | 27.74 | −0.4% | +10.4% | 0.60 | 0.74 | monotone, not asymptotic |
| lift amplitude (N/m) | 229.1 | 174.1 | 177.0 | 153.9 | +13.2% | +15.0% | −0.05 | – | sign change within variability; noise-inflated |
| lift amplitude, smoothed | 228.8 | 172.3 | 165.6 | 153.9 | +11.9% | +7.6% | 0.12 | 3.1 | monotone, order > formal: not demonstrably asymptotic |
| `u_y`, lift means | | | | | | | | – | near zero: undefined |

Late-time variability on 4x: `u_x` mean 1.0%, `u_x` amplitude 0.8%, `u_y`
amplitude 0.6%, frequency 0.1%, drag mean 0.8%, smoothed lift amplitude
1.0%, smoothed drag amplitude 4.2%.

The 2x→4x changes in displacement (6–11%) are resolved: they are about 10×
the variability. The three-level sequence shows that the two-level
"apparent orders" of 2.6–3.8 in `tab:hronturek` were artefacts of the
crossing near the Featflow value.

## 5. Temporal-error assessment

Δt is halved with h, so the observed order is that of the combined
space–time path. Four runs on the 2x mesh separate the contributions there:
two time steps (Δt = 0.0005 and 0.00025 s), each at two interface
tolerances (1e-5 and 1e-6). The changes below are relative, computed on the
smoothed statistics, and for means are changes in magnitude.

| QoI | Halve Δt (1e-5) | Halve Δt (1e-6) | 1e-5→1e-6 (Δt 0.0005) | 1e-5→1e-6 (Δt 0.00025) | 2x→4x |
|---|---:|---:|---:|---:|---:|
| `u_x` mean | +1.4% | −0.7% | +1.7% | −0.4% | +10.6% |
| `u_x` amplitude | +1.4% | −0.4% | +1.7% | −0.2% | +9.8% |
| `u_y` amplitude | +0.8% | −0.0% | +0.2% | −0.7% | +6.0% |
| frequency | −0.1% | +0.4% | −0.4% | +0.0% | −0.9% |
| drag mean | +0.1% | −0.2% | +0.2% | −0.1% | +0.7% |
| drag amplitude | +3.3% | +2.5% | −1.4% | −2.2% | +10.9% |
| lift amplitude | +0.3% | +1.2% | −0.5% | +0.4% | −3.9% |

The run at Δt = 0.00025 s with tolerance 1e-6 stopped at t = 6.73 s (see
Section 7). It is analysed over 5.73–6.73 s, and its 18 923 completed steps
all converged.

- **Displacements.** The temporal and coupling contributions (≤1.7%) are a
  minor fraction of the 2x→4x change (6–11%), so the path is
  mesh-dominated.
- **Frequency.** Δt/coupling effects (≤0.4%) are as large as the 2x→4x
  change (0.9%), so its order is not interpretable as spatial.
- **Smoothed lift amplitude.** Up to about 30% of the 2x→4x change may be
  temporal or coupling. Smoothed drag amplitude: about 25%.

## 6. Comparison with other references (context only)

- **Tuković et al. (2018)** reported a lift amplitude of 155.9 N/m.
- **Featflow's own lift amplitude** rises monotonically with level (146.0,
  149.9, 153.9) without converging.
- **This sequence** converges *downwards* from above, towards about 165 N/m.

The exact solution is therefore bracketed only loosely, somewhere between
about 154 and 165 N/m. No literature survey of other codes was done here.

## 7. Coupling/iterative error

- **Tolerance met throughout.** Every coupled step on every level, including
  every step of the analysis window, ended with the interface residual below
  1e-5. An unconverged step is fatal, and the residual files confirm this.
  Iterations do not deteriorate with refinement (about 7 per step on all
  levels).
- **Coupling noise in the forces at the 4x time step.** At Δt = 0.00025 s
  the 1e-5 tolerance leaves step-to-step force noise:
  - lift: 7.3 N/m rms (peaks ≈ 35 N/m) on 4x;
  - drag: 1.1 N/m rms on 4x;
  - for comparison, 1.7 and 0.3 N/m on 2x at Δt = 0.0005 s.

  The same noise appears on 2x when only Δt is halved, and a 1e-6 tolerance
  reduces it 8×. It is coupling error that grows as Δt falls. It inflates
  the extrema-based amplitudes: on 2x at Δt = 0.00025 s, tightening to
  1e-6 changes the benchmark drag amplitude by −9.8% and the lift
  amplitude by −3.2%, but the smoothed values only by −2.2% and +0.4%.
  Displacements are unaffected. This is the same mechanism as the FSI2
  tolerance finding in `sec:iterative_errors`, now shown for FSI3 at the
  finest time step.
- **Consequence for acceptance.** The driver's (unchanged) 15% tolerance
  fails on the 4x drag amplitude (18.1%).
- **New limitation: IQN-ILS does not reach 1e-6 at Δt = 0.00025 s.** With
  `outerCorrTolerance 1e-6`, IQN-ILS twice stalled at 5–9 × 10⁻⁶ (at
  `nOuterCorr 30`) and then diverged within a step, after which the solid
  SNES failed:
  - on 2x at Δt = 0.00025 s, at t = 6.73 s;
  - on 4x, restarted from the 1e-5 run at t = 4 s, at t = 4.73 s.

  The fluid solver tolerances (1e-6, absolute) appear to set a floor on the
  relative interface residual, as reported for FSI1. A clean 1e-6 4x run
  would need tighter fluid tolerances. The 4x coupling error is therefore
  bounded from the 2x pair at the same Δt, as given above.

## 8. Remaining uncertainties

- **Offset or reference error?** Whether solids4foam converges to values
  different from the exact solution (displacements and lift amplitude
  ≈ 4–8% from Featflow), or Featflow level 4 is itself several percent off,
  cannot be decided from three levels on a path that is not asymptotic.
- **Possible cause of the displacement overshoot.** The plate has 6, 12 and
  24 cells through its thickness. The rising displacements are consistent
  with a cell-centred FV plate in bending that is still too stiff at 2x.
  This is not demonstrated.
- **Frequency.** Resolved only to about 0.4%.
- **Late-time wander.** The ≈1% slow wander limits single-window
  statistics. A longer window, or several windows, would sharpen the
  smaller differences.

## 9. Consequences for the manuscript

- **`tab:hronturek`.** Replace the two-level FSI3 apparent orders
  (1.1–3.8). Add the 4x column and three-level ratios and orders, and add
  the smoothed lift and drag amplitudes alongside the benchmark-defined
  values.
- **Text following the table** ("FSI2/FSI3 displacements and drag … within
  1.2–3.1% … no order can be claimed"):
  - FSI3 displacements are *not* converging onto Featflow: the 4x mesh
    overshoots it.
  - The observed order is about 1, not "high".
- **Lift-amplitude statement** ("falls from 48.9% to 13.3%"): add 4x
  (+15.0% raw, +7.6% smoothed), and say that the gap narrows but does not
  close.
- **Reference uncertainty** ("references uncertain by ≈1%", in the
  reference-category paragraph and the closing paragraph of the Hron–Turek
  subsection): for FSI3 this holds only for drag mean and frequency; it is
  3–5% for the amplitudes and `u_x`.
- **Reference provenance** (`manuscript_todos.md` items on `tab:hronturek_ref`
  and the citable source): the 2006 values are the 2006 proceedings table
  (level 4, Δt = 0.0005 s), retained by Featflow only as commented-out "old
  values"; table A is the later Featflow table.
- **`sec:iterative_errors`.** Add the FSI3 evidence: force noise at 1e-5
  grows as Δt falls and inflates force amplitudes; IQN-ILS has a residual
  floor near 5e-6 at Δt = 0.00025 s.
- **Remove the `\todo{ESSENTIAL: run FSI3 at 4x …}`** and the S3 Δt-study
  item; both are done.

## 10. Recommended remaining work

1. **CSM3 and CFD3 sub-benchmarks** (solid-only and rigid-plate) at the
   1x/2x/4x resolutions, to tell whether the displacement overshoot comes
   from the plate or the fluid discretisation. This is the most informative
   next step. Neither exists in solids4foam yet.
2. **Force amplitudes at fine Δt.** Either tighten the fluid solver
   tolerances so that a 1e-6 coupling tolerance is reachable, or adopt and
   justify a noise-robust amplitude measure.
3. **An 8x level**, for a genuine asymptotic order. Expensive: the case does
   not strong-scale beyond 32 ranks, so it would be well over a week.
4. **A longer analysis window, or several windows**, to beat down the ≈1%
   late-time wander.

## 11. Suggested manuscript wording

> For FSI3, three meshes (5 336, 21 344 and 85 376 fluid cells, with Δt
> halved with h) give monotone sequences for every primary quantity except
> the extrema-based lift amplitude. None is in the asymptotic range: the
> observed order along this space–time path is about 1 for the
> displacements, half the formal order. The frequency converges to within
> 0.3% of the Featflow level-4 value and the mean drag to 0.3%. The
> displacements, which agreed with Featflow to 2–3% on the intermediate
> mesh, overshoot it on the finest mesh (`u_y` amplitude +3.5%, `u_x`
> about 7%). The lift amplitude, once the coupling noise at the finest time
> step is filtered, falls from 11.9% to 7.6% above Featflow and appears to
> approach a value about 7% above it. These differences are comparable to
> the reference's own level-3-to-4 changes (2.6% in lift amplitude, ≈4% in
> `u_x` and drag amplitude), so agreement with the reference cannot be
> asserted to better than about 5% for amplitude quantities, and FSI3
> remains a category C benchmark.
