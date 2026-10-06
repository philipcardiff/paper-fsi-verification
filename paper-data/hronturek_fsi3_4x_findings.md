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

## 12. CSM3 structural convergence diagnosis

Structural-only study of the FSI3 plate, run to test whether the solids4foam
solid alone explains the p ≈ 1 displacement convergence of sections 4 and 9.
No fluid, no coupling, no FSI run is involved. Code and compact results:
solids4foam branch `verification/hronturek-csm3` (commit `9a1c5e8ee`,
directory `tutorials/fluidSolidInteraction/HronTurek/verification/csm3`);
machine-readable tables: `hronturek_csm3_levels.csv`,
`hronturek_csm3_observed_order.csv`, `hronturek_csm2_static_levels.csv`.
OpenFOAM v2512 on xenosim, PETSc SNES, the FSI3 solid settings unchanged
(`nonLinearGeometryTotalLagrangianTotalDisplacement`, `leastSquaresS4f`
gradient, `diffStencilLaplacian` stabilisation with `scaleFactor 0.5`).

### CSM3 definition

Turek–Hron CSM tests (Featflow "CSM tests"; structure alone, no fluid):

| | CSM2 (static) | CSM3 (transient) |
|---|---|---|
| geometry | the FSI flag, x ∈ [0.24899, 0.6] m, thickness 0.02 m, 2D plane strain | same |
| support | left face (x = 0.24899 m) clamped, others traction-free | same |
| ρ, ν | 1000 kg/m³, 0.4 | same |
| E (μ) | **5.6e6 Pa (2.0e6)** = the FSI3 plate | 1.4e6 Pa (0.5e6) |
| constitutive law | St. Venant–Kirchhoff, finite strain | same |
| load | gravity (0, −2, 0) m/s² on the plate only | same |
| time | steady | from rest in the undeformed state, `backward` (BDF2), dt = 1e-3 s, 2 s |
| QoIs at A = (0.6, 0.2) | u_x, u_y | mean ± amplitude [frequency] of u_x, u_y |

Important: **CSM3 is not the FSI3 plate.** CSM3 has μ = 0.5e6 (the FSI2
plate stiffness), FSI3 has μ = 2e6, which is the steady CSM2 plate. The
transient CSM3 was run because it is the dynamic structural benchmark that
was asked for; CSM2 was added because it is the same plate as FSI3 and has a
reference that is converged to five digits. The loading is gravity, not the
fluid load, so this isolates the structural operator, not the FSI3 response.
CSM1 (static, μ = 0.5e6) was not run: the single-step Newton solve stalls at
2x and above for the soft plate.

### Reference provenance

Featflow "CSM tests" page (read from the page source, not digitised):

- **CSM2** steady, level 5+1 (22 772 elements): u_x = −0.469000 mm,
  u_y = −16.9739 mm. The Featflow levels 4+2 → 4+3 → 5+1 change by
  < 1e-5 relative, so this value is effectively converged *for Featflow's
  discretisation*.
- **CSM3** transient, level 4+0, dt = 0.005 s: u_x = −14.305 ± 14.305 mm,
  u_y = −63.607 ± 65.160 mm, f = 1.0995 Hz. The spatial levels 2+0 → 4+0
  change by ≤ 0.05%, but the *time step* does not: for dt = 0.02, 0.01, 0.005
  s the u_x mean is −14.404, −14.645, −14.305 mm (non-monotone, 2.4% spread),
  the u_y mean −64.371, −64.766, −63.607 mm (1.8%), the frequency 1.0956,
  1.0978, 1.0995 Hz. The CSM3 reference is therefore a ≈ 1–2% reference, not
  an exact one.

### Mesh family

One block, 105k × 6k cells (the FSI3 solid mesh is k = 1, 2, 4), square
cells, orthogonal and with aspect ratio 1.003; h = 3.343/k mm. Geometry,
material, load, formulation and tolerances are identical on all levels; no
physical parameter is scaled.

| level | cells (length × thickness) | total cells | h (mm) | CSM3 runtime (s) |
|---|---|---:|---:|---:|
| 1x | 105 × 6 | 630 | 3.343 | 52 |
| 2x | 210 × 12 | 2 520 | 1.672 | 156 |
| 4x | 420 × 24 | 10 080 | 0.836 | 628 |
| 8x | 840 × 48 | 40 320 | 0.418 | 3 320 |

16x (161 000 cells) diverged in the single-step Newton solve of the static
case and was not pursued.

### Temporal/iterative control

All times are for the transient CSM3, window = closing 1 s of a 2 s run (one
full period, T ≈ 0.9 s), extrema refined by a parabola.

- **Time step.** 2x with dt = 1e-3, 5e-4, 2.5e-4 s: u_y mean −60.946,
  −60.984, −60.995 mm; u_x mean −12.842, −12.846, −12.847 mm; frequency
  1.1304, 1.1303, 1.1303 Hz. The dt → dt/2 change is 0.06% (u_y mean) and
  the sequence is second-order (ratio 3.5). On 4x, dt 1e-3 → 5e-4 changes
  u_x mean by 0.04% and u_y mean by 0.04%. Against a 1x→2x change of 14–28%
  and a 4x→8x change of 1.3–2.7%, the temporal error is negligible.
- **Iterative error.** SNES `rtol = stol = 1e-6` → 1e-11 on 2x gives
  identical results to every printed digit for both the transient run and
  the static CSM2 (2x, 4x). No linear/nonlinear tolerance effect.
- **Point extraction.** The monitored value is the `pointD` value at the mesh
  vertex (0.6, 0.2), which exists on every level (even thickness cell count).
  Written `pointD` at that vertex and the function-object output agree to
  all printed digits (1x, 2x, 4x). The high-order reconstruction, which uses
  an independent point evaluation, converges to the same limit
  (17.02 mm), so the extraction is not what limits the order.

### Results

CSM3, dt = 1e-3 s. Errors are `(value − Featflow)/|Featflow|`; for the means
a positive error means a smaller deflection than Featflow.

| QoI | 1x | 2x | 4x | 8x | Featflow (dt = 0.005) | err 1x | err 2x | err 4x | err 8x |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| u_x mean (mm) | −9.270 | −12.842 | −14.188 | −14.578 | −14.305 | +35.2% | +10.2% | +0.8% | −1.9% |
| u_x amplitude (mm) | 9.270 | 12.842 | 14.189 | 14.578 | 14.305 | −35.2% | −10.2% | −0.8% | +1.9% |
| u_y mean (mm) | −52.096 | −60.946 | −63.885 | −64.700 | −63.607 | +18.1% | +4.2% | −0.4% | −1.7% |
| u_y amplitude (mm) | 52.127 | 60.907 | 63.895 | 64.732 | 65.160 | −20.0% | −6.5% | −1.9% | −0.7% |
| frequency (Hz) | 1.2251 | 1.1304 | 1.1028 | 1.0954 | 1.0995 | +11.4% | +2.8% | +0.3% | −0.4% |

Successive differences and observed order (mesh ratio 2, p = log₂(d₁/d₂)):

| QoI | d(1x→2x) | d(2x→4x) | d(4x→8x) | ratio, 1x–4x | p, 1x–4x | ratio, 2x–8x | p, 2x–8x |
|---|---:|---:|---:|---:|---:|---:|---:|
| u_x mean (mm) | −3.572 | −1.346 | −0.389 | 2.65 | 1.41 | 3.46 | 1.79 |
| u_x amplitude (mm) | +3.572 | +1.347 | +0.389 | 2.65 | 1.41 | 3.47 | 1.79 |
| u_y mean (mm) | −8.850 | −2.939 | −0.815 | 3.01 | 1.59 | 3.61 | 1.85 |
| u_y amplitude (mm) | +8.780 | +2.988 | +0.836 | 2.94 | 1.56 | 3.57 | 1.84 |
| frequency (Hz) | −0.0947 | −0.0276 | −0.0074 | 3.43 | 1.78 | 3.73 | 1.90 |

Relative last change (4x→8x): u_x 2.7%, u_y mean 1.3%, u_y amplitude 1.3%,
frequency 0.7%.

Steady CSM2 (FSI3 plate), the same four levels, u_y (mm) and its error
against Featflow 5+1 (−16.9739 mm; positive = less deflection):

| variant | 1x | 2x | 4x | 8x | p 1x–4x | p 2x–8x |
|---|---:|---:|---:|---:|---:|---:|
| baseline (as FSI3) | −13.339 (+21.4%) | −15.829 (+6.7%) | −16.694 (+1.65%) | −16.938 (+0.21%) | 1.53 | 1.83 |
| stabilisation 0.25 | −14.209 | −16.175 | −16.797 | −16.966 | 1.66 | 1.88 |
| stabilisation 1.0 | −12.109 | −15.283 | −16.522 | −16.890 | 1.36 | 1.75 |
| high order, p = 2 | −15.884 | −16.840 | −16.993 | −17.018 | 2.65 | 2.56 |
| high order, p = 3 | −16.940 | −17.011 | −17.019 | −17.023 | 3.12 | – |

(u_x gives the same picture: baseline −0.2891, −0.4075, −0.4534, −0.4667 mm
against −0.4690; p = 1.37 and 1.78.) The u_y limit of the high-order
solutions is ≈ 17.02 mm, 0.29% above Featflow's value; the standard scheme
approaches the same limit from below.

### Observed order

- The solids4foam CSM3 displacement QoIs converge **monotonically on all
  four levels**, from the stiff side: the coarse plate is too stiff, the
  amplitudes are too small and the frequency too high.
- The observed order is **not 1**. On the same 1x–2x–4x triple used for FSI3
  it is 1.4 (u_x) and 1.6 (u_y), and on 2x–4x–8x it is 1.8–1.9 for every QoI,
  including the frequency (1.78 → 1.90). It rises monotonically towards 2.
  The steady CSM2 plate repeats this (u_y: 1.53 → 1.83; relative to the
  high-order limit the errors are 21.6%, 7.0%, 1.9%, 0.50%, with local
  orders 1.63, 1.85, 1.97).
- This is *evidence consistent with entry into the asymptotic range*, with
  p still below 2 at the 2x–8x triple. It is not a demonstration of
  asymptotic second order: no 16x level exists, and the high-order
  reference for the limit is itself a solids4foam result. A p = 2
  extrapolation of 4x→8x gives u_x mean −14.71, u_y mean −64.97, u_y
  amplitude 65.01 mm and f = 1.093 Hz; these are indicative only.
- **Against Featflow.** The CSM3 8x values sit 0.7–1.9% from the dt = 0.005
  s reference and the extrapolated values 2–3% on the means, but within the
  reference's own time-step spread: Featflow at dt = 0.01 s gives −14.645,
  −64.766, 64.948 mm, i.e. our extrapolated u_x mean is 0.4%, the u_y mean
  0.3% and the u_y amplitude 0.1% from that column. The CSM3 reference
  cannot discriminate at the 1–2% level, so the CSM3 *level-to-level
  convergence* is the stronger evidence, not agreement with the table.

### Diagnosis

Does the solid alone explain the p ≈ 1 FSI behaviour? **No for the order,
yes in large part for the size and direction.**

1. *The solid is not first order.* On the FSI3 levels the structural order is
   1.4–1.6 and rises to 1.8–1.9 one level later; no first-order plateau is
   seen on either the transient CSM3 or the steady CSM2.
2. *The coarse FSI3 solid levels carry a large pre-asymptotic stiffness
   error*, 21% (1x), 6.7% (2x), 1.7% (4x), 0.2% (8x) of the CSM2 deflection.
   As a proxy, normalising to 4x, the CSM2 u_y is 0.799 : 0.948 : 1 and the
   FSI3 u_y amplitude (section 4) is 0.824 : 0.943 : 1; for u_x the CSM2
   values are 0.638 : 0.899 : 1 against 0.742 : 0.911 : 1 in FSI3. The
   structural stiffness error alone therefore reproduces the sign,
   monotonicity and rough size of the FSI3 1x→2x→4x displacement changes.
   This is a scaling argument, not a decomposition: the FSI3 load is
   fluid-dynamic, its frequency is 5.5 Hz, not 1.1 Hz, and the fluid and
   coupling errors add to the structural one. It also predicts that the FSI3
   displacements still have ≈ 1.5% (u_y) to ≈ 3% (u_x) to gain from the solid
   alone between 4x and 8x.
3. *Where the structural error comes from* (CSM2, baseline discretisation):
   - **Bending / thickness resolution dominates.** With 105 cells along the
     plate, refining only the thickness 6 → 12 → 24 → 48 gives u_y = −13.339,
     −15.088, −15.684, −15.854 mm (errors +21.4, +11.1, +7.6, +6.6%).
     Refining only the length 105 → 210 → 420 → 840 with 6 thickness cells
     gives −13.339, −13.920, −14.075, −14.114 mm and saturates at +16.9%.
   - **The stabilisation is active and not small.** `diffStencilLaplacian`,
     `scaleFactor 0.5`, stiffens the plate. Halving the factor to 0.25 raises
     |u_y| by 6.5% (1x), 2.2% (2x), 0.6% (4x), 0.17% (8x); doubling to 1.0
     lowers it by 9.2%, 3.5%, 1.0%, 0.3%. The effect falls by about 3–4× per
     refinement, consistent with an O(h²) consistent term, not a first-order
     one, and extrapolating the factor to 0 recovers ≈ 1.7 of the 3.7 mm
     deficit at 1x. The unstabilised solid (factor 0) does not converge in
     the linear solve (`DIVERGED_ITS`), so the term cannot simply be removed.
   - **Reconstruction order.** The default linear least-squares displacement
     reconstruction leaves the rest. The existing high-order option
     (`highOrderResidual`, moving least squares) ran cleanly with no changes
     to the code: p = 2 polynomials give 15.88 (1x), 16.84, 16.99, 17.02 mm
     (order 2.6) and p = 3 gives 16.94, 17.01, 17.02, 17.02 mm (order 3.1),
     0.2% from the limit already on the 1x mesh.
4. *Clamp and extraction:* the point extraction is not the cause (above); a
   clamp-localised first-order error is excluded by the observed order on
   the monitored tip displacement itself. Geometric nonlinearity is not the
   cause either: the formulation is the same total-Lagrangian one in all
   runs and the order is unaffected by the stabilisation or reconstruction
   settings. No updated-Lagrangian comparison was made.
5. *Bug or limitation?* No defect was found: this is a discretisation
   limitation (second-order scheme with large pre-asymptotic error on
   6-cell-thick bending) plus a conservative default stabilisation. No code
   change was made. The single non-bug oddity is that the single-step Newton
   solve of the steady problem stalls for E = 1.4e6 at 2x and above; the
   transient solver is unaffected.

### Consequence for FSI3

**CSM3 shows nominal convergence (p = 1.4–1.6 on the FSI3 levels, 1.8–1.9
on 2x–8x), so the FSI3 observed order of ≈ 0.9–1.1 does not follow from the
asymptotic order of the solid.** The solid does contribute a large,
monotone, pre-asymptotic stiffness error on the 1x and 2x FSI3 meshes that
reproduces much of the FSI3 displacement change, and which was misread as
slow convergence. The remaining step from p ≈ 1.5 to p ≈ 1 has to come from
the fluid, the interface or the interaction of the structural error with
them. The statement "the solids4foam solid converges only first order" is
not supported and should not be made.

### Recommended next focused task

One FSI3 diagnostic with the 2x fluid and coupling settings fixed:
**decoupled refinement of the solid only**, i.e. FSI3 on the existing 2x
fluid mesh with the solid mesh at 2x (existing), 4x and 8x (or the 4x solid
with `highOrderResidual`), compared on u_y/u_x amplitude and frequency. If
u_y rises by ≈ 6% and 1.5% as the CSM2 proxy predicts, the solid is the
dominant source of the mesh dependence and the FSI3 level sequence only has
to be re-run with a refined solid; if it does not, the first-order trend lies
in the fluid or the interface.
