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
  orders 1.63, 1.86, 1.95).
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
u_y rises by ≈ 5.5% from the 2x to the 4x solid (and a further ≈ 1.5% to
8x), as the CSM2 proxy predicts, the solid is the dominant source of the mesh dependence and the FSI3 level sequence only has
to be re-run with a refined solid; if it does not, the first-order trend lies
in the fluid or the interface.

## 13. FSI3 solid-only refinement at fixed fluid resolution

Error-decomposition study asked for by the section 12 recommendation. Code
and compact results: solids4foam branch
`verification/hronturek-fsi3-solid-refinement` (commit `7c15e5bcd`, driver
option `--solid-refinement`, `scripts/fsi3_solid_refinement_analysis.py`,
`reference/fsi3_solid_refinement/`). Machine-readable:
`hronturek_fsi3_solid_refinement.csv`, `.json`. OpenFOAM v2512 on xenosim.
No `.tex` change. **This is not a spatial-order study (the fluid is fixed);
no order is derived.**

### Setup

The 2x fluid mesh (21 344 cells, 168 interface faces), `dt = 0.0005 s`,
IQN-ILS with `outerCorrTolerance 1e-5`, mesh motion, materials, end time
7 s and the 1 s analysis window are **identical in all runs and identical to
the existing 2x run**, which is reused as the solid-2x baseline (its case
directory differs from the new ones only in `blockMeshDict` of the solid and
the decomposition: 8 ranks instead of 14). Only the solid block mesh
changes: 2x = 210 × 12, 4x = 420 × 24, 8x = 840 × 48 cells (2 520 / 10 080 /
40 320). Solid formulation, `diffStencilLaplacian` stabilisation
(`scaleFactor 0.5`) and reconstruction order are unchanged. QoIs use the
same extraction as the earlier tables (last full u_y period of the closing
1 s). Forces: extrema-based amplitude and, in parentheses, the 4 ms
moving-average amplitude.

### Results

| Quantity | fluid 2x / solid 2x | solid 4x | solid 8x | Δ 2→4 | Δ 4→8 | Δ 2→8 |
|---|---:|---:|---:|---:|---:|---:|
| u_x mean (mm) | −2.7878 | −2.7970 | −2.7705 | +0.33% | −0.95% | −0.62% |
| u_x amplitude (mm) | 2.6568 | 2.6544 | 2.6292 | −0.09% | −0.95% | −1.04% |
| u_x frequency (Hz) | 11.046 | 11.002 | 10.982 | −0.39% | −0.19% | −0.58% |
| u_y amplitude (mm) | 34.175 | 34.285 | 34.289 | +0.32% | +0.01% | +0.33% |
| u_y frequency (Hz) | 5.5227 | 5.5013 | 5.4907 | −0.39% | −0.19% | −0.58% |
| drag mean (N/m) | 459.74 | 459.44 | 459.08 | −0.07% | −0.08% | −0.14% |
| drag amplitude, extrema (N/m) | 28.03 | 27.01 | 26.49 | −3.64% | −1.91% | −5.49% |
| drag amplitude, smoothed | 27.63 | 26.28 | 26.11 | −4.89% | −0.65% | −5.50% |
| lift amplitude, extrema (N/m) | 174.14 | 157.32 | 152.67 | −9.66% | −2.96% | −12.33% |
| lift amplitude, smoothed | 172.27 | 155.78 | 150.95 | −9.57% | −3.10% | −12.37% |

(Featflow L4: u_x −2.88 mm, u_x ampl. 2.72, u_y ampl. 34.99 mm, f 5.46 Hz,
drag ampl. 27.74, lift ampl. 153.91 N/m.) Cycle-to-cycle spread over the
window is ≤ 0.02 mm in u_y amplitude (0.06%) and 0.7–1.2 N/m in the lift
amplitude, so the displacement changes of 0.3–1% are resolved but small; the
late-time variability (1 s scale wandering) is ~0.3% in u_y amplitude, so the
u_x changes of ~1% are marginal.

Comparison with the matched path (1x → 2x → 4x, fluid and solid refined
together, dt halved): u_x mean +24.9% / +10.6%, u_y amplitude +14.5% /
+6.0%, frequency −1.17% / −0.87%.

### Interface cleanliness

| solid | solid interface faces | fluid interface faces | mean / max outer iterations | steps above tol. | force balance (rms(F_f+F_s)/rms F_f) |
|---|---:|---:|---:|---:|---:|
| 2x | 432 | 168 | 7.01 / 13 | 0 | 5.9e-4 |
| 4x | 864 | 168 | 6.99 / 13 | 0 | 5.4e-4 |
| 8x | 1 728 | 168 | 6.98 / 13 | 0 | 5.0e-4 |

The solid interface is already 2.6× finer than the fluid one on the baseline
and 10× finer at 8x. The AMI force transfer stays conservative to 5e-4 (it
improves slightly), the iteration count is unchanged, every step meets the
tolerance and no instability occurs. Wall time grows (44 → 93 → 154
core-hours) because of the solid solve. The experiment is clean. A fine-solid /
coarse-fluid non-matching interface shows no material artefact in these
metrics; a possible effect on the traction interpolation itself was not
separately probed.

### Comparison with the CSM prediction

**It did not respond as predicted.** Section 12 suggested ~5–6% in u_y from
solid 2x → 4x (CSM2 tip error 6.7% → 1.7%) and ~1–2% from 4x → 8x. Observed:
u_y amplitude +0.32% and +0.01%, u_x mean +0.33% and −0.95% (u_x does not
even keep one sign). The frequency does respond monotonically and
decelerates (−0.39%, −0.19%), consistent with a softer plate, but its
magnitude (−0.58% over 2x → 8x) is far smaller than the ~3% a 6.7% static
stiffness error would produce for a pure structure (f ∝ √k). The
fluid-loaded plate evidently has a much lower sensitivity to structural
stiffness than the CSM2/CSM3 structure alone (the plate here is not close
to a static, gravity-loaded cantilever; a plausible but untested reason is
that the dynamic response is set by the fluid forcing and the added mass,
with the stiffness entering mostly through the frequency). The 4x → 8x
change is small relative to 2x → 4x for the frequency and u_y, i.e. the
solid has converged in the sense of these QoIs; for u_x it does not
shrink (−0.09% → −0.95%), but stays ≈ 1%.

### Error decomposition

- **u_y amplitude:** the full 1x → 2x → 4x path moves +14.5% and +6.0%.
  Solid alone moves +0.3% (2x → 8x), i.e. **≈ 5% of the 2x → 4x change**.
  The remainder is the fluid: fluid 2x → 4x at the 4x solid (the existing
  matched 4x run against this study's solid 4x) is +5.7%, which reproduces
  the whole matched change (0.3% + 5.7% = 6.0%).
- **u_x mean:** solid alone −0.6% against a matched +10.6%; the fluid
  accounts for +10.3%. The solid explains none of it, and has the wrong sign.
- **Frequency:** solid alone −0.58% (2x → 8x), against −0.87% matched; the
  solid explains about **two thirds** of the matched frequency change. The
  fluid-only change at solid 4x is −0.50% (not additive with the above
  because the two sequences share the solid 2x → 4x step, but the order of
  magnitude agrees).
- **Forces:** the solid refinement *does* change the loads: lift amplitude
  −12.3% (174.1 → 152.7 N/m, essentially onto Featflow's 153.9) and drag
  amplitude −5.5% (28.0 → 26.5 against 27.74), with drag mean unchanged. In
  the matched path these two shifts are hidden by the opposite fluid
  contribution (lift +12.5% fluid-only, 1x → 2x −24%, so the matched
  sequence is non-monotone in the extrema amplitude). The lift amplitude of
  the matched 4x run (+15%, reported as inflated by coupling noise in §4) is
  thus partly a coarse-solid effect, not only noise. This is a force result;
  it was not pursued (no force-noise campaign).
- The additivity is approximate: the matched 4x run also halves `dt` (a
  ≤ 1.7% displacement effect, §5/§7), and the pair f2/s4x → f4/s4x differs
  from fluid-only by that time-step change.

### Consequence for the apparent p ≈ 1 behaviour

The solid is **not** the origin of the displacement mesh dependence, nor of
the p ≈ 1 path: ~95% of the u_y change and all of the u_x change on the
matched 2x → 4x step come from the fluid (with the interface/coupling it
drives). The coarse-solid stiffness error identified in section 12 is real
and visible in the frequency (about two thirds of its 2x → 4x drift) and in
the lift/drag amplitudes (−10 … −12%), but it does not set the displacement
amplitudes in this FSI. The statement of section 12 that the solid "likely
explains much of the size and sign of the FSI3 displacement changes" is
**not supported**: the coincidence of the normalised CSM2 deflection ratios
with the u_y amplitude ratios (0.80 : 0.95 : 1 vs 0.82 : 0.94 : 1) was a
coincidence. No formal order is inferred; the displacement p ≈ 1 along the
matched path is, on this evidence, a fluid/interface property.

### Answers to the open questions

1. u_x mean +0.33% (2→4), −0.95% (4→8); amplitude −0.09%, −0.95%.
2. u_y amplitude +0.32%, +0.01%.
3. Frequency −0.39%, −0.19% (−0.58% total).
4. No, for the displacements (predicted 5–6%, seen 0.3%); partly yes for
   the sign and monotone decay of the frequency.
5. Yes for u_y and frequency (0.01%, 0.19%); u_x stays ~1%.
6. ≈ 5% of u_y, none of u_x, ~⅔ of the frequency drift, and the solid *does*
   explain a large part of the lift/drag amplitude shifts.
7. No. The p ≈ 1 sequence is fluid-dominated.
8. No material artefact (conservation 5e-4, iterations unchanged).
9. Yes: the standard solid is adequate for the displacements already at 2x
   against these QoIs; for force amplitudes a 4x–8x solid is needed.
10. Not tested: the optional `highOrderResidual` run was skipped because the
    structural error does not dominate the displacements (the instruction made
    it conditional on that).
11. See below.

### Classification

**C — little structural effect** on the displacement QoIs (the quantities of
primary interest), with a clearly visible structural effect on the force
amplitudes and a partial one on the frequency. Not D: the experiment is
clean.

### Recommended next single task

The CFD3 sub-benchmark on the 1x/2x/4x fluid meshes (rigid plate, no
coupling), to test whether the first-order behaviour of the loads and, by
extension, the fluid-driven displacement change comes from the fluid
discretisation alone. Run it with the dt of the matched path first, and only
then with a fixed dt to separate space and time.

## 14. CFD3 fluid-only convergence diagnosis

Asked for by the section 13 recommendation: does the fluid discretisation
alone show sub-nominal convergence in the loads that drive FSI3? Code and
compact results: solids4foam branch `verification/hronturek-cfd3` (commit
`00065a6fb`; `scripts/hron_turek_cfd3.py`, `scripts/cfd3_analysis.py`,
`reference/cfd3/`). Machine-readable: `hronturek_cfd3_paths.csv`,
`hronturek_cfd3_runs.csv`, `hronturek_cfd3_study.json`. OpenFOAM v2512,
xenosim. No `.tex` change, no FSI or solid run.

### CFD3 definition

CFD3 of Turek and Hron (2006) is the **rigid-flag version of the FSI3 problem
on the same geometry**: channel 2.5 × 0.41 m, cylinder of radius 0.05 m at
(0.2, 0.2), flag 0.35 × 0.02 m attached to the cylinder with its tip at
(0.6, 0.2), parabolic inflow of mean 2 m/s (maximum 3 m/s, `Re = 200`),
ρ = 1000 kg/m³, ν = 0.001 m²/s, no-slip walls, cylinder and flag, zero-stress
outflow. It is **unsteady**: periodic vortex shedding (St ≈ 0.22). Forces act
on the cylinder and flag together. Nothing was substituted: the case uses the
tutorial's own fluid `blockMeshDict` (the FSI3 fluid domain with the flag
cut-out), inlet and outlet, with the plate patch made a fixed wall; the flag
is the undeformed FSI3 flag. The published CFD3 reports drag and lift only
(mean ± amplitude [frequency]); there is **no pressure-difference quantity**,
so none is reported here. Initial condition: fluid at rest, impulsively
started; the flow takes ≈ 15 s (1x, `dt = 10⁻³`) to saturate (lift
envelope grows from the small asymmetry of the channel), so each run goes to
25 s (1x matched to 40 s) and the closing 1 s (3–4 lift periods) is analysed
with the FSI3 extraction (mean = (max+min)/2, amplitude = (max−min)/2 over
the last full lift period, 4 ms moving-average amplitude as the smoothed
variant). Periodicity: per-period lift amplitude varies by ≤ 0.19 N/m
(0.07%) over the closing 3 s, and its drift over those 3 s is ≤ 0.12 N/m, in
every run.

### Reference provenance

Featflow CFD3 table (`fsi_cfd_tests.html`, Q2/P1 discontinuous FEM, levels
1+0 … 4+0, 576 … 36 864 elements), best values level 4+0, `dt = 0.005`:
drag 439.45 ± 5.6183 [4.3956 Hz], lift −11.893 ± 437.81 [4.3956 Hz]. At
`dt = 0.01` the same level gives 439.38 ± 5.4639 [4.3825], −9.9868 ± 434.79.
The reference is **not temporally converged**: halving dt changes the lift
amplitude by 0.7%, the drag amplitude by 2.8% and the frequency by 0.3%; the
level sequence at `dt = 0.005` is monotone for the drag mean (437.41,
439.05, 439.45; ratio 0.24) but **not** for the drag amplitude (5.586,
5.580, 5.618) and the lift amplitude steps (434.74, 436.17, 437.81) do not
decrease (1.4, 1.6). I therefore treat the reference as good to ≈ 0.1% in
the drag mean, ≈ 0.5–1% in the lift amplitude and frequency, and 2–3% in
the drag amplitude. No reference is exact.

### Mesh and timestep families

| level | fluid cells | interface (plate+cylinder) | cells vs Featflow | `dt` matched / fixed | mesh quality (1x) | core-h |
|---|---:|---|---|---|---|---:|
| 1x | 5 336 | same blocks as FSI3 1x | between Featflow 1+0 and 2+0 (576, 2 304 quads, but Q2) | 10⁻³ / 2.5·10⁻⁴ | max non-orth 24.5°, skew 0.41, aspect 17 | 0.5 / 1.4 |
| 2x | 21 344 | ×2 per direction | – | 5·10⁻⁴ / 2.5·10⁻⁴ (+ 1.25·10⁻⁴) | – | 11.4 / 17.2 / 51.6 |
| 4x | 85 376 | ×4 | – | 2.5·10⁻⁴ (matched = fixed) | – | 187.6 |

The meshes are the FSI3 fluid meshes (blocks multiplied by 1, 2, 4, one
spanwise cell, graded as in the tutorial): `h` halves each level. Schemes
and solver settings are those of the FSI3 fluid (backward, linearUpwind
cell-limited, `nOuterCorrectors 3`, `nCorrectors 3`). Because the matched
4x is the fixed-dt 4x, six runs cover everything: 1x at 10⁻³ and 2.5·10⁻⁴,
2x at 5·10⁻⁴, 2.5·10⁻⁴ and 1.25·10⁻⁴, 4x at 2.5·10⁻⁴.

### Matched space-time results

| QoI | 1x | 2x | 4x | Featflow L4 | err 1x / 2x / 4x | ratio | order | verdict |
|---|---:|---:|---:|---:|---|---:|---:|---|
| drag mean (N/m) | 446.01 | 441.81 | 440.15 | 439.45 | +1.49 / +0.54 / +0.16% | 0.40 | 1.34 | monotone, sub-nominal |
| drag amplitude | 4.062 | 5.527 | 5.732 | 5.618 | −27.7 / −1.6 / +2.0% | 0.14 | 2.84 | monotone, order > 2 |
| drag amplitude, smoothed | 4.029 | 5.507 | 5.721 | – | | 0.15 | 2.79 | as above |
| lift mean | −67.9 | −19.0 | −13.8 | −11.89 | | – | – | near zero: undefined (converging) |
| lift amplitude (N/m) | 296.0 | 424.4 | 438.0 | 437.81 | −32.4 / −3.05 / +0.05% | 0.106 | 3.24 | monotone, order > 2 |
| lift amplitude, smoothed | 295.1 | 423.8 | 437.7 | – | | 0.108 | 3.21 | as above |
| frequency (Hz) | 4.415 | 4.454 | 4.448 | 4.3956 | +0.44 / +1.33 / +1.20% | −0.15 | – | non-monotone |
| Strouhal `f D/Ū` | 0.2208 | 0.2227 | 0.2224 | 0.2198 | | | – | non-monotone |
| drag mean, pressure part | 401.5 | 395.8 | 394.3 | | | 0.25 | 2.00 | monotone |
| drag mean, viscous part | 44.29 | 45.94 | 45.79 | | | −0.09 | – | non-monotone |
| lift amplitude, pressure | 294.5 | 421.7 | 434.8 | | | 0.10 | 3.28 | monotone |
| lift amplitude, viscous | 3.69 | 5.53 | 5.74 | | | 0.11 | 3.14 | monotone (1% of the total) |

### Fixed-dt spatial results

`dt = 2.5·10⁻⁴` on all three meshes (the 4x column is the same run).

| QoI | 1x | 2x | 4x | ratio | order |
|---|---:|---:|---:|---:|---:|
| drag mean | 445.49 | 441.78 | 440.15 | 0.44 | 1.18 |
| drag amplitude | 3.926 | 5.520 | 5.732 | 0.13 | 2.91 |
| lift amplitude | 284.6 | 423.7 | 438.0 | 0.103 | 3.27 |
| frequency (Hz) | 4.411 | 4.454 | 4.448 | −0.13 | non-monotone |
| drag mean, pressure | 401.0 | 395.8 | 394.3 | 0.28 | 1.83 |
| drag mean, viscous | 44.33 | 45.90 | 45.79 | −0.07 | non-monotone |
| lift amplitude, pressure | 283.1 | 421.0 | 434.8 | 0.101 | 3.31 |

Matched against fixed dt: the 1x→2x lift-amplitude change is +43.4% matched
and +48.9% fixed (the 1x matched run is the coarsest in time, which partly
masks the coarse-mesh deficit), and the 2x→4x change is +3.2% against +3.4%.
Orders are 1.34/1.18 (drag mean), 2.84/2.91 (drag amplitude) and 3.24/3.27
(lift amplitude): **fixing dt does not change the conclusion**.

### Temporal control

2x at `dt = 5·10⁻⁴, 2.5·10⁻⁴, 1.25·10⁻⁴`: lift amplitude 424.44, 423.67,
421.99 (−0.18%, −0.40%), drag mean 441.81, 441.78, 441.73 (−0.01%),
drag amplitude 5.527, 5.520, 5.505, frequency 4.4542, 4.4538, 4.4525 Hz
(inside the cycle spread). **The differences grow as dt falls** (ratios
1.6–2.2): there is no temporal convergence to quote an order from, and the
time error is not formally resolved. Its total (0.6% in the lift amplitude,
0.01% in the drag mean) is nevertheless 5–100× smaller than the 2x→4x
spatial change (3.2%, 0.4%). The persistent downward drift suggests a
dt-dependent splitting/tolerance effect of the PIMPLE loop (not
investigated; recorded as a negative result). **Temporal error is subordinate
for the lift amplitude and the drag mean, comparable to the 2x→4x change in
the drag amplitude only if the smoothed change (0.4%) is compared with the
spatial one (3.9%), i.e. still subordinate.**

### Observed order

Defensible: drag mean p ≈ 1.2–1.3 (matched 1.34, fixed 1.18; pressure part
1.8–2.0; viscous part not monotone); lift amplitude p ≈ 3.2 and drag
amplitude p ≈ 2.8–2.9, both **above the formal order** and therefore
pre-asymptotic error reduction (error −32% → −3% → +0.05% against the
reference), not asymptotic orders; no Richardson extrapolation is
meaningful. Not reported: frequency/Strouhal (non-monotone, differences
within 0.1–0.4% of the reference's own dt sensitivity) and the lift mean
(near zero). The 4x lift amplitude agrees with Featflow to 0.05%, but the
reference is uncertain by ≈ 0.7% in time, so this is agreement within the
reference uncertainty, not a demonstration that 4x is converged.

### Comparison with FSI3

| quantity | FSI3 1x→2x, 2x→4x (section 4) | CFD3 matched 1x→2x, 2x→4x |
|---|---|---|
| lift amplitude (smoothed) | 228.8 → 172.3 → 165.6 (−24.7%, −3.9%; ratio 0.12; p 3.1) | 295.1 → 423.8 → 437.7 (+43.6%, +3.3%; ratio 0.108; p 3.2) |
| drag amplitude (smoothed) | 22.6 → 27.6 → 30.6 (+22%, +11%; p 0.74) | 4.03 → 5.51 → 5.72 (+37%, +3.9%; p 2.8) |
| drag mean | 456.6 → 459.7 → 462.1 (within variability) | 446.0 → 441.8 → 440.2 (−0.9%, −0.4%; p 1.3) |
| frequency | 5.588 → 5.523 → 5.474 (p 0.4) | 4.415 → 4.454 → 4.448 (non-monotone) |
| `u_y` amplitude | +14.5%, +6.0% (p 1.07) | not a CFD3 quantity |
| `u_x` mean | +24.9%, +10.6% (p 0.91) | not a CFD3 quantity |

Directions of change differ (the FSI3 lift amplitude falls, the CFD3 one
rises: the fixed flag sheds much more strongly, 438 against 154 N/m), so
there is no one-to-one mapping. What carries over is the **shape**: the FSI3
smoothed lift amplitude reproduces the CFD3 error-reduction ratio (0.12 vs
0.108) and the order above 2 (3.1 vs 3.2), i.e. the lift-amplitude behaviour
in FSI3 looks like a fluid-resolution effect: a rapid pre-asymptotic
reduction with a 2x error of 3–12%. The CFD3 2x→4x load changes are small
(0.4% drag mean, 3.3% lift amplitude, 3.9% drag amplitude) and converging
faster than the FSI3 displacements (6% and 10.6%, p ≈ 1). The CFD3 drag
mean, the load that drives `u_x`, has a sub-nominal order (1.2–1.3) like
`u_x` (0.91), but a change 25× smaller than the `u_x` change. The fluid is
**pre-asymptotic at 1x and 2x** (2x errors: lift amplitude −3%, drag
amplitude −1.6%, drag mean +0.5%); at 4x it is within the reference
uncertainty. Whether load changes of this size can produce 6–10%
displacement changes in a flexible plate is not answered by a rigid
calculation.

### Diagnosis

**Classification: C, mixed.** The rigid fluid load converges cleanly
(monotone, error reduction faster than second order, 2x error of a few
percent for the amplitudes) except for the drag mean (p ≈ 1.2–1.3, but
only 0.4% in the 2x→4x change) and the frequency (non-monotone, ≈ 0.1%). It
does **not** reproduce the p ≈ 1 sequence of the FSI3 displacements, nor
their magnitude. Pressure dominates: the lift amplitude is 99% pressure and
converges like the total; the viscous drag part is non-monotone (±0.15
N/m) and 10% of the drag. The static fluid operator therefore does not by
itself explain the p ≈ 1 `u_x`/`u_y` behaviour. Per the limitation of a
rigid calculation, this does not clear the fluid in FSI3: CFD3 does not test
the ALE mesh motion, the moving no-slip boundary or the interface transfer.
The FSI3 lift amplitude behaves like the CFD3 one, so the fluid part of
FSI3 *force* convergence looks CFD3-like; the *displacement* convergence is
what departs from it. No optional fluid-side diagnostic was run (the
behaviour is not clearly sub-nominal).

### Answers to the questions

1. **CFD3:** the rigid-flag FSI3 problem, Re = 200, unsteady, drag/lift on
   cylinder + flag (see above); no pressure-difference quantity.
2. **Trust in the references:** good to ≈ 0.1% (drag mean), 0.5–1% (lift
   amplitude, frequency), 2–3% (drag amplitude); not temporally converged,
   non-monotone in the amplitudes.
3. **Monotone on 1x/2x/4x:** yes for the drag mean, drag and lift amplitudes
   and the pressure parts; no for the frequency/Strouhal and the viscous
   drag.
4. **Defensible orders:** drag mean 1.2–1.3; amplitudes 2.8–3.3 (above
   formal, pre-asymptotic).
5. **Fixed dt:** no change of conclusion (orders within 0.15).
6. **Temporal error:** subordinate (0.6% vs 3.2% for the lift amplitude) but
   not itself convergent (differences grow as dt falls).
7. **Drag and lift:** drag mean slowly (p ≈ 1.2–1.3, 0.4% change), amplitudes
   rapidly (3–4% change at 2x→4x).
8. **Pressure/viscous:** pressure drag p ≈ 1.8–2.0, viscous drag
   non-monotone (±0.15 N/m); lift pressure p ≈ 3.3, viscous lift is 1.3% of
   the amplitude.
9. **Reproduces the sub-nominal FSI3 displacement behaviour?** Not in
   magnitude or order; it reproduces the FSI3 lift-amplitude shape.
10. **Fluid discretisation the leading explanation for the displacement
    trend?** No, not on this evidence; it is a contributor (2x loads are 3%
    off in amplitude) that cannot be excluded for the displacement.
11. **ALE / interface transfer next?** Yes.
12. **Next task:** below.

### Recommended next single task

A fluid-only moving-mesh test without coupling: prescribe the 4x FSI3 plate
motion (the interface displacement history of one period, taken from the
committed 4x run or from a fixed periodic shape, e.g. the first bending
mode scaled to the tip amplitude) on the fluid meshes 1x/2x/4x with the
ALE velocity-Laplacian motion solver, and compare the plate
force/pressure-difference history across levels. This separates the
moving-boundary/ALE fluid response from the coupling and the structure
using the same machinery as CFD3. First step: check that the
solids4foam fluid model can prescribe the boundary motion (e.g. through
`pointMotionU`/`fixedValue` displacement patch); if not, the equivalent
alternative is a single FSI3 2x/4x pair with tight interface tolerance on
the interface transfer only.

## 15. Prescribed-motion fluid-only ALE test

Task recommended in section 14. Code and compact results: solids4foam branch
`verification/hronturek-ale-prescribed` (commit `27c9bfdf7`, stacked on
`verification/hronturek-cfd3`; `scripts/hron_turek_ale.py`,
`scripts/ale_analysis.py`, `reference/ale/`). Machine-readable:
`hronturek_ale_paths.csv`, `hronturek_ale_runs.csv`, `hronturek_ale_study.json`.
OpenFOAM v2512, xenosim. No solid, no coupling, no `.tex` change.

### Setup

Each run restarts the developed rigid CFD3 state of the same level (25 s,
section 14) and bends the flag with a prescribed periodic motion on the FSI3
fluid mesh with the FSI3 mesh motion (`velocityLaplacian`, quadratic
inverse-distance diffusivity on plate and cylinder, `newMovingWallVelocity` on the
plate, same schemes/solvers/PIMPLE). Prescribed transverse displacement
`η(x,t) = A s(t) φ(x) sin(2π f (t − 25))`, `A = 35 mm` (Featflow FSI3 tip
amplitude 34.99 mm), `f = 5.5 Hz` (FSI3), `φ` the clamped-free beam mode 2
normalised to unit tip value (FSI3 flaps at ≈ 2.5× the first in-vacuum
frequency, i.e. near the second mode), `s` a smooth 1 s ramp. The point
velocity is `(η(t) − η(t − dt))/dt`, so the Euler-integrated mesh follows η
exactly. Eight seconds are run (1 s ramp + 7 s); the analysis uses the last
five forcing periods (0.91 s). Mesh quality stays good at the extreme
deflection (4x: max non-orthogonality 31°, skewness 0.44; at rest 25°, 0.41).
Levels and time steps are those of section 14: matched 1x/2x/4x with
`dt = 10⁻³/5·10⁻⁴/2.5·10⁻⁴`, and fixed `dt = 2.5·10⁻⁴` (the 4x run is common). Five
runs, 0.2 to 109 core-hours. **Flag kinematics are an assumption** (a
standing mode-2 shape at the FSI3 tip amplitude and frequency, not the actual
FSI3 motion, which was not stored); the test therefore gauges the fluid
response to a given FSI3-like motion, not the FSI3 loads.

Quantities (per unit depth): force on cylinder + flag ("total") and on the flag
("plate"); the Fourier components at `f` split into the part in phase with the
displacement (`a1`) and with the velocity (`b1`); extrema-based amplitude over
the last period; the **generalised force** `Q`, the pressure force on the flag
projected on φ (viscous part excluded: < 1% of the load), and the work per
cycle `π A b1`. The flow locks to the forcing in every run (variance outside
the first six harmonics < 10⁻⁴; fundamental changes by ≤ 0.8 N/m, ≤ 0.1%,
between the last two five-period windows).

### Results

| quantity | 1x (dt 10⁻³) | 2x (5·10⁻⁴) | 4x (2.5·10⁻⁴) | 1x→2x | 2x→4x | order |
|---|---:|---:|---:|---:|---:|---|
| `Q` in phase with displacement `a1` (N/m) | 1391.5 | 1445.8 | 1443.4 | +3.90% | −0.16% | non-monotone |
| `Q` in phase with velocity `b1` (N/m) | −291.0 | −269.6 | −273.9 | +7.4% | −1.6% | non-monotone |
| work per cycle `π A b1` (J/m) | −31.99 | −29.64 | −30.11 | +7.4% | −1.6% | non-monotone |
| total lift amplitude, extrema (N/m) | 3572.4 | 3514.1 | 3498.0 | −1.63% | −0.46% | 1.85 (monotone) |
| plate lift amplitude, extrema | 3174.9 | 3266.5 | 3281.3 | +2.9% | +0.45% | 2.6 (above formal) |
| plate lift, fundamental `a1` | −3261.6 | −3355.2 | −3326.9 | +2.9% | −0.8% | non-monotone |
| plate lift, `b1` | −353.0 | −505.3 | −457.9 | +43% | −9.4% | non-monotone |
| total drag mean | 638.2 | 643.6 | 635.4 | +0.8% | −1.3% | non-monotone |
| total drag, 2nd harmonic | 169.0 | 181.3 | 185.4 | +7.3% | +2.3% | 1.6 (monotone) |

Fixed `dt = 2.5·10⁻⁴` (1x, 2x, 4x): `a1` 1398.4, 1445.9, 1443.4 (+3.4%, −0.17%);
`b1` −291.7, −271.8, −273.9 (+6.8%, −0.75%); total lift extrema 3567.1,
3514.2, 3498.0 (p 1.7); plate lift extrema 3189.2, 3264.7, 3281.3 (p 2.2);
total drag 2nd harmonic p 1.4. Temporal control (2x, `dt` 5·10⁻⁴ → 2.5·10⁻⁴):
`a1` 1445.77 → 1445.87, `b1` −269.6 → −271.8 (0.8%), total lift extrema 3514.1 → 3514.2.
Scale of the load: the modal inertia force of the FSI3 flag for this motion,
`ρ_s t L (∫φ²=¼) ω² A`, is **73 N/m**, against `Q ≈ 1.4 kN/m`.

### Observed order

No defensible order for the key modal load components: `a1` and `b1` of `Q`
change sign between differences (the 1x level is the outlier; the sequence
overshoots at 2x and falls back at 4x), as do the plate-lift quadrature
and the drag mean. Monotone sequences with orders in the second-order
range or above exist only for amplitudes that are dominated by the smooth
inertial part: the lift amplitudes (1.7–2.6) and the drag second harmonic (1.4–1.6).
No order is near 1 for a load that is changing by more than its
window-to-window noise (0.1%).

### Comparison with FSI3

2x → 4x, the modal pressure load changes by −2.4 N/m (`a1`, −0.16%) and
−4.3 N/m (`b1`, −1.6%), i.e. 3% and 6% of the 73 N/m structural scale;
1x → 2x by +54 and +21 N/m, i.e. 74% and 29% of it. The FSI3 `u_y` amplitude
changes by +14.5% (1x→2x) and +6.0% (2x→4x). The load differences
decrease much faster (ratio 0.03–0.2) than the FSI3 displacement
differences (0.41–0.48): the prescribed-motion fluid response is
converging near or above second order in sequence, not along the
p ≈ 1 sequence of the displacements. The load amplitudes that matter to
the structure are **small differences of large forces**: the flag is only
a ≈ 20th of the fluid force in inertia, so a relative load error of 1%
(≈ 14 N/m) is ≈ 20% of the structural inertia force. That amplification
(≈ 20 at the inertia scale; the nearest FSI3 evidence is the 5.7% `u_y` change from
refining only the fluid with the solid at 4x, section 13) is a plausible
route from sub-percent fluid-load changes to the 6–10% displacement
changes, but this test does not demonstrate it: the actual coupled
response, and the actual FSI3 kinematics, are not involved.

### Diagnosis

**Classification: B/C — the moving-mesh fluid response is not the p ≈ 1
source.** (B: ALE load converges cleanly and quickly for a given motion;
not A; not D: runs periodic, temporal error ≤ 0.8% in `b1`, 0.007% in `a1`.)
The mesh-motion, moving-wall and ALE flux machinery of FSI3 reproduces
smooth, fast load convergence on 1x/2x/4x for a prescribed flag motion;
the slow, first-order-like displacement sequence appears only when the
loads are fed back through the structure. The 1x level is
clearly pre-asymptotic (4–7% off in the modal loads); the 2x level is within
≈ 0.2% (`a1`) and 1.6% (`b1`) of the 4x one. The remaining candidates are
the interface/coupling treatment (AMI mapping between 168 fluid and 432–1728
solid faces, IQN-ILS) and the amplification of small load differences by
the coupled system, not the static or moving fluid operator.
Limitations: standing mode-2 kinematics, pressure-only generalised force, no
viscous part, no x-displacement, one amplitude, one frequency; the test says
nothing about the coupled response.

### Answers to the questions

1. **Did the ALE test show sub-nominal convergence?** No: the key modal
   loads are non-monotone with 2x → 4x changes of 0.16% and 1.6%; the
   monotone amplitudes have orders 1.4–2.6.
2. **Lock-in and periodicity?** Yes (variance outside harmonics < 10⁻⁴).
3. **Matched vs fixed dt?** Same conclusion (differences ≤ 0.8%).
4. **Is the p ≈ 1 FSI3 behaviour reproduced?** No.
5. **Fluid ALE the leading explanation?** No, on this evidence.

### Recommended next single task (superseded by section 16)

Coupled sensitivity (amplification) test, 2x only: in the existing FSI3 2x
configuration (dt 0.0005, solid 4x, IQN-ILS), scale the fluid traction passed
to the solid by (1 ± 0.01) and (1 ± 0.03) — two or four 2x runs of 7 s, about 10 core-hours
each — and compare the `u_y` and `u_x` amplitude changes with the 5.7% / 10%
changes of a fluid 2x → 4x refinement. If a 1–3% load change reproduces them, the
displacement sequence is a load-error amplification and the diagnosis
passes to the fluid load error at the interface (pressure/viscous traction
accuracy on the moving flag); if not, the interface transfer or the
coupling iteration is implicated. This requires a one-line traction scaling
in the interface (a driver-level change on a throwaway branch).

## 16. Fluid-only replay of the coupled FSI3 trajectory

Section 15 recommended a traction-scaling test. An independent review
(Codex Astra, one query) disagreed: a uniform scaling tests one error
direction only, the ratio of fluid force to structural inertia is not a
response amplification, and a failure to reproduce the trend would not
implicate the interface. It recommended replaying the **actual** converged
coupled interface trajectory on fluid 2x and 4x with identical motion and
comparing phase-resolved loads and cycle work. I adopted that. The tighter
tolerance and time-step controls it suggested as a cheap discriminator exist
already (sections 5 and 7). Code and compact results: solids4foam branch
`verification/hronturek-trajectory-replay` (commit `05af89c4c`, stacked on
the ALE branch; `hron_turek_replay_source.py`, `hron_turek_replay.py`,
`replay_analysis.py`, `reference/replay/`). Machine-readable:
`hronturek_replay_paths.csv`, `hronturek_replay_runs.csv`,
`hronturek_replay_study.json`. No solid in the replay runs, no `.tex` change.

### Setup

1. **Trajectory.** The completed FSI3 run fluid 2x / solid 2x (IQN-ILS, dt
   5·10⁻⁴) was restarted from t = 7 s for 1.2 s on 8 ranks, writing the
   positions of all plate-patch points at every step. The solid is total
   Lagrangian, so the restart uses `restart no` (the IQN-ILS history is rebuilt
   in a few steps); the restarted run reproduces the original periodic state
   (frequency 5.5236 against 5.523 Hz, tip amplitude 34.18 against 34.175 mm).
   The fluid 4x / solid 4x run could not be used: it did not write the restart
   state, and restarting costs ≈ 160 core-hours.
2. **Fit.** The last five periods (7.248–8.153 s, f = 5.5236 Hz) are fitted
   per point of the flag's surface polyline with a mean and four harmonics
   in x and y (maximum fit residual 0.08 mm, 0.2% of the amplitude); the
   coefficients are interpolated along the polyline (exactly on the 2x
   points, linearly in between for 4x, nested subset for 1x).
3. **Replay.** From the developed rigid CFD3 state of each level the plate
   point velocity is the finite difference of the replayed displacement
   (smooth 1 s ramp from 25 s), same mesh motion and fluid settings as FSI3
   and sections 14–15; 8 s run, last five periods analysed. Five runs:
   matched (1x/2x/4x with dt 10⁻³, 5·10⁻⁴, 2.5·10⁻⁴) and fixed dt 2.5·10⁻⁴
   (1x, 2x). 0.2–140 core-hours.
4. **Quantities.** Force on the flag ("plate") and on cylinder + flag
   ("total"), per unit depth; first harmonic split into the part in phase with
   the tip displacement (`inphase`) and with the tip velocity (`quad`);
   phase of the force; pressure work on the flag per cycle.

**Validation.** The 2x replay reproduces the coupled 2x total forces over the
same window: drag mean 459.5 against 458.4 N/m (0.2%), lift fundamental
166.9 against 167.3 (0.2%), drag second harmonic 25.7 against 25.6 (0.5%).
Replay runs are locked (variance outside the harmonics < 10⁻⁴); the
fundamental changes by ≤ 0.8 N/m between the last two five-period windows.

### Results (same motion, matched path, N/m)

| quantity | 1x | 2x | 4x | 1x→2x | 2x→4x | order |
|---|---:|---:|---:|---:|---:|---|
| total lift, fundamental | 138.9 | 166.9 | 168.4 | +20% | +0.9% | 4.3 (pre-asymptotic) |
| total lift, in phase with displacement | −101.7 | −82.6 | −49.2 | +19 | +33 | differences grow |
| total lift, in phase with velocity | 94.6 | 145.1 | 161.0 | +50 | +16 | 1.7 |
| total lift, phase (deg) | 137.1 | 119.7 | 107.0 | −17.4 | −12.7 | 0.45 |
| plate lift, fundamental | 176.1 | 178.7 | 164.0 | +2.6 | −14.7 | non-monotone |
| plate lift, in phase with displacement | −174.0 | −162.2 | −135.5 | +11.8 | +26.7 | differences grow |
| plate lift, in phase with velocity | 27.2 | 75.0 | 92.5 | +48 | +17.5 | 1.45 |
| plate lift, phase (deg) | 171.1 | 155.2 | 145.7 | −15.9 | −9.5 | 0.74 |
| plate drag mean | −13.5 | −16.2 | −18.4 | −2.8 | −2.2 | 0.31 |
| total drag mean | 464.8 | 459.5 | 456.0 | −5.4 | −3.4 | 0.64 |
| pressure work on flag per cycle (J/m) | −2.08 | −0.130 | −0.080 | +1.95 | +0.05 | 5.3 |

Fixed dt = 2.5·10⁻⁴ (1x, 2x, 4x): total lift in phase with displacement
−112.6, −81.7, −49.2; in phase with velocity 83.7, 146.7, 161.0 (order 2.1);
total lift phase 143.4°, 119.1°, 107.0° (differences −24.3°, −12.1°); plate
lift phase 174.4°, 154.6°, 145.7° (−19.8°, −8.9°); plate drag mean −13.2,
−16.3, −18.4 (order 0.47); total drag mean order 0.61. The 2x matched and
fixed values differ by ≤ 2% (e.g. plate quadrature 75.0 against 76.7).
The 2x→4x change is the same in both paths (4x is common).

### Observed order

- The **amplitude** of the total lift fundamental converges quickly (+0.9%
  2x→4x), and so does the pressure work per cycle in absolute terms
  (0.05 J/m).
- The **phase** of the load against the flag motion does not: total lift
  137° → 120° → 107° (matched; order 0.45), plate lift 171° → 155° → 146°
  (order 0.74); fixed dt 143° → 119° → 107° (order 1.0) and 174° → 155° → 146°
  (order 1.2). The in-phase component changes by +19 and +33 N/m (total), +12 and
  +27 N/m (plate): the differences grow with refinement, so no order is
  defined. The quadrature component converges at order 1.5–2.1 but changes by
  11–23% 2x→4x.
- The **mean streamwise force on the flag** changes by −14% 2x→4x (order
  0.3–0.5) and the total drag mean by −0.7% (order 0.6).

### Comparison with FSI3

The replay isolates the fluid response to a fixed motion, so any change is a
fluid-resolution effect. For the 2x→4x step the in-phase load on the flag
changes by 27 N/m and the quadrature one by 17.5 N/m. The structural modal
inertia scale for this motion is 73 N/m (section 15), so these are 37% and
24% of it, and the flag's mean streamwise load changes by 14%. These are
the quantities that fix the flag's equilibrium shape (`u_x` mean: mean
streamwise load, −14% against the FSI3 `u_x` change of +10.6%), frequency
(in-phase load: added-mass-like) and amplitude (quadrature load: damping or
negative damping). The sign and size of those changes are consistent with
the 6–10% FSI3 displacement changes, and their orders (0.3–1.5, differences
that do not shrink) are consistent with the observed p ≈ 0.9–1.1. The lift
amplitude, by contrast, is converged to 1% at 2x, which hides the phase
change.

### Diagnosis

**Classification: A, fluid-dominated, for the flag loads.** With
the **same imposed coupled motion**, the fluid load on the flag is
sub-nominally convergent in phase and in the mean streamwise component, while
the load amplitude and the CFD3 loads converge quickly. This reproduces,
without any solid or coupling, the character of the FSI3 displacement
sequence, and replaces the earlier hypothesis that the coupled system
amplifies small load differences with a direct fluid-side observation. Caveats:
one source trajectory (fluid 2x / solid 2x, 4x–1x interpolated along the
surface); the true coupled 4x motion differs by ≈ 6% in amplitude and will
differ in phase, so the replay does not give the FSI3 4x load, only the
fluid's response to a fixed motion; the 1x level is coarse; the work per cycle
is pressure only and is the difference of large terms; no formal order is
quoted for non-monotone or growing differences. This does not exclude a
coupling contribution on top.

A concrete candidate for a first-order source is visible in the setup, not
tested here: the FSI3 tutorial fluid uses `zeroGradient` pressure on the plate
(`p.dirichletNeumann`), whereas the library has `movingWallPressure`
(`∂p/∂n = −ρ n·a_wall`, used by other FSI tutorials). For a wall acceleration
`a_n ≈ ω² A ≈ 40 m/s²`, a zero-gradient wall pressure is in error by
`ρ a_n h/2` (first order in the near-wall cell size h), about 80 Pa at 2x or
≈ 30 N/m over the flag: the size of the observed in-phase load change. This
is a hypothesis from a size estimate, not a result.

### Answers to the questions

1. **Did the replay confirm that the fluid load on the moving flag is
   sub-nominal?** Yes, in phase and in the mean streamwise force; no in the lift amplitude.
2. **Does the replay reproduce the coupled forces?** Yes (0.2–0.5% at 2x).
3. **Matched vs fixed dt?** Same conclusion (≤ 2% at 2x).
4. **Is the p ≈ 1 FSI3 behaviour reproduced in the fluid?** Qualitatively, in
   the phase and mean-load orders (0.3–1.2).
5. **Is the coupled system needed to explain it?** Not required; not excluded.

### Recommended next single task

Test the wall-pressure condition on the replay: repeat the 1x/2x/4x replay (same
trajectory, same table, 8 s, matched dt) with `movingWallPressure` on the
flag instead of `zeroGradient`, and compare the same quantities. If the
phase and in-phase/mean loads then converge at order ≈ 2 and the 2x→4x
changes shrink, the first-order source is the wall-pressure boundary
condition of the FSI3 tutorial and the next FSI3 run is a matched 2x
sensitivity with the corrected condition. If they do not, the cause is in the
ALE flux/mesh-motion discretisation (e.g. the wall-flux and the
`backward` ddt on the moving mesh) and the next task is to switch those one at
a time on the replay. Cost: about 140 core-hours for the 4x run, 15 for
1x + 2x.
