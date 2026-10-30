# RING: window and estimator dependence of the wet/dry frequency ratio (last-two-period fit and resolution)

Status: post-processing of STORED displacement histories only (no new simulation, no source result changed). Independent evidence, draft.
Tools: `tools/ring_lastperiods_fit.py` (windows), `tools/ring_lastperiods_diag.py` (residual, harmonics, amplitude), `tools/ring_lastperiods_ext.py` (nuisance terms), `tools/ring_lastperiods_fig.py`.
Data: `results/ring_lastperiods_reproduction.csv`, `results/ring_lastperiods_windows.csv` (75 rows: strong/moderate/weak M2 at 100 steps/period, strong M2 at 50 and 25, moderate/weak M2 at 50, strong/weak M1; 13 windows each where available), `results/ring_lastperiods_extmodels.csv`, `results/ring_lastperiods_diagnostics.csv`. Figures (inspected): `figures/ring_lastperiods_windows.png`, `figures/ring_lastperiods_residual.png`.
Exact references (official JSON): dry 2.5550421958; wet 2.4356148807 / 1.8045131724 / 0.7686544910; ratio 0.9532581829 / 0.7062557226 / 0.3008382767.

## 1. The original procedure and its reproduction
Quantity: D(t) = (u_x(theta = 0) - u_y(theta = 90 deg))/2 at the outer-surface probes, 401 samples (every step), t_end = 4 T_ref. Estimator (tools/analyze_ring.py, via run_report.py, `0-4T` window): damped sinusoid d = exp(-s t)(a sin wt + b cos wt), (w, s) found by nested-grid variable projection over 0.9..1.1 w_ref; the ratio is omega_fit(wet) / omega_fit(dry) with the dry run of the same mesh (and the same steps/period); windows are given in periods of each run's own exact reference period.
Reproduction: calling the unchanged `AR.fit_damped` on the stored probes reproduces the stored `omega_fit` of all 13 runs with relative difference **0.0 (bit for bit)**. A fast scipy implementation of the same model (used for the many windows) agrees with the original to 3e-14 .. 7e-12 relative (max 7.1e-12, moderate M2 s50), i.e. 3 to 6 orders below any effect discussed here. The ratio error recomputed with exact JSON values: weak +0.00413 %, moderate +0.00723 %, strong +0.03987 % (the stored reports use rounded reference omegas: +0.03987 % strong; the PR tables +0.0399).

## 2. Last two periods and every other window (M2, 100 steps/period; ratio error in %, damped-sinusoid fit)
| window (periods) | weak | moderate | strong | strong omega_wet (rad/s) |
|---|---|---|---|---|
| 0-4 (original) | +0.00413 | +0.00723 | +0.03987 | 0.773787885 |
| 0.5-4 | +0.00413 | +0.00723 | +0.03989 | 0.773787992 |
| 1-4 | +0.00413 | +0.00724 | +0.03994 | 0.773788363 |
| 1.5-4 | +0.00414 | +0.00724 | +0.03998 | 0.773788723 |
| **2-4 (last two periods)** | **+0.00414** | **+0.00726** | **+0.04013** | 0.773789822 |
| 2.5-4 | +0.00416 | +0.00730 | +0.04034 | 0.773791482 |
| 3-4 (one period) | +0.00420 | +0.00741 | +0.04127 | 0.773798712 |
| 0-3 | +0.00413 | +0.00723 | +0.03991 | 0.773788192 |
| 0-2 | +0.00414 | +0.00726 | +0.04004 | 0.773789161 |
| 1-3 | +0.00414 | +0.00725 | +0.04006 | 0.773789340 |
| 0-1 / 1-2 / 2-3 (one period) | +0.00419 / +0.00419 / +0.00420 | +0.00736 / +0.00736 / +0.00740 | +0.04067 / +0.04074 / +0.04092 | |
- **Last two periods vs all four: +0.00026 pp in the strong ratio error (0.04013 vs 0.03987 %), +0.00003 pp moderate, +0.00001 pp weak.** The floor is the same: dropping the first two periods does not remove it. Start-up content is therefore not the main source of the floor (a start-up transient that decayed over the first periods would make the 2-4 fit differ by the floor itself, 0.04 pp, not by 0.0003 pp).
- Frequency change, strong wet, 4 -> last 2 periods: +1.94e-6 rad/s absolute, +2.5e-6 relative (+0.00025 %). Ratio change identical to the table. Dry (all windows 0-4 .. 2-4): omega_dry 2.571080630 .. 2.571080615, within 1e-8 relative.
- Window spread (all 6 windows of >= 2 periods: 0-4, 1-4, 2-4, 0-3, 1-3, 0-2): strong 0.03987..0.04013 (width 0.00026 pp, 2.6e-6 of omega), moderate 0.00723..0.00726 (0.00003), weak 0.00413..0.00414 (0.00001). One-period windows are ~0.0008 pp (strong) off, so they are not used.
- Drift with window start: the strong wet omega rises monotonically with the window start (0.7737879, ..., 0.7737898 at start 2T, 0.7737987 at start 3T with the single last period); a smooth drift of 1e-6 relative per period, not oscillating. The same direction in all regimes and at all dt and meshes (strong M1: 0.03545 -> 0.03585 %, strong M2 s50: 0.03971 -> 0.03995, s25: 0.03908 -> 0.03924).
- Time-study cases (M2, strong): 50 steps/period 0-4 +0.03971, 2-4 +0.03995; 25 steps/period +0.03908, +0.03924; as in the main study the ratio error is nearly dt independent and the window effect is small (<= 0.0002 pp for >= 2 period windows).

## 3. Estimators (strong M2, 100 steps/period)
| estimator | window 0-4 | window 2-4 | window 0-2 |
|---|---|---|---|
| damped sinusoid (original) | +0.03987 | +0.04013 | +0.04004 |
| undamped sinusoid | +0.03987 | +0.04011 | +0.04004 |
| zero crossings (least-squares line through crossing times, as in `analyze_ring.py`) | +0.03988 | +0.04044 | +0.03997 |
| extremum times (3-point parabolic interpolation) | +0.03973 | +0.04104 | +0.03914 |
Zero crossings and extrema are coarser (sampling dt = T/100, interpolation error), with a window spread of 0.0006 and 0.0019 pp respectively; the damped fit is the most stable. Damped and undamped fits agree (the decay is 2.2e-6 /s relative, 1.8e-3 % amplitude loss per period), so damping does not matter for the frequency. All estimators give the floor +0.040 % +/- 0.001 pp for windows >= 2 periods.
Nuisance terms (`results/ring_lastperiods_extmodels.csv`, strong M2, six windows): adding a constant, a linear drift and/or a third harmonic to the fit model moves the 4-period result from +0.03987 to +0.03979 and **collapses the window dependence**: with DC + linear + 3 omega all six windows give +0.03979 +/- 0.00001 %. The window drift above is therefore explained by slow and third-harmonic content in the wet displacement (section 4), and the remaining uncertainty of the floor is +/- 0.0001 pp. The floor is +0.0398 %, whatever the model.

## 4. Transient, second mode, damping, amplitude (`results/ring_lastperiods_diagnostics.csv`, `ring_lastperiods_residual.png`)
| M2 | fit residual rms / A (fundamental only) | 3 omega harmonic / A | 2 omega / A | DC / A | slow (0.05 omega) peak in residual / A | amplitude loss per period |
|---|---|---|---|---|---|---|
| dry | 2.0e-6 | 2.9e-8 | 7.9e-9 | 1.4e-9 | (negligible; peak at 23 omega, 2.8e-6: hex8 ring higher mode / numerical) | 8.6e-8 |
| weak | 2.6e-6 | 3.3e-6 | 3.9e-8 | 3.8e-7 | - | 4.1e-5 % |
| moderate | 6.0e-6 | 6.0e-6 | 1.6e-7 | 3.6e-6 | 7e-6 | 3.2e-4 % |
| strong | 4.6e-5 | 3.4e-5 | 1.8e-6 | 3.4e-5 | 6.7e-5 | 1.8e-3 % |
- Second mode / start-up: the residual of the single damped sinusoid is at most 4.6e-5 of the amplitude (strong), with a component at 3 omega (3.4e-5) and a slow drift/DC component (3.4e-5 to 6.7e-5), both **growing in time rather than decaying** (figure, panel a), so they are not a start-up transient. The dry ring has none of it (3e-8). Both are wet-only. The 3 omega component scales with rho_f exactly like the floor: h3/A = 3.3e-6, 6.0e-6, 3.4e-5 (x1.8, x5.6 weak -> moderate -> strong) against floors 0.0041, 0.0072, 0.0399 % (x1.8, x5.5). The 3 omega content is of the size of a very weak fluid-side nonlinearity, but a third harmonic does not by itself identify one (a different n = 2 mode of the wet system with omega about 3 omega_1 could be responsible); this is not resolved.
- Amplitude: the envelope sqrt(d^2 + (v/omega)^2) is constant over the four periods to the printed precision (e.g. strong 0.00194 in each period); amplitude A = 5.8e-4 (dry), 6.1e-4 (weak), 8.3e-4 (moderate), 1.94e-3 m (strong), i.e. A/a = 6e-4 .. 1.9e-3. Damping per period as in the table (decay 2.2e-6 /s relative strong), irrelevant for omega at the 1e-7 level.
- If the 3 omega content were a Duffing-type response, the implied relative frequency shift would be 12 x (h3/A) = 3.9e-5, 7.2e-5, 4.0e-4 (weak, moderate, strong). The measured ratio-error floors expressed the same way are 4.1e-5, 7.2e-5, 4.0e-4 (ratio error 0.0041, 0.0072, 0.0399 %). The two agree within about 5 % in all three regimes. This is a numerical coincidence worth stating, not a demonstrated cause: the 12 x h3 value is the shift a cubic oscillator with that third-harmonic content would have, whereas the third harmonic could have other origins (a different wet mode near 3 omega_1, a time-discretisation or enforcement effect), and the fit with an explicit 3 omega term removes only 0.0001 pp of the floor (so the 3 omega term is a symptom of the same physics, not a bias of the estimator). For a cubic fluid-side nonlinearity both the third harmonic (A^3) and the shift (A^2) would drop on halving the amplitude (shift x 1/4, h3/A x 1/4); this is exactly what the halved-amplitude run tests.

## 5. Resolution
- Data are deterministic (no sampling noise); resolution is limited by the fit-model (window/estimator) dependence, not by round-off. 4 periods at 100 steps/period (401 samples) with the damped-sinusoid fit reproduce the same omega to 1e-12 relative across implementations, and window choices (>= 2 periods) and nuisance terms change the strong ratio error by at most 0.0003 pp (2.6e-6 of omega) and 0.0001 pp respectively. A time-discretisation error that is common to wet and dry cancels in the ratio (ratio error varies by < 0.001 pp over 25..100 steps/period).
- Required: discriminate +0.0399 % from +0.0100 % (an A^2 effect scaling by 4), a change of **0.030 pp = 3.0e-4 in omega_wet**, and a null. The fit uncertainty is about **100x smaller** (0.0003 pp) than the step to be discriminated; a 10 % share of the floor that scales with A^2 (delta = 0.003 pp) is still 10x above the uncertainty.
- **Decision: the frequency resolution is sufficient** (figure panel c). The half-amplitude run is worth running; the pre-registered sigma is 0.0003 pp (RING_AMPLITUDE_DIAGNOSTIC.md, section 1.6).
- Caveat: the half-amplitude wet run is compared with the baseline through the same estimators; the last-two-period fit is the check that start-up content does not drive any difference.

## 6. Interpretation (facts vs hypotheses)
Facts: (i) the floor is reproduced for every window of >= 2 periods, every estimator (damped, undamped, zero crossings, extrema) and every fit model (+0.0398 .. +0.0401 %); (ii) start-up content (first periods) is not the source of the floor; (iii) the wet displacement contains small wet-only components (3 omega, slow drift), growing with time and with rho_f, of size 3e-6 .. 7e-5 of A; the dry ring has none.
Hypotheses (not established): a fluid-side amplitude effect would make those components and the floor scale with A^2; weak w.n = 0 enforcement would not. Unresolved: whether the 3 omega and drift components share a cause with the floor.

## 7. Follow-up (added after the run)
The half-amplitude run (job 14084) confirmed that the resolution was sufficient: delta = -0.0302 pp against a fit uncertainty of 0.0003 pp, identical for every window and estimator (RING_AMPLITUDE_DIAGNOSTIC.md, section 3).
