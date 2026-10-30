# RING AMPLITUDE DIAGNOSTIC (halved initial amplitude, strong M2)

Status: independent FEBio v4.13 evidence, draft. Section 1 (pre-registration) was written and saved BEFORE the job was submitted. Sections 2 onward are filled after the run (see the end of the file).

## 1. Pre-registration (written before submission, 2026-10-30)

### 1.1 Question
Philip Cardiff (PR #26 comment 6528902475) excluded K, extra mesh stiffness and pure ring geometric nonlinearity as the cause of the mesh-independent wet/dry ratio-error floors (+0.0043 / +0.0076 / +0.041 % weak/moderate/strong). Remaining candidates: (a) a fluid-side amplitude effect (convection, moving domain; wet runs only), (b) weak w.n = 0 enforcement, (c) start-up content. Halving V0 discriminates (a) from (b)/(c): (a) predicts an amplitude-dependent floor, (b)/(c) do not (at fixed geometry they are amplitude-independent in a linear solve).

### 1.2 Case (the one remaining authorised ring simulation, 13 of 13)
- `R_AMP_strong_M2_V0half`: strong added mass (rho_f = 1400 kg/m3), mesh M2 (ring 8x96, fluid 24x96), K = 1e8 Pa, E_mesh = 1e-3 Pa, rhoi = 1, 100 steps/period (dt = 0.08174263586 s, from the exact period), 4 periods (400 steps), free-slip wall and interface, per-node IC workaround, pristine v4.13 binary, 4 cores. ONLY change: V0 = 0.75e-3 m/s (baseline 1.5e-3). Baseline: `runs/phaseD/D3_strong_M2` (Slurm 13720, 1 core).
- Generator: `tools/make_ring_feb.py` got a new option `--v0` (default 1.5e-3). Regenerating the baseline input with the modified generator (without `--v0`) reproduces the stored baseline hash bit for bit (df4b891a...fcef).
- Command: `python3 tools/make_ring_feb.py --out runs/ringNext/R_AMP_strong_M2_V0half --tag R_AMP_strong_M2_V0half --rho-f 1400 --k 1e8 --e-mesh 0.001 --rhoi 1 --mesh-factor 2 --steps-per-period 100 --periods 4 --ic potential --ic-field compatible --ic-mode pernode --wall slip --interface slip --v0 0.75e-3`
- Input sha256: `97f4926170f4ddd614050c83d05ef75c4ee0a4c9aa4ad828221a62dbab936d40` (baseline D3: `df4b891af2cb203ecd8a6c6e0fbd4714a249c7c7ecadbf72d45f99ae8f01fcef`).
- `diff` against the baseline .feb (18973 lines each): 4656 differing lines, all of them (i) the title comment, (ii) the ring velocity IC line (0.0015 -> 0.00075 in the expression) and (iii) the 4654 per-node fluid IC lines, whose values are exactly 0.5 x baseline (max relative deviation from 0.5x: 0.0, checked in python). Geometry, mesh, materials, BCs, time integration, solver settings and output are byte-identical.
- Expected cost: about 0.9 core-hours (R4 strong M2 on 4 cores: 786 s wall, 0.87 allocated core-h); expected wall 13 to 15 min.

### 1.3 Dry reference: reused, not re-run (decision and reason)
The dry M2 solution (`runs/phaseB/D1_dry_f2`, V0 = 1.5e-3) is reused for the half-amplitude ratio. Reason, from the actual model rather than an assumption:
- Constitutive law: FEBio `isotropic elastic` is St Venant-Kirchhoff in a total-Lagrangian, geometrically nonlinear solve, so amplitude independence is NOT guaranteed by the model.
- Measured bound from the stored data: the dry D(t) contains a third harmonic of 2.9e-8 relative to the fundamental (second harmonic 7.9e-9, DC 1.4e-9; `results/ring_lastperiods_diagnostics.csv`). For a Duffing-type softening/stiffening the relative frequency shift is 12 x (3rd-harmonic ratio) = 3.4e-7 at the baseline amplitude (5.8e-4 m), i.e. 3.4e-5 %, and it would change by at most that times 3/4 on halving. The ratio effect to be discriminated is 0.030 % (= 3e-4), about 900x larger. A quadratic nonlinearity gives a shift of order 15 x (h2)^2 = 1e-15. This is a **supporting estimate, not a proof**: it assumes that the third harmonic measures the cubic nonlinearity of the fundamental frequency (a Duffing-type argument); the dry fundamental frequency was NOT recalculated at half amplitude, so the reuse of the dry reference remains an assumption that a matched dry run (seconds) would remove.
- The dry ring's own dt (T_dry/100, from the exact dry period) does not depend on V0.
- A matched dry run at V0 = 0.75e-3 would take seconds but counts as a simulation and the bound above makes it unnecessary; this is a choice, and the bound is stated so that a reader can reject it. If someone wants it, the dry M2 run is `--v0 0.75e-3` without `--rho-f`, seconds on 2 cores.

### 1.4 Resolution of the measurement (Part 1, `RING_LAST_TWO_PERIOD_FIT.md`)
The original estimator reproduces the stored omega exactly (0.0 relative difference). The strong M2 ratio error is +0.0399 % with a fit-model uncertainty (windows, estimator, nuisance terms) of at most +/-0.0002 pp (all windows of length >= 2 periods, base model: +0.03987 to +0.04013 %; with nuisance terms: +0.03979 to +0.03986 %). The change to be discriminated (0.0399 -> about 0.0100 %, i.e. 0.030 pp) is about 100x larger than that, so a result is resolvable.

### 1.5 Quantities compared with the baseline (same estimator, same windows, same post-processing)
omega_wet, omega_wet^2 (and the equivalent additive shift of omega_wet^2 vs the exact value), wet/dry ratio and its error (original 4-period fit and last-two-period fit), the n = 2 pressure harmonic (scaled by 2), interface leakage (relative to the ring normal velocity), amplitude decay per period, fit residual and window spread, and solver convergence (all 400 steps converged, Newton iterations per step, normal termination from the log, no max_refs events, no negative Jacobian).

### 1.6 Decision rule (fixed now, not tuned afterwards)
Let delta = ratio_err(half) - ratio_err(baseline) (same estimator, percentage points). The uncertainty on delta from window/estimator choices is sigma = 0.0003 pp (conservative, 1.5x the largest spread observed in Part 1; the deterministic data have no statistical noise).
1. Amplitude-squared effect supported: delta < -3 sigma = -0.001 pp. The share f of the floor that scales with amplitude squared is then f = -delta / (0.75 x 0.03987); a full A^2 floor would give delta = -0.0299 pp (ratio error +0.0100 %), but the fourfold value is NOT assumed: the measured f and its uncertainty are reported, and a partial f is a partial result.
2. Amplitude hypothesis unsupported: |delta| <= 3 sigma. The floor does not depend on amplitude beyond the measured fit uncertainty (any A^2 share f is bounded by f < 0.001 / (0.75 x 0.03987) = 3.3 %); candidates left: weak w.n = 0 enforcement, start-up content, something else.
3. delta > +3 sigma: unexpected, reported as such (amplitude dependence of the opposite sign, or a confounder such as iterative-error scaling); not interpreted as support.
4. A result that depends on the estimator (the original fit and the last-two-period fit disagree about delta by more than 3 sigma) is reported as inconclusive. Nothing further is run.
Confounder noted in advance: the stopping tolerances (vtol, dtol, ftol ...) are partly absolute, so the relative iterative error could differ at half amplitude. R4 (100x looser tolerances) changed omega_wet by less than 1e-15 relative, so this is not expected to matter; the half-amplitude run does not test it.

### 1.7 Resources and process
Submitted through `tools/submit_guarded.sh` (v2 guard, 4 cores, partition main) after `squeue -u $USER` (empty) and `ps` (no febio4 of this user); ledger row pre-registered with `tools/ledger_add.sh` (status `submitting`) and updated with `tools/ledger_update.py` afterwards.

### 1.8 Secondary observable added before the result was read (job 14084 still running at this time)
Part 1 found a wet-only third-harmonic content h3/A = 3.4e-5 in the strong M2 baseline (3.3e-6 / 6.0e-6 weak / moderate), whose Duffing-equivalent frequency shift 12 x h3/A = 4.0e-4 equals the floor in all three regimes (RING_LAST_TWO_PERIOD_FIT.md, section 4). Prediction for the half-amplitude run if that content is a quadratic-in-amplitude fluid-side effect: h3/A falls to 1/4 (8.4e-6) together with the floor; if it is a linear (amplitude independent) feature, h3/A and the floor stay unchanged. h3/A is computed with `tools/ring_lastperiods_diag.py` logic (same model). This is recorded as a secondary check; the decision rule in 1.6 is unchanged and uses the ratio error only.

## 2. Run record (read from the log, `time.txt`, `report.json`; not from the exit code)
| | baseline D3 (V0 = 1.5e-3) | half amplitude (V0 = 0.75e-3) |
|---|---|---|
| Slurm job / cores / partition | 13720 / 1 / main | **14084 / 4 / main** (host xenosim) |
| input sha256 | df4b891a...fcef | **97f49261...6d40** |
| termination | normal | **`N O R M A L   T E R M I N A T I O N` in the log; exit code 0 only recorded** |
| steps converged | 400 / 400 | **400 / 400 (400 `converged at time` lines)** |
| Newton iterations per step | 4.55 (1819 residual evals) | 4.05 (1621); per step min/median/max 4 / 4 / 5 |
| max_refs events / Jacobian / NaN messages | 0 / 0 | 0 / 0 |
| convergence margins (worst achieved/tolerance) | see D3 margins | residual 1.3e-4, energy 6.9e-10, displacement 2.0e-3, velocity 2.8e-4, dilatation 0.945 (n = 400) |
| wall / CPU (user+sys) / allocated CPU-h / RSS | 2033 s / 1957 s / 0.54 / 346 MB | **741 s (12:21) / 2524 s (3.4 cores avg, 340 % CPU) / 0.82 allocated (0.70 used) / 359 MB** |
The pre-registered estimate (about 0.9 core-hours, 13 to 15 min) was met (0.82 allocated core-hours, 12.3 min). This was the 13th and last authorised ring simulation; nothing else was run. The dry reference was reused, not re-run (section 1.3).

## 3. Results (same estimators for both; dry = `D1_dry_f2`)
Data: `results/ring_amplitude_comparison.csv`, figure `figures/ring_lastperiods_amplitude.png` (inspected), tool `tools/ring_amplitude_compare.py`.
| quantity | baseline | half amplitude | change |
|---|---|---|---|
| omega_wet, original 4-period fit (rad/s) | 0.773787885 | **0.773554498** | -3.02e-4 relative (-0.0302 %) |
| omega_wet^2 (rad^2/s^2) | 0.598743 | 0.598382 | -3.6e-4 |
| wet error vs exact omega (0.768654491) | +0.66784 % | +0.63748 % | -0.0304 pp |
| **wet/dry ratio error, original 4-period fit** | **+0.03987 %** | **+0.00970 %** | **delta = -0.03017 pp** |
| ratio error, last two periods (2-4) | +0.04013 % | +0.00977 % | -0.03036 pp |
| ratio error, 0-2 window | +0.04004 % | +0.00974 % | -0.03030 pp |
| ratio error, DC + linear + 3 omega fit (all windows) | +0.03979 % | +0.00968 % | -0.03011 pp |
| equivalent additive shift of omega_wet^2 (2 omega_exact^2 eps) | 4.71e-4 | 1.15e-4 | ratio 4.1 |
| window spread of the fit (percent of omega) | 1.66e-4 % | 4.1e-5 % | |
| amplitude A of D(t) (m) | 1.9356e-3 | 9.680e-4 | 2.000 |
| n = 2 pressure harmonic (Pa) | 1.3752 | 0.6874 | x 0.4998 (linear) |
| 3 omega content h3 / A | 3.37e-5 | 8.43e-6 | **x 0.2503** |
| 2 omega / A ; DC / A | 1.77e-6 ; 3.38e-5 | 4.7e-7 ; 9.7e-6 | x 0.27 ; x 0.29 |
| interface leakage rms / ring normal velocity | 0.758 % | 0.438 % | x 0.58 |
| max abs w.n / V0 | 1.05 % | 0.28 % | x 0.26 (absolute x 0.13) |
| n = 2 part of the leaked flux / n = 2 ring flux | 0.310 % | 0.307 % | unchanged |
| numerical amplitude loss per period | 1.82e-5 | 1.23e-5 | x 0.68 |
| p phase (deg) / layer-0 added-mass ratio | 178.197 / 0.98886 | 178.197 / 0.98886 | unchanged |

Decision rule (section 1.6): delta = -0.0302 pp, sigma = 0.0003 pp, so delta = -100 sigma: **rule 1 holds: an amplitude-squared effect is supported.** All estimators and windows agree on delta (-0.0301 to -0.0304 pp, a spread of 0.0003 pp, inside the registered sigma, so rule 4 does not apply). The measured A^2 share of the floor is f = 0.03017 / (0.75 x 0.03987) = **1.01 +/- 0.01** (not assumed; measured). Fitting floor(A) = c + k A^2 through the two strong M2 points gives c = -0.0004 % (about zero within the +/-0.0003 pp window uncertainty) and k A_baseline^2 = 0.0402 %; i.e. the strong M2 floor is entirely consistent with a pure amplitude-squared term and no amplitude-independent remainder larger than about 0.001 pp (that is, the whole +0.040 % floor, and the 3 omega content, disappear as A -> 0 within resolution; this extrapolation rests on two amplitudes only).

## 4. Interpretation (facts, hypotheses, unresolved)
Facts (measured):
1. Halving V0 at fixed everything else lowered the strong M2 ratio error from +0.0399 % to +0.0097 %, a factor 4.1 (with a measured, not assumed, A^2 share 1.01 +/- 0.01). The result does not depend on the estimator (original fit, last two periods, 0-2, nuisance-term fit).
2. The wet third harmonic h3/A fell by 4.00 (0.2503), as an A^2 effect predicts for a cubic response; the n = 2 pressure harmonic scaled exactly linearly (x 0.4998).
3. Across the three regimes of the baseline, the ratio-error floor divided by the square of the wet amplitude is nearly constant: 0.00413 % / (0.611 mm)^2, 0.00723 % / (0.825 mm)^2, 0.03987 % / (1.936 mm)^2 = 1.11e4, 1.06e4, 1.06e4 % per m^2 (within 5 %). This explains why the floor grows about 10x from weak to strong and why it is a constant shift of omega_wet^2: A = V0 / omega_wet (an observation made after the run; the half-amplitude run is the independent confirmation).
4. The dry ring (no fluid) has no comparable content (h3/A = 2.9e-8); the floor appears only with the fluid.
5. The exact reference (official JSON) is the linearised, infinitesimal-amplitude solution. With V0 = 1.5e-3 m/s the wet amplitudes are 6e-4 .. 1.9e-3 m (6e-4 .. 1.9e-3 of the mean radius a = 1 m), and FEBio's finite-amplitude solution differs from the linear reference by +0.0043 / +0.0076 / +0.041 % in the frequency ratio.
Hypotheses (consistent with the data, not separated by it):
- A finite-amplitude fluid-side effect present in the FEBio formulation (convective term and/or moving-domain/ALE geometry nonlinearity, in the wet runs only) shifts omega_wet by an amount proportional to A^2. The data fix the scaling (A^2) and the location (wet only: the dry ring's small third-harmonic content suggests an amplitude effect of order 3e-7 on its frequency, an estimate, not a measurement at half amplitude), not the term.
- The weak w.n = 0 enforcement is not itself the cause of an amplitude-independent floor; its n = 2 content is amplitude independent (0.31 %), while its total leakage falls faster than linearly with amplitude (rms leak x 0.58 relative, x 0.29 absolute; max x 0.13), so the non-n = 2 leakage is an amplitude-nonlinear by-product and could be a symptom of the same nonlinearity. Whether the interface condition feeds the effect cannot be decided from this run.
- The Duffing-type equivalence 12 x h3/A = 4.04e-4 (baseline), 1.01e-4 (half) against the measured ratio-error offsets 3.99e-4 and 0.97e-4 is quantitative agreement (within 4 %) for a cubic-response mechanism; this is a consistency statement, not a mechanism identification.
Unresolved:
- Which fluid-side term (convection, ALE mesh motion, geometry update of the interface, weak interface condition) carries the A^2 effect; no run isolates it (and none was authorised).
- Whether solids4foam's finite-amplitude behaviour differs: its own ratio error also contains the amplitude content (not tested; its IQN-ILS errors are dominated by temporal error).
- Weak and moderate were not re-run; their amplitude dependence is inferred from the cross-regime scaling (fact 3), not measured. Only two amplitudes were run for strong M2; an A^2 law is the simplest description of two points plus the three-regime scaling, not a proof of exponent 2 (an exponent from the two points: ln(0.03987/0.00970)/ln 2 = 2.04).
- M1/M3 and the dt dependence at half amplitude were not run; the floor's mesh and time-step independence were established at the baseline amplitude only.
- The remaining +0.0097 % is consistent with the same A^2 law and with the extrapolated A -> 0 value -0.0004 +/- 0.0003 %; its sign and size below the half amplitude are an extrapolation from two amplitudes (the exponent 2.04 and the A -> 0 value are two-point quantities, not a rigorously verified zero-amplitude limit).
Implication for the benchmark statement (for the user to decide, not an action): the FEBio ratio error is within the official 0.1 % band at V0 = 1.5e-3 and its mesh-independent plateau is an amplitude effect of the problem definition (finite-amplitude response vs a linear reference), not a discretisation error that mesh or time refinement would remove. Because K, E_mesh, mesh, dt and tolerances were excluded earlier, and because the dry ring is unaffected, this supports the reading Philip proposed ("belongs to the problem definition, not to FEBio's discretisation"), with the caveat that the exact term is not identified.

## 5. Files
`runs/ringNext/R_AMP_strong_M2_V0half/` (input, logs, xplt, probes, reports, `slurm_14084.out`), `results/ring_amplitude_comparison.csv`, `figures/ring_lastperiods_amplitude.png`, `tools/ring_amplitude_compare.py`, `tools/make_ring_feb.py` (new `--v0`, default unchanged; baseline hash reproduced), ledger row `ringAddedMass,R_AMP_strong_M2_V0half`.
