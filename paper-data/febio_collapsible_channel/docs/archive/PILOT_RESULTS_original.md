# collapsibleChannel pilot: RESULTS (CC agent, 2026-10-31)

Stock pristine FEBio v4.13 (sha256 f5ae1ed6f347f875c392f93d847f2369f437345a9b785296d46b0bd6fbb72f02), solids4foam `development` b200718, oomph-lib reference CSVs copied to `reference_data/cc/`. Pre-registration: `COLLAPSIBLE_CHANNEL_PILOT_PLAN.md` (incl. amendments A1, A2, both written before the run they affect). Audit: `COLLAPSIBLE_CHANNEL_EQUIVALENCE_AUDIT.md`.

## 1. Outcome in one paragraph
Like-for-like physics is achievable in stock FEBio for everything except two discretisation items, and BOTH were measured. (1) hex8 wall: the tutorial's aspect-20 wall (80x8) is 11.6 % too stiff in midpoint deflection (PROVEN incompatibility at the tutorial wall mesh); hex8 converges at second order to the beam (Richardson +0.01 %), 160x4 misses its 3 % gate by 0.054 pp (-3.054 %), 320x8 passes 1 % (-0.757 %). (2) Fluid: the fixed-collapsed-shape fluid test passes mass, steadiness and streamwise-parity gates but FAILS the pre-registered F1-vs-F2 force gate (first-mode load -8.63 %, above the 8 % stop value) through a wall-normal row-parity mode, and the single permitted diagnostic (20-step coupled run through the first trough) REFUTED the hypothesis that this is negligible in the dynamic collapse: the nodal wall pressure that FEBio applies to the wall exceeds the element-mean wall pressure by +19 % to +65 % (net force and first-mode load) during the loading phase, the wall-adjacent pressure rows alternate by about +-50 Pa on a smooth pressure of 130 Pa (36 %), and the wall midpoint lags oomph-lib by up to 12.8 % of the peak. By the pre-registered rule NO full CC-C1 (140 steps) was run. Coupled comparison with the reference: NOT ESTABLISHED; applied-load parity: DIAGNOSTIC DISCREPANCY (PROVEN mechanism and effect on the load actually applied). Allocated cost 0.195 of 4.0 CPU-h.

## 2. Job table (sacct, allocated CPU-h; my jobs only)
| job | case | cores | elapsed | CPU-s | CPU-h | outcome |
|---|---|---|---|---|---|---|
| 14214 | S1_80x8 | 2 | 2 s | 4 | 0.0011 | normal; midpoint -0.126920 m (-11.633 %) |
| 14215 | S1_160x4 | 2 | 2 s | 4 | 0.0011 | normal; -0.139242 (-3.054 %) |
| 14216 | S1_160x8 | 2 | 4 s | 8 | 0.0022 | normal; -0.139320 (-3.000 %) |
| 14217 | S1_320x8 | 2 | 5 s | 10 | 0.0028 | normal; -0.142542 (-0.757 %) |
| 14218 | F1_unpinned | 4 | 31 s | 124 | 0.0344 | normal, 40/40 steps, mean Newton 1.6 |
| 14219 | D1_fsi_20steps | 4 | 2 m 18 s | 552 | 0.1533 | normal, 20/20 steps, mean Newton 5.25 (max 6), no negative Jacobian |
| total | | | | 702 | **0.195** | cap 4.0; campaign tool: charged 2.930 CPU-h total all agents at the end of my work |
No job failed or was cancelled. Heavy lock `cc` held 04:20-end of D1 and released. Ledger rows (benchmark `collapsibleChannel`) S1_*, F1_unpinned, D1_fsi_20steps. Run dirs `runs/campaignCC/` (inputs, logs, Mon_* files, xplt, *_result.json). Input sha256: S1_80x8 5c782a2d..., S1_160x4 d4d97e82..., S1_160x8 b049c145..., S1_320x8 1d364a8b..., F1_unpinned 7da9da05..., D1 2d57045c0de9.... Generators: `tools/make_cc_feb.py` (0e0348346f1a...), analysis `tools/cc_s1_analyze.py` (9302d42d...), `cc_f1_analyze.py` (73f43ce4...), `cc_fsi_analyze.py` (9d472cae..., the version after one fix: volume-conservation sign and trough bracketing, made before interpreting the results).

## 3. CC-S1 static clamped wall (solid only; reference -0.14362832 m midpoint, -0.09645925 m quarter points)
| mesh (cols x rows, hex8) | aspect | midpoint (m) | error | quarter-point error | gate | result |
|---|---|---|---|---|---|---|
| 80x8 (tutorial) | 20 | -0.126920 | -11.63 % | -19.46 % | 3 % (decides the tutorial conformal mesh) | FAIL |
| 160x4 | 5 | -0.139242 | -3.054 % | -6.27 % | 3 % | FAIL by 0.054 pp (strict) |
| 160x8 | 5 | -0.139320 | -3.000 % | -6.17 % | informative | |
| 320x8 | 5 | -0.142542 | -0.757 % | -1.71 % | 1 % | PASS |
Order: error ratios 80->160 3.9, 160->320 4.0 (second order in the in-plane spacing; thickness rows 4 vs 8 change the result by 0.05 pp); Richardson limit -0.143642 (+0.01 % of the beam). Labels: tutorial-mesh incompatibility of stock hex8 PROVEN; convergence to the beam SUPPORTED (three meshes); the solid mesh used for the later stages (320x8) PASS 1 %. Amendment A1 documents the uncovered branch (160x4 marginal fail, 320x8 pass): coupled stages use 320x8 conformal with fluid 40/320/80 x 16 = 7,040 cells (the tutorial's 3,200 cells could not be used with a passing hex8 wall).
Free-body route used later: the S1 load-step curve (20 points, `Mon_bottom.txt` every step) gives the wall's quasi-static restoring pressure p_el(u_mid).

## 4. CC-F1 fixed collapsed shape, fluid only (7,040 cells, shape of S1_320x8, K = 1e7, dt 2.5 s x 40, rhoi 0; `f1_result.json`)
| quantity | observed | pre-registered gate | result |
|---|---|---|---|
| normal termination / Newton | 40/40, mean 1.6, no negative Jacobian | all steps | PASS |
| mass balance Q_out/Q_in - 1 | 6.0e-7 (Q_in 0.8810) | 1e-3 | PASS |
| steady drift F1, F2 (last 20 %) | 1.2e-11, 4.7e-12 | 1e-3 | PASS |
| streamwise parity a_t (wall / row1 / row2) | 2.4e-4 / 1.6e-4 / 2.2e-4 Pa = 6.8e-5 of V (3.51 Pa) | 5 % of V | PASS |
| F1 / F2 net F_y | 31.68 / 34.07 N/m (-6.99 %) | 3 % (3-8 conditional; > 8 stop) | conditional |
| F1 / F2 first-mode load G | 15.33 / 16.77 (-8.63 %) | same | FAIL (> 8 %) |
| junction pins | rule (c) false (drift 4e-12); rule (b) literally true by 7e-4 / 5e-4 Pa (monotone gradient, 1e-4 V), judged a non-event | pin if triggered | no pin run (recorded) |
Mechanism PROVEN from the stored nodal data: wall-normal y row parity, smooth in x: p(wall) = p(row 2) to 1e-3 Pa, p(row 1) differs by 0.47 .. -1.34 Pa (D = p0 - 2 p1 + p2: +1.09 .. -2.7 Pa, rms 1.75 Pa); F1 samples the even rows (applied nodal load), F2 averages even and odd. Absolute spread 2.4 N/m = 0.06 % of the dynamic scale (4.2e2 Pa x 10 m; eps_abs 5.7e-4), which is why a diagnostic of the dynamic regime was permitted (A2); the stop rule itself was not relaxed. Definitions of F1, F2, G, a_t, D: `tools/cc_f1_analyze.py` header.

## 5. Diagnostic D1: coupled 20 steps (A2 hypothesis H_D), 9,600 cells, dt 0.025, rhoi 0.9, K 1e7, E_mesh 1, 4 cores, 2 m 18 s
Numerical health: 20/20 steps, Newton 5-6 per step (mean 5.25), no negative Jacobian; volume conservation: max |int (Q_in - Q_out) dt - area change| = 0.24 % of the peak area change (PROVEN consistent); interface x-flux through the ends sign-reverses at t = 0.17 s (inlet flow reversal), as expected for a pressure-driven channel.
Applied-load spread (F1 nodal wall-row pressure vs F2 element mean, force per unit depth on the moving interface, N/m; `fsi_result.json` rows):
| t (s) | F1_y | F2_y | (F1-F2)/F2 | G spread | y-parity D (mean, Pa) |
|---|---|---|---|---|---|
| 0.10 | 791 | 663 | +19.4 % | +19.4 % | 51 |
| 0.15 | 1342 | 1078 | +24.5 % | +24.6 % | 106 |
| 0.20 | 1702 | 1301 | +30.8 % | +31.2 % | 161 |
| 0.25 | 1659 | 1172 | +41.6 % | +42.4 % | 196 |
| 0.30 | 1313 | 808 | +62.4 % | +65.1 % | 203 |
| 0.35 | 522 | 116 | (F2 crosses zero) | | 164 |
| 0.40 | -717 | -865 | -17.0 % | -18.6 % | 60 |
| 0.45 | -2093 | -1843 | +13.6 % | +14.4 % | -101 |
| 0.50 | -2842 | -2234 | +27.2 % | +29.3 % | -245 |
At t = 0.25 s the pressure at x = 10 on rows (wall, 1, 2, floor) is 193.4, 77.0, 192.9, 191.7 Pa: even rows 193 Pa, odd row 77 Pa (parity amplitude about 58 Pa on a mean of 135 Pa). Streamwise a_t stays 1e-3 Pa (the mode is purely wall-normal). H_D criterion (max spread <= 3 % over t >= 0.10 s with |F2_y| >= 100 N/m): maximum 349 % (zero crossing) and 19-65 % during the loading phase: H_D REFUTED (PROVEN by the data).
Solid free-body identification (route used because reaction forces are not available for FEBio fluid-FSI supports: HRON_TUREK_FSI1_STAGE_BC_RESULTS.md): with inertia negligible (rho_s h = 0.05 kg/m^2, 0.2 Pa) the wall's quasi-static balance gives the uniform-equivalent fluid pressure p_f = p_ext - p_el(u_mid), p_el from the S1 static curve; weighted with the static shape psi it is compared with the psi-weighted equivalents of F1, F2 and the odd rows (valid while |u_mid| <= 0.1425 m, i.e. t <= 0.35 s; shape error of order 10 % is expected at small amplitude):
| t | identified | F1 equiv | F2 equiv | odd-row equiv |
|---|---|---|---|---|
| 0.20 | 173.7 | 182.2 | 139.0 | 95.7 |
| 0.225 | 183.3 | 184.0 | 134.7 | 85.3 |
| 0.25 | 179.8 | 178.5 | 125.5 | 72.4 |
| 0.30 | 141.3 | 141.2 | 85.8 | 30.4 |
| 0.35 | 58.9 | 59.4 | 14.0 | -31.4 |
(Pa; sign: positive pushes the wall up). The load that moves the wall equals F1 to within 0.1-1 % for t = 0.225-0.35 s (at 0.075-0.15 s the quasi-static shape assumption is poor: 129 vs 143/115). So the applied load is the parity-contaminated nodal wall-row pressure, SUPPORTED/PROVEN (consistent with the source-level algebra of the earlier campaign). Which of F1/F2 is physically right is NOT resolved by this test (the true pressure is smooth across the 0.0625 m cell: physical cross-row variation = rho a_y h = 0.3 Pa, so neither the even nor the odd family is exact at the nodes); F2 (the element mean) is the cavity-validated better estimate (F2 = the mean of the two families) and the data says the applied load exceeds it by 20-60 % in the loading phase, with sign-dependent behaviour (it amplifies both the positive and the negative fluid load).
Wall response vs the reference (oomph-lib converged; wall midpoint, differences in % of the 0.216 m peak; informational only, 20 steps are NOT a validation): FEBio -0.0523 m at 0.25 s vs -0.0560 (oomph), -0.0874 vs -0.0979 at 0.30 s, -0.1248 vs -0.1451 at 0.35 s, -0.1764 vs -0.2041 at 0.425 s (maximum lag -12.8 % of peak), -0.1975 vs -0.2157 at 0.475 s, -0.2016 vs -0.2101 at 0.50 s; the FEBio wall is still descending at 0.5 s (first trough not bracketed: trough_reached = False), so its trough is later than the reference's 0.4685 s (> +6.7 %) and the amplitude at 0.5 s is 4 % of peak short. History error through 0.5 s: 12.8 % (max) / 6.5 % (rms) of the peak vs the converged reference; 11.3 % / 6.0 % vs oomph at the same step. For comparison the s4f tutorial's own history errors are 4.2-4.8 %. Direction and size agree qualitatively with the pre-registered load-balance law for a one-signed amplification of the applied load by 25-60 % (trough error 0.23 eps = 6-14 %, time/period 0.5 eps = 12-30 %); the early collapse (t <= 0.2 s, amplitudes < 0.025 m) agrees within 0.7 % of the peak, where the load is small. This is evidence consistent with the parity-bias mechanism, not a proof that it is the only cause (hex8 wall is 0.76 % too stiff, a few percent of dynamic effect from rhoi 0.9, dt 0.025 and K = 1e7 are not separated).

## 6. Classification (five terms)
| quantity | class | basis |
|---|---|---|
| static clamped wall, hex8 tutorial mesh 80x8 | DIAGNOSTIC DISCREPANCY (hex8 locking, -11.6 %) | S1 |
| static clamped wall, hex8 320x8 | QUALIFIED AGREEMENT (-0.757 %, 1 % gate; second-order convergence on 3 meshes, a second independent solver/reference is the oomph beam) | S1; not VERIFIED (one code, mesh sequence only in-plane; beam value is a different model, 0.05 % continuum offset) |
| steady frozen-shape fluid: mass balance, steadiness, streamwise parity | QUALIFIED AGREEMENT with its own declared gates (one mesh, no external reference; not VERIFIED) | F1 |
| applied wall load (nodal wall-row pressure vs element mean), steady and dynamic | DIAGNOSTIC DISCREPANCY (F1-F2 -8.6 % steady first-mode; +19..+65 % dynamic) | F1, D1 |
| first trough value and time, trough-to-trough period, rebound, second trough | NOT ESTABLISHED (no full run; the 20-step diagnostic does not bracket the trough) | stop rule |
| FEBio equal-order nodal-dilatation parity in a wall-driven added-mass flow | DIAGNOSTIC DISCREPANCY (mechanism PROVEN at the data level; first-principles cause as in FEBIO_COMMON_PRESSURE_MECHANISM.md, SUPPORTED) | D1 |

## 7. Limitations and uncertainties
- One mesh level, one time step, no K sweep, rhoi 0.9 only (one value), E_mesh 1 only; the 20-step diagnostic does not separate the parity bias from temporal/spatial/solid errors. The mapping from the applied-load amplification to the trough error uses the AE sensitivity law of the reassessment, not re-derived.
- The F1 gate's F_y/G measures are relative to a small steady force (a design weakness stated in advance); the dynamic test resolved the question in the other direction (the effect is LARGE in the dynamic regime, so the strict stop was not conservative-by-accident).
- The free-body identification uses the static shape and valid only at |u| <= 0.1425 m; it identifies the load that the wall sees (F1) with the S1 static curve (320x8 hex8, 0.76 % stiff).
- Inputs: the fluid mesh differs from the tutorial's (4x finer wall-region columns, 7,040 vs 3,200 cells) because the only passing hex8 wall is conformal with 320 columns.
- The pin rule (b) literal trigger at the two junctions was deemed a non-event without a run (1e-4 V). Not tested whether pins change F1-F2 (FSI1 showed they do not).
- hex20 wall + hex8 fluid (TestSuite fs19/fs20 route) was NOT run: S1 shows hex8 at 320x8 suffices for the wall, so it would not address the dominant fluid-side parity. It cannot remove the parity (hex20 fluid in `fluid-FSI` does not fix it: FEBIO_HEX20_FORMULATION_TEST.md).

## 8. What a full study would cost (NOT run, recommended only)
- C1 at 9,600 cells x 140 steps, 4 cores: from the measured D1 speed (6.9 s per step at 4 cores) about 16 min = 1.1 CPU-h; level 2 (4x cells, dt/2): about 8-10x = 9-11 CPU-h; time-step study (dt 0.0125, 0.00625): 2.2 and 4.5 CPU-h; K sweep 1.1 each. Given the measured parity bias a full run would only quantify the discrepancy, it cannot be a verification of the solver by itself. A scientifically meaningful coupled comparison needs an operator without the wall-normal parity (a stabilised or mixed-order fluid; none in stock v4.13: BENCHMARK_REASSESSMENT_AFTER_FSI1.md sec. 1-2), or a plugin, which is outside this pilot.
- Cheapest remaining information within stock FEBio (not authorised here): C1 full run to 3.5 s (1.1 CPU-h) to record the trough, period and decay as the parity-biased FEBio answer, labelled DIAGNOSTIC, plus one 2x finer-in-y fluid mesh F1 (0.3 CPU-h) to see whether the y-parity amplitude falls under refinement.

## 9. Reproducibility
```
python3 tools/make_cc_feb.py --mode solid --out runs/campaignCC/S1_320x8 --tag S1_320x8 --nsc 320 --ns 8          # likewise 80x8, 160x4, 160x8
python3 tools/make_cc_feb.py --mode fluid --out runs/campaignCC/F1_unpinned --tag F1_unpinned --nu 40 --nw 320 --nd 80 --nr 16 --k 1e7 --dt 2.5 --t-end 100 --rhoi 0 --shape runs/campaignCC/shape_S1_320x8.json
python3 tools/make_cc_feb.py --mode fsi --out runs/campaignCC/D1_fsi_20steps --tag D1_fsi_20steps --nu 40 --nw 320 --nd 80 --nr 16 --nsc 320 --ns 8 --k 1e7 --dt 0.025 --t-end 0.5 --rhoi 0.9 --e-mesh 1 --plot-stride 5
tools/submit_guarded.sh runs/campaignCC/<case> <case>.feb <cores> main <limit> cc_<case>   # 2 cores S1, 4 cores F1/D1
python3 tools/cc_s1_analyze.py runs/campaignCC/S1_*; python3 tools/cc_f1_analyze.py runs/campaignCC/F1_unpinned; python3 tools/cc_fsi_analyze.py runs/campaignCC/D1_fsi_20steps
```
`shape_S1_320x8.json` = dx, dy of the S1_320x8 final interface nodes from `Mon_bottom.txt` (the F1 shape file is part of the run record). Full C1 (140 steps): the same fsi command with `--t-end 3.5` (not run).
