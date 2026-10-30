# PHASE B: ringAddedMass pilot, INTERIM REPORT (stopped at the user's gate)

Date: 2026-10-29. **Gate status: NO-GO to the mesh/Δt/parameter matrix.** The wet-ring smoke test completes with valid solver residuals, but its numerical agreement with the exact reference is not acceptable, and the pressure field is non-physical. B5 and Phase C were deliberately not started. Acceptance thresholds are unchanged from `PHASE0_AUDIT.md` §9 (S1–S6); nothing was tuned.

## 1. The segmentation fault

**Symptom.** Wet-ring input (`runs/phaseB/B2_wet_K1e8_f1/`, job 13598): exit 139 after 0.13 s wall; log empty because the runner used `-silent`.

**Complete backtrace** (gdb on the pristine v4.13 binary, job 13601, `runs/phaseB/B2_wet_K1e8_f1_diag1/BACKTRACE_complete.txt`; 12 frames, release build so no line numbers):
```
#0 MSimpleExpression::value(...)                      libfecore.so
#1 FEMathExpression::value(FEModel*, FEMaterialPoint const&)
#2 FEMathValueVec3::operator()(FEMaterialPoint const&)
#3 FEInitialFluidVelocity::GetNodalValues(int, vector<double>&)   libfebiofluid.so
#4 FENodalIC::Activate()
#5 FEModel::Activate()
#6 FEModel::Init()
#7 FEBioModel::Init()
#8 febio::SolveModel(...)   #9 FEBioApp::RunModel()  #10 FEBioApp::Run()  #11 main
```
**When.** During **model initialisation** (`FEModel::Init` → `Activate`), at **physical time t = 0, before any time step and before any nonlinear iteration**; wall clock ≈ 0.13 s after the input was read ("Reading file … SUCCESS!" is the last line printed).

**Isolation** (`runs/phaseB/bisect/`, identical input apart from the fluid IC): constant-vector fluid IC (V1) runs; no fluid IC (V3) runs; trivial `type="math"` fluid IC `X*0.001,0,0` (V2) crashes. So the fault is independent of my expression, mesh, interface, material and BCs. The same `type="math"` mechanism works for the **solid** `velocity` IC (B1 ran).

**Root cause (tagged source).** `FEBioFluid/FEInitialFluidVelocity.cpp`, `FEInitialFluidVelocity::Init()` ends with `return true;`. The solid counterpart `FEInitialVelocity::Init()` ends with `return FENodalIC::Init();`. `FENodalIC::Init()` validates the node set and calls `FEInitialCondition::Init()`, which initialises the parameter valuators. Skipping it leaves the vec3 math valuator uninitialised, and evaluating it dereferences a null. A constant vector needs no initialisation, so it works. The same code is present on `master` (`67cc26e`). This is a bug in FEBio v4.13, not in the input.

**Exact fix (verified).** One line, `return true;` → `return FENodalIC::Init();` (`patches_febio_initialfluidvelocity.diff`). A copy of the tag with this patch (`upstream/febio-v4.13-patched`, binary sha256 `4b83bd96…`) runs V2 (the input that segfaults the pristine build) past initialisation and through 239+ converged steps (`runs/phaseB/bisect/V2_mathX_patched/`).

**Workaround used for all results below (pristine v4.13, no source change).** The generator evaluates the potential-flow velocity per node and writes one constant-vector `initial fluid velocity` IC per node (`--ic-mode pernode`, default). The values are the same analytic formula; the parser is not involved. `--ic-mode math` needs the patched binary. **All pilot results use the pristine binary** (sha256 `f5ae1ed6…`).

Original failed cases and logs are preserved: `B2_wet_K1e8_f1`, `_diag1` (gdb), `_diag2`, `_diag3` (two expression-syntax fixes that did not help: `sqrt()^3` → product; `--x`/`+-x` → parenthesised), and `bisect/`. The expression changes were a wrong lead, recorded rather than hidden.

## 2. Official fluid-FSI example on this exact build

Done in Phase A with the same binary: `fs01` (Couette; matches the analytic start-up solution to 1e-9 at steady state), `fs15` (inviscid, `rhoi=1`), `fs19`, and restart reproducibility (1.5e-10). See `PHASE_A_BUILD_AND_FSI_SMOKE.md`. Those did not use the fluid-velocity IC, which is why they could not reveal this bug.

## 3. Model audit (item 4)

Parsed from the generated `.feb` by an independent script (`tools/audit_ring_model.py`, output `runs/phaseB/B2_wet_K1e8_f1_pernode/model_audit.txt`):
- 192 ring + 576 fluid hex8; all Gauss-point detJ > 0 (min 9.8e-6 and 1.6e-5 m³).
- Element volumes 0.0157052 (ring) and 0.0901084 (fluid) m³ vs analytic annulus 0.0157080 and 0.0901244: **−0.018 %**, the chordal faceting of hex8 on arcs (a small geometric discretisation error).
- 98 nodes shared by ring and fluid domains, all at r = 1.05; the 48 interface quads are each a face of both a ring and a fluid element (conformal).
- BCs: u_z = 0 and w_z = 0 on all (plane strain); mesh normal displacement fixed and fluid normal relative velocity fixed on x=0 and y=0; wall nodes (r = 1.5): u = 0, w = 0; interface nodes w = 0 (no-slip, as in the reference); `fluid-FSI traction` on the interface surface.
- IC kinematics: expression vs polar formula 8.7e-19; radial velocity continuity at the interface 5e-20; zero normal velocity at the wall. **Measured from the run log at t = 0: ring KE = 1.0001 × analytic; fluid KE = 0.923 × potential-flow value** (lumped-mass estimate; plausible for a coarse mesh with wall and interface nodes forced to w = 0, not separately verified).

## 4. Dry ring (B1): independent post-processing

`tools/analyze_ring.py` on `runs/phaseB/B1_dry_hex8_f1/probes.txt` (QoI D = (u_r(0°) − u_r(90°))/2 at the outer surface; Slurm 13597, 400/400 steps converged, 12 s, 1 core; `rhoi=1`, 100 steps/period, 4 periods):

| window | ω zero-crossing | ω damped-harmonic LSQ fit | error vs 2.555042 |
|---|---|---|---|
| 0–4T | 2.621063 | 2.621063 | **+2.584 %** |
| T–4T | 2.621064 | 2.621063 | +2.584 % |
| 0–3T | 2.621063 | 2.621063 | +2.584 % |
| 0–2T | 2.621063 | 2.621063 | +2.584 % |

The two diagnostics agree to 1e-7 and the window sensitivity is ~1e-6. Amplitude loss 2e-9 per period (no numerical damping at `rhoi=1`); amplitude 5.715e-4 m = V0/ω as expected. **S1 (≤0.5 %) FAILS.** The ring is 2.6 % too stiff. The likely cause is hex8 bending/shear locking with four elements through the thickness (the solids4foam standard solid is 20.7 % too stiff in bending on its coarsest mesh). **This is a hypothesis; no element-type or refinement test has been run**, per your instruction not to start the matrix.

## 5. First wet-ring result (B2) and K runs (B3, B4)

Pristine v4.13, per-node potential-flow IC, `rhoi=1`, ρ_f = 140, ν_f = 1e-7 (μ = 1.4e-5 Pa·s), E_mesh = 1 Pa, tolerances dtol=vtol=ftol=1e-6, etol=1e-8, rtol=1e-10, full Newton, 1 core, 3.5 to 4 min each. All three: 400/400 steps converged, normal termination, no Jacobian or NaN messages.

| run | K (Pa) | ω fit 0–4T | error vs 1.804513 | wet/dry ratio (vs B1 dry) | error vs 0.706256 | amp loss/period |
|---|---|---|---|---|---|---|
| B3 | 1e7 | 1.729264 | −4.170 % | 0.659757 | −6.58 % | +7.7e-5 |
| **B2** | 1e8 | 1.728888 | −4.191 % | 0.659613 | −6.60 % | +5.9e-5 |
| B4 | 1e9 | 1.727008 | −4.295 % | 0.658896 | −6.71 % | −2.7e-4 (growth) |

The T–4T window changes ω by ≤ 0.012 % (B2). The B1 and B2 meshes use the same ring elements, so the wet/dry ratio cancels the ring's 2.6 % stiffening error.

**Per gate:** S2 (ω within 0.5 %, ratio within 0.3 %) **FAIL** for all three. S3 (K-convergence): ω changes by 0.022 % from 1e7 to 1e8 and −0.11 % from 1e8 to 1e9, but converges to a value 4.2 % from the exact reference, so it cannot support convergence to the incompressible value: **INDETERMINATE**. S4 amplitude loss ≤ 0.5 %: PASS (growth of 2.7e-4 in B4 is far below the s4f growth tolerance of 2e-3). S5 (no acoustic/pressure ringing): **FAIL** (below). S6 (archived): PASS (`run_meta_*.json`, logs, `time.txt`, input hash; `results/pilot_summary.csv`).

**What the K sweep rules out.** The −6.6 % ratio error is **not** a weak-compressibility effect. My Phase 0 estimate, a compressibility bias of order (ωb/c)² ≈ 1e-5 at K = 1e8, was a hypothesis; the measured K-dependence is ≤ 0.11 % over two decades, and the disagreement with the reference is ~40× larger and does not shrink with K. An effective added-mass ratio μ_eff = (0.6596)⁻² − 1 = 1.30 against the exact 1.008 means the discrete fluid inertia appears ≈ 29 % too large on this mesh.

## 6. Evidence of non-physical pressure and velocity fields (B2, B3, B4)

Interface-adjacent element at θ ≈ 0, harmonic fit at each run's own ω:

| run | K | fitted |p| amplitude | residual rms after removing the harmonic | max |p| in any element |
|---|---|---|---|---|
| analytic added-mass pressure m_a·ω·V0 | – | **0.325 Pa** | – | – |
| B3 | 1e7 | 0.27 | **11.1 Pa** | 36.2 |
| B2 | 1e8 | 1.29 | **1.07 Pa** (83 % of amp) | 4.1 |
| B4 | 1e9 | 13.4 | **6.6 Pa** | 28.8 |

The pressure is dominated by element-to-element oscillation 3 to 35 times the physical signal and is not monotone in K, so this is a pressure-field instability, not simply K × round-off. (Max |J−1| ≈ 4e-8 at K = 1e8.) Consistent with Phase A (pressure of order K·1e-9 in a pure shear flow). Also in B2 the ring-node velocity departs from the time-derivative of its own displacement by **6.7 % rms** (dry B1: 0.1 %), with a 2.3 % period-2 component, and the total kinetic energy at displacement zero-crossings wanders between 0.82 and 1.00 of its initial value. So the nodal velocity field carries non-modal content that the displacement-based frequency fit does not see. **S5 FAILS for B2, B3, B4.**

## 7. What is and is not established

Established: root cause and verified one-line fix for the crash; a working workaround on the pristine binary; the dry-ring frequency with two diagnostics (+2.584 %); a converged, undamped wet-ring run at three K; K-independence of ω; non-physical pressure noise in all three.

**Not established (no run has tested these):** why the added mass is ~29 % large. Candidate causes, ranked by my judgement, all untested:
1. Equal-order (Q1 velocity / Q1 dilatation) pressure instability at a coarse 12-cell gap (supported by the pressure noise; refinement or a different fluid element would test it).
2. The strong no-slip constraint on the interface and wall nodes, which over-constrains the tangential fluid motion in an inviscid limit (free-slip variant would test it).
3. ALE mesh-velocity handling with the massless regularising solid and `rhoi=1` (velocity inconsistency, §6; E_mesh and `rhoi<1` variants would test it).
4. Spatial discretisation of the potential flow on the coarse mesh (mesh factor 2 would test it).
These are the next diagnostics. They are *not* the B5/C matrix; each is one or two single-core runs of ~4 min, and I have not started them.

## 8. Proposed next step (needs your go-ahead; nothing started)

A short, bounded diagnostic set on 1 core each (≤ 5 runs, ~4 min each): (a) free-slip vs no-slip interface; (b) mesh factor 2 wet; (c) `rhoi` 0.9 and E_mesh 1e-3/1e3 on the coarse mesh; (d) a dry hex20-vs-hex8 structural check to test the locking hypothesis for S1. Stop criterion: the wet/dry ratio within 0.3 % *without* tuning, or a documented cause. If none reproduces the added mass, the honest result is that native FEBio fluid-FSI hex8 does not reproduce the ringAddedMass reference at these resolutions, which is itself reportable.

## 9. Files

`results/pilot_summary.csv`; `runs/phaseB/{B1_dry_hex8_f1,B2_wet_K1e8_f1_pernode,B3_wet_K1e7_f1,B4_wet_K1e9_f1}/` (inputs, logs, `probes.txt`, `analysis.json`, `run_meta_*.json`); failed/diagnostic: `B2_wet_K1e8_f1*`, `bisect/`; tools: `tools/{make_ring_feb.py,analyze_ring.py,audit_ring_model.py,run_febio_job.sh,submit_ring.sh}`; patch `patches_febio_initialfluidvelocity.diff`.

Parallel scaling (item 9): measured in Phase A (best 2.9–3.8× at 4–8 threads, 16 threads slower). The ring models are smaller than those test cases, so 1 core per run remains the plan; 16–32 cores per simulation are not justified and will not be used unless a measurement on a larger ring mesh shows otherwise.

---

## ADDENDUM 2026-10-29 (after the coupling audit and five-run diagnostic campaign)

The text above is the interim report as written and is **kept unchanged**. Updates that supersede or qualify it:

- **Superseded hypothesis (section 4 and 7).** The dry-ring error is *not* a locking plateau: the mesh-factor-2 dry ring gives +0.628 % (from +2.584 %), observed order 2.04 (`runs/phaseB/D1_dry_f2`). Compressibility was also ruled out (section 5).
- **Root cause of most of the wet-ring error found.** The strong nodal no-slip condition on the rigid outer wall: R1 (corrected IC, no-slip wall) wet/dry ratio error −6.01 %; D5 (same but free-slip wall) −1.11 %; spurious n = 0 pressure at ω 3.9 Pa → 0.000 Pa. A real input defect (IC incompatible with no-slip BCs, audit F1) explained only ~0.6 points.
- **Corrected input version v2** (`--ic-field compatible`) and the wall option (`--wall slip`) are in `tools/make_ring_feb.py`; the originals (B1–B4) are untouched.
- Gate status is unchanged in substance: **the 0.3 % screening tolerance is still not met (best −1.11 %)**, S1/S2 fail, S5 passes only for the free-slip wall run. Phase C has not been started.
- See `PHASE_B_COUPLING_AUDIT.md`, `PHASE_B_DIAGNOSTIC_FINDINGS.md`, `results/diagnostic_matrix.csv`, `NEXT_DECISION.md`.
