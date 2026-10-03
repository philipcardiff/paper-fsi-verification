# FSI verification benchmark inventory

Evidence inventory for a paper on verification benchmarks for continuum
fluid-structure interaction (FSI). Machine-readable twin:
[`benchmarks.yaml`](benchmarks.yaml).

## Provenance and ground rules

- **Source of evidence:** solids4foam `development` at commit `7c928c55`
  (2026-10-03, merge of PR #506), cloned fresh from GitHub, plus the merged and
  open pull requests listed below. Every number in this file comes from a
  case README, a `verification/README.md`, a reference JSON/CSV stored with the
  case, a PR description or comment, or `simulationData/` in this paper
  repository. Values marked *(derived here)* were computed in this inventory
  from tabulated values in those files; values marked *(read from plot)* were
  read off a versioned PNG and are approximate.
- **No new simulations were run** for this inventory.
- **Self-convergence and external agreement are kept separate.** "Observed
  order" means the order of solids4foam's own differences (or its error
  against an exact solution, where one exists). Agreement with a published or
  independent code is reported separately and never counted as convergence.
- **Reference categories:**
  - **A**: analytical or semi-analytical (no discretisation error of its own,
    or one controlled to round-off).
  - **B**: an independently converged numerical solution of *the same*
    problem, with its own convergence evidence.
  - **C**: a classical/published benchmark whose reference values carry
    unquantified or only partly quantified error. Weak, ambiguous or
    same-code-lineage references are flagged as **C-weak**.

### Pull requests examined

| PR | State | Content relevant here |
|---|---|---|
| #510 | merged 2026-10-02 | ringAddedMass tutorial + verification (exact added-mass reference) |
| #513 | merged 2026-10-03 | womersleyTube tutorial + verification (exact Womersley reference) |
| #507 | merged 2026-10-02 | collapsibleChannel tutorial + oomph-lib verification |
| #506 | merged 2026-10-03 | HronTurek FSI1/FSI2/FSI3 in one tutorial + verification |
| #501 | merged 2026-10-02 | 3dTube verification (Robin vs IQN-ILS; Tuković, Lozovskiy, Eken) |
| #497 | merged 2026-10-02 | mokLidDrivenCavity tutorial + verification |
| #498 | merged 2026-10-02 | blobInTreacle tutorial + verification |
| #412, #518 | merged 2026-09-09 / 10-02 | beamInCrossFlow mesh studies (both forms); IQN-ILS robustness |
| #439 | merged 2026-09-24 | migration of solid-benchmarks studies, incl. cavityFlexibleBottom |
| #450, #452, #492 | merged | Robin convergence criteria, automatic Robin coefficient, FSI predictor default |
| #491, #487 | merged 2026-10-02 | immersed-boundary FSI (immersedHronTurekFsi2); not body-fitted |
| #499, #308 | merged | membraneRoof (self-declared "not a verification"); flexibleDamBreak fixes |
| #233 | **open** | quasi-monolithic PETSc SNES FSI solver (the solver of the parent paper) |
| #503 | **open** | BDF3–6 `d2dt2` record (not for merging); not used by any benchmark |
| #528, #249, #116 | **open** | solidRobin BC, newtonIcoFluid SNES helper, damage models: no FSI verification content |

No **open** PR adds or extends an FSI verification case. All verification
evidence is on `development`.

---

## Compact comparison table

| Case | Dim. | Regime | Solid kinematics | Main feature tested | Ref. cat. | Mesh study | Δt study | Coupling study | Observed orders | Best agreement with reference | Recommendation |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ringAddedMass | 2-D plane strain | transient free vibration | small (linear geometry) | added mass, μ = 0.1–10 | **A** | 3 levels × 4 regimes | 3 levels | IQN-ILS vs Robin | Δt: 1.90–2.14; ratio (mesh): 1.68–1.77 | ω within 0.13% (−0.12% is BDF2 phase error); Richardson-in-time within +0.018% | **Main, 4.1** |
| womersleyTube | axisymmetric (1° wedge) | time-periodic | small (linear geometry) | viscous wave propagation in compliant tube | **A** | 3 levels | 3 levels | IQN-ILS vs Robin | mesh ≈2 (profile, Q, wall amp, c); ≈1 (wall phase, attenuation); Δt 1.3–1.8 | ≤0.12% on finest except attenuation −0.45% | **Main, 4.1** (after time-order diagnosis) |
| collapsibleChannel | 2-D | transient (damped) | large (TL, StVK) | thin bending wall, very light, strong added mass | **B** (oomph-lib) | solid 4 + fluid 4 levels | 4 levels | IQN-ILS vs Robin; density sweep | fluid mesh 1.7, 2.4 then floor; Δt 1.8, 2.0 | 0.50% of peak (reference precision ≈0.3%) | **Main, 4.2** |
| HronTurek FSI1 | 2-D | steady | large (StVK) | steady FSI, small loads | **C → B-grade** (Featflow 6 levels) | 2 levels | n/a (pseudo-transient) | steady IQN vs Robin | ≈1.5–1.7 (vs Featflow 7+0) | drag 0.24%, u_y 4.6% on 2x | **Main, 4.2 or 4.3** (short) |
| HronTurek FSI3 | 2-D | periodic (self-excited) | large (StVK) | low mass ratio, vortex-induced oscillation | **C** (Featflow L4) | 2 levels (Δt scaled) | none separate | IQN-ILS vs Robin (transient) | 2-level vs Featflow: 1.1–3.8 *(derived here)* | ≤3.1% except lift amp 13.3% | **Main, 4.3** (needs 4x) |
| HronTurek FSI2 | 2-D | periodic (self-excited) | large (StVK) | large-amplitude flapping, mesh motion | **C** (Featflow L4) | 2 levels (Δt scaled) | none separate | tolerance 1e-5/1e-6/1e-7 | 2-level: errors fall 4–10× except lift amp | ≤3.6% except lift amp 7.8% | **Main, 4.3** (with FSI3) |
| 3dTube | **3-D** | transient pulse | small (Hookean) | 3-D pressure-wave propagation, added mass | **C-weak** (digitised figures) | 2 levels | 4 levels | Robin vs IQN-ILS | not defined (2 mesh levels); Δt changes 1.74→0.83→0.016% | u_r,max 1.8% from Tuković 2018 | **Main, 4.3** (needs 3rd level) |
| mokLidDrivenCavity | 2-D | periodic (forced) | large (TL), membrane | enclosed-ish cavity, tension-dominated membrane | **C-weak** (two inconsistent groups) | 4 levels | 3 levels | IQN-ILS vs Aitken (iterations only) | peak ≈1.2; trough <0.5 | trough, mean in group-B envelope; peak 3.3% above | **Main, 4.3** ("revisited") |
| beamInCrossFlow (original) | **3-D** (half) | steady | StVK (small strains) | 3-D external flow, fluid-force resolution | **C** (Richter 2015) | 4 levels | none | Robin vs IQN-ILS (≤0.4%) | not asymptotic | u_x(A) 0.07%, F_x −1.4% vs Richter | **Supplementary** |
| beamInCrossFlow (modified) | **3-D** (half) | steady | **large** (StVK) | 3-D large deformation + mesh motion | **C-weak** (Tuković 2018 Fig. 28, same code lineage) | 5 levels | none | Robin vs IQN-ILS (≤0.4%) | not asymptotic (finest change ≈4.6%) | u_x −2.3%, u_y −0.7%, 2u_z +2.3% | **Supplementary** (main only if converged + independent ref.) |
| cavityFlexibleBottom | 2-D | steady (pseudo-transient) | large (TL, StVK) | internal flow over flexible cavity floor | **C-weak** (Tuković 2018, digitised) | 3 levels | pseudo-Δt (paper data, monolithic) | none recorded | u_y net 1.96 | F_y ≤0.2%; u_y converges to ≈5% offset | **Exclude until discrepancy resolved** |
| blobInTreacle | 2-D | transient → steady | large (TL, StVK) | viscous, strongly damped; 2nd-order time accuracy | **C-weak** (preprint figures, coarse) | 4 levels | 4 levels | none | Δt 2.39–2.54; mesh 2.40 | t = 1 s shape 1.6%; steady u_x −7.4% | **Supplementary** |

Immersed-boundary FSI2, membraneRoof, flexibleDamBreak, fillingElasticContainer,
perpendicularFlap, oneWayCavity and cerebralAneurysm were reviewed and
**excluded** (see [Cases reviewed and excluded](#cases-reviewed-and-excluded)).

---

## Detailed notes

### 1. ringAddedMass

- **Files:** `tutorials/fluidSolidInteraction/ringAddedMass/{README.md, verification/README.md, verification/reference/ringAddedMass_verification_references.json, verification/scripts/ring_exact_frequencies.py}`; PR #510.
- **Dimension / regime:** 2-D plane strain, quarter model with symmetry; transient free vibration (n = 2 ovalling mode) over four periods.
- **Fluid:** `pimpleFluid`, incompressible Navier-Stokes with ν_f = 1e-7 m²/s, so effectively inviscid (Stokes layer < 0.12% of gap, unresolved by design).
- **Solid:** `linearGeometryTotalDisplacement`, PETSc SNES, cubic high-order residual (default) or standard second-order residual. Small deformation (amplitude ≈ 5e-4 a).
- **Feature tested:** partitioned coupling under controllable added mass. The mass ratio is μ = 0.101, 1.008 and 10.08 (weak, moderate, strong) at fixed solid. This is a clean isolation of the added-mass mechanism.
- **Reference (cat. A):** exact linear eigenvalue of a plane-strain continuum ring with exact potential-flow added-mass traction, by Chebyshev collocation (24 intervals, frequencies to ≈1e-8; `collocationChange` ≈ 4.6e-8 in JSON). Thin-ring closed form (Brennen 1982 added mass) checks it to O((h/a)²). Linearised viscous solution also stored: frequency shifts −0.004/−0.026/−0.073% and damping ratios of the same size.
- **Reference values:** ω_dry = 2.555042; ω_wet = 2.435615 / 1.804513 / 0.768654 rad/s; ratios 0.953258 / 0.706256 / 0.300838.
- **QoIs:** wet and dry frequency (zero crossings), wet/dry frequency ratio, amplitude loss per period.
- **Mesh study:** refinement factors 1, 2, 4 (solid 4×48 → 16×192; fluid 12×48 → 48×192), matching interfaces, 100 steps/period. High-order solid frequency changes ≤0.06% between meshes ("converged in space on the coarsest mesh"). The ratio error converges at 1.77/1.73/1.68 to +0.001/+0.003/+0.006% on the finest mesh.
- **Time-step study:** 25, 50, 100 steps/period on mesh 2. Observed orders 1.90 (dry), 1.92, 2.02, 2.14. Richardson extrapolation gives every frequency within +0.018% of exact.
- **Solid-discretisation study:** standard solid is +20.66/+6.17/+1.58% too stiff (dry), order 1.66. Its wet/dry ratio still agrees with the high-order one to 0.003% (coupling verified independently of solid accuracy).
- **Coupling study:** IQN-ILS vs Robin-Neumann on mesh 2 agree to 0.03% in frequency. Mean iterations: IQN-ILS 2.99/3.95/4.91 vs Robin 3.87/4.85/12.29 (weak/moderate/strong).
- **Best values (finest, 100 steps/period, IQN-ILS):** ω = 2.55173 (dry), 2.43248, 1.80224, 0.76770 rad/s; errors −0.13% (BDF2 phase error).
- **Independent-code evidence:** none needed (exact). The collocation script is independent of solids4foam but by the same authors.
- **Cost:** tutorial ≈75 s serial; full study 35 runs ≈30 min with 32 concurrent serial runs; longest run ≈25 min.
- **Limitations:**
  - Viscous damping not verified (layer unresolved by design).
  - IQN-ILS shows mesh-dependent amplitude loss (up to 2.1%/period, coarsest, strong). It halves with each refinement and is absent with Robin; the mechanism is unexplained.
  - Recorded on v2412 only. The v2412 regression references fail on v2512 by 0.05% in half-period while passing against the exact value (PR #510 caveats).
  - Frequency mesh order not measurable (already converged).
- **Unique contribution:** the only case with an exact reference *and* a swept added-mass ratio over two decades; separates coupling, solid and time errors cleanly.
- **Would strengthen it:** v2512 rerun; explanation of IQN-ILS amplitude decay; optionally n = 3 and a thinner ring (h/a) to show the O((h/a)²) thin-ring trend; monolithic solver (PR #233) on the strong level.

### 2. womersleyTube

- **Files:** `womersleyTube/{README.md, verification/README.md, verification/reference/womersleyTube_verification_references.json, verification/scripts/womersley_exact.py}`; PR #513.
- **Dimension / regime:** axisymmetric (1° wedge); time-periodic travelling wave (6 periods run, last 2 analysed).
- **Fluid:** `pimpleFluid`, viscous incompressible, Womersley number α = 5.01.
- **Solid:** `linearGeometryTotalDisplacement`, PETSc SNES, standard second-order residual; the high-order residual fails on wedges (#512). Small deformation (5e-4 R). Thick wall h/R = 0.1, free (untethered) wall.
- **Feature tested:** viscous pressure-wave propagation in a compliant tube; non-reflecting finite segment via exact travelling-wave boundary data at both ends; wall compliance plus viscous attenuation.
- **Reference (cat. A):** exact linear solution for a thick elastic cylinder (Chebyshev collocation, 16 intervals) coupled to full linearised Navier-Stokes. Collocation checks: k8 vs k24 agree to ~1e-11. Checked against Womersley's thin-wall equation in the h/R → 0 limit (Filonova et al. 2020 Eq. 41) and against the Lamé-stiffness tethered long-wave limit to 2e-5. Womersley's classical thin-wall speed is 4.7% too high here.
- **Reference values:** k = 0.143978 − 0.020727i m⁻¹; c = 0.872796 m/s; λ = 43.64 m.
- **QoIs:** velocity-profile max error, flow-rate amplitude/phase at L/2, wall radial displacement amplitude/phase, wave speed and attenuation from a log-pressure fit along the tube.
- **Mesh study:** factors 1, 2, 4 (base 32 axial, 10 fluid radial, 4 wall radial) at 200 steps/period, Robin coupling. Orders: profile 2.06, flow_amp 2.16, flow_phase 1.92, wallMid_amp 2.15, speed 2.48, wallMid_phase 0.97, attenuation 0.63.
- **Time-step study:** 50, 100, 200 steps/period on mesh 2. Orders 1.29–1.84 (flow_amp 1.64, flow_phase 1.84, wallMid_amp 1.71, wallMid_phase 1.74, speed 1.29, attenuation 1.51), **below the nominal 2; cause not identified**.
- **Coupling study:** IQN-ILS vs Robin agree to 3e-4 on every quantity; 14.8 vs 8.7 iterations/step.
- **Best values (m4, n200):** profile 1.18e-3; flow_amp −7.15e-4; flow_phase 1.5e-5 rad; wallMid_amp 8.6e-4; wallMid_phase −1.88e-3 rad; speed 8.9e-4; attenuation −4.47e-3.
- **Independent-code evidence:** none needed (exact).
- **Cost:** tutorial ≈4 min serial (IQN-ILS), half with Robin; runs 2.5–25 min; whole study ≈30 min with 5 cores.
- **Limitations / unresolved:**
  - Time order 1.3–1.8.
  - An unexplained ≈0.5% start-up transient that grows as Δt shrinks.
  - The mesh study at n200 shares a time error of the same size as the finest mesh error.
  - Wall phase and attenuation converge at about first order in space.
  - Untethered wall differs from Womersley's classical constrained tube (#511).
  - Robin result moves 0.26% of amplitude between v2412 and v2512 (≈10× IQN-ILS sensitivity, unexplained; PR #513 comment).
  - OpenFOAM.com only; recorded on v2412 macOS.
  - *Checked here:* the BDF start-up fix `fdf9b856` (branch `feature/bdf34-time-schemes`, PR #503 line) only changes the new `bdfD2dt2Scheme`. Its commit message states `backward` already used the correct start-up, so it is not the explanation.
- **Unique contribution:** exact, viscous, time-periodic, wave-propagation reference with distinct amplitude and phase QoIs. Complements ringAddedMass (inviscid, standing mode). The only axisymmetric case.
- **Would strengthen it:**
  - **Diagnose the time order (essential):** time study on mesh 4 and/or n = 400, interface-velocity condition order, start-up transient.
  - v2512 rerun.
  - Tethered variant once #511 is fixed.

### 3. collapsibleChannel

- **Files:** `collapsibleChannel/{README.md, verification/README.md, verification/reference/*.csv, collapsibleChannel_verification_references.json, reference/oomph-lib/collapsible_channel_ref.cc}`; PR #507.
- **Dimension / regime:** 2-D; transient (ramped external pressure, collapse, damped oscillation, period ≈0.8 s), analysed to t = 1.5 s.
- **Fluid:** `pimpleFluid`, Re = Re·St = 50, Poiseuille inflow driven by pressure.
- **Solid:** total-Lagrangian St Venant-Kirchhoff, PETSc SNES (`snes_mf_operator`), cubic high-order (default) or linear least-squares reconstruction. **Large deformation** (collapse ≈ 21% of channel width at first trough, h/L = 0.005). Regularising density ρ_s = ρ_f = 1 kg/m³.
- **Feature tested:** extremely light, thin, bending wall; very strong added mass; internal flow. Partitioned robustness limits are mapped by a density sweep.
- **Reference (cat. B):** oomph-lib monolithic Newton FE, Q2P-1 fluid plus geometrically nonlinear Kirchhoff-Love beam, BDF2/Newmark. The **driver source is stored in the repo** and regenerable. The converged reference is Richardson-in-time (Δt = 0.003125, 0.0015625).
  - Reference precision ≈0.3% of peak, dominated by oomph-lib spatial error: refinement 1 vs 2 = 0.26%; Richardson correction 0.06%; CR vs Taylor-Hood < 0.001%.
  - The reference **modifies** the published Jensen & Heil / Heil problem: σ₀ = 0 (no pre-stress), h/a = 1/20, clamped ends, pressure ramp, regularising density (shift 0.28% of peak).
  - Model-form difference beam vs continuum is O(h/L) (static: 0.05%).
- **QoIs:** wall-midpoint (and quarter points) vertical displacement history; first trough; static deflection.
- **Studies and results:**
  - **Static wall:** cubic within 0.25% on all meshes, converging to 0.05% (continuum-beam difference). Linear solid −16.7% → −0.22% (≈2nd order in in-plane spacing).
  - **Solid-in-FSI:** cubic 0.44% solid error on the tutorial mesh vs 17.4% for linear.
  - **Fluid mesh:** 4 levels, Robin, Δt = 0.00625, vs oomph-lib at the same Δt. Error 4.65 → 1.43 → 0.53 → 0.50%; orders 1.7, 2.4, then the reference floor.
  - **Time step:** Δt 0.025 → 0.003125 on the tutorial mesh. Self-differences 4.56, 1.33, 0.34%; orders 1.8, 2.0 (oomph-lib itself 1.8, 1.9).
- **Coupling / iteration study:**
  - Density sweep ρ_s/ρ_f = 0…3 maps IQN-ILS and Robin failure regions: IQN-ILS never completes the whole study; Robin completes from ρ_s/ρ_f = 1.
  - Robin vs IQN-ILS agree to 0.08–0.23% of peak, but **the gap does not shrink with refinement** (≈ coupling tolerance level).
- **Best values:** first trough −0.21610 m (fluid level 4, Δt 0.00625) vs converged oomph-lib −0.2160 m.
- **Independent-code evidence:** **yes**: oomph-lib (independent FE, monolithic). Strongest independent-code comparison in the suite.
- **Cost:** tutorial ≈10 min (cubic, IQN-ILS); full sweep ≈2.5 h one case per core; finest fluid mesh 70 min. **Serial only** (high-order Jacobian alpha scheme not parallelised).
- **Limitations:**
  - Not the published (pre-stressed, massless) problem.
  - The temporal error on the finest mesh was not measured separately.
  - Mesh and time studies are separate refinements (no combined path).
  - Parallel solid fails; a 4-rank linear-reconstruction run gave a *wrong* displacement, to recheck after the fix.
  - Tested on v2412 only; the verification refuses other variants.
- **Unique contribution:** independently converged monolithic reference for a large-deformation, bending-dominated, strongly added-mass internal flow. Also a documented map of where partitioned couplings fail.
- **Would strengthen it:** time study on fluid level ≥3; repeat after the parallel Jacobian fix; monolithic solver (PR #233) run with ρ_s = 0 (the published massless case) as a robustness contrast.

### 4. HronTurek (FSI1, FSI2, FSI3)

- **Files:** `HronTurek/{README.md, verification/README.md, verification/reference/HronTurek_verification_references.json, TurekHron_fsi2/3_reference_history.csv}`; PRs #506, #472, #492, #219.
- **Dimension / regime:** 2-D (one spanwise cell). FSI1 steady (pseudo-transient to t = 30 s); FSI2 and FSI3 periodic self-excited.
- **Fluid:** `pimpleFluid`; Re = 20 (FSI1), 100 (FSI2), 200 (FSI3).
- **Solid:** St Venant-Kirchhoff (`StVenantKirchhoffElastic`), nonlinear geometry; **large deformation** in FSI2 (u_y ≈ ±80 mm).
- **Features tested:** FSI1 steady coupling with small loads; FSI3 low mass ratio (ρ_s = ρ_f), high frequency; FSI2 heavy, soft plate, large-amplitude motion and fluid-mesh distortion.
- **Reference:**
  - FSI1: Featflow table, six levels 2+0 → 7+0 (≈1M elements), spread ≤0.73% (u_x), 0.19% (u_y), 0.14% (drag). Effectively an **independently converged (B-grade)** reference.
  - FSI2 / FSI3: Featflow level 4 tables (FSI3 Δt = 0.00025; FSI2 Δt = 0.0005). FSI2 level 3→4 changes 0.9–2.5%, so the reference is uncertain by ≈1%. Category **C**.
  - **Ambiguity flagged:** the commonly quoted Turek & Hron (2006) summary values differ from the Featflow tables. FSI3 drag amplitude 22.66 vs 27.74 N/m (−18%), u_y frequency 5.3 vs 5.46 Hz; FSI2 u_y frequency 2.0 vs 1.93 Hz. The study uses the tables; the driver reproduces them from the published histories to ≤0.7% (FSI3, drag amplitude 1.1%) and ≤0.3% (FSI2).
- **QoIs:** tip point A u_x, u_y (mean ± amplitude, frequency); drag and lift on cylinder + plate per unit depth.
- **Mesh studies** (levels 1x, 2x; Δt halved with mesh; 4x defined but not run):

  | FSI3 (IQN-ILS) | 1x error *(derived here)* | 2x error | 2-level order vs Featflow *(derived here)* |
  |---|---:|---:|---:|
  | u_x mean | 21.3% | 3.1% | 2.8 |
  | u_x amplitude | 19.2% | 2.1% | 3.2 |
  | u_y amplitude | 14.4% | 2.3% | 2.6 |
  | u_y frequency | 2.5% | 1.2% | 1.1 |
  | drag mean | 0.7% | 0.3% | 1.5 |
  | drag amplitude | 15.5% | 1.2% | 3.8 |
  | lift amplitude | 48.9% | 13.3% | 1.9 |

  - FSI2 2x errors: u_x mean 3.1%, u_x amp 1.9%, u_y amp 1.3%, frequencies 1.8%, drag mean 0.5%, drag amp 2.0%, lift amp 7.8% (same 7.9% on 1x, i.e. not converging; 2.4% of it is coupling noise).
  - FSI1 2x errors: u_x 0.97%, u_y 4.64%, drag 0.24%, lift 1.55%; observed orders 1.5–1.7 against Featflow 7+0.
- **Time-step study:** none independent of the mesh (Δt scaled with h). FSI1: Δt 0.025 → 0.01 changes u_y 0.04%, u_x 0.26%.
- **Coupling / tolerance studies:**
  - FSI3 transient Robin vs IQN-ILS from a common restart: the difference decreases 1x→2x and is smaller than the 1x→2x change. It is attributed to a different discrete interface condition (Robin pressure gradient from solid acceleration), not to coupling error. Robin needs 152 (1x) and 52 (2x) iterations/step vs 9 for IQN-ILS (#493).
  - FSI1 steady Robin vs IQN-ILS agree to ≤0.01%.
  - FSI2 `outerCorrTolerance` sweep: 1e-5 inflates lift amplitude by 13% (1x) and 30% (2x); 1e-6 vs 1e-7 shows no change. **The only explicit coupling-tolerance study in the suite.**
  - FSI1 fluid-solver tolerance 1e-6 → 1e-9 changes steady values ≤2e-4.
- **Cost:** FSI3 1x 1.1 h serial, 2x 4.1 h on 8 ranks; FSI2 1.2 h / 4.0 h (7 ranks); FSI1 781 s / 3062 s (2 ranks); FSI1 4x estimated >14 h on 3 ranks.
- **Independent-code evidence:** Featflow (independent monolithic FE). Tuković et al. (2018) values stored for information (same code lineage).
- **Limitations:**
  - Only two mesh levels, so no reference-free order.
  - FSI3 lift amplitude 13% off.
  - FSI2 lift amplitude stuck at 8%.
  - FSI2 fluid mesh degrades (non-orthogonality 26° → 58° by 10.5 s; IQN-ILS failure at 11.4 s on 2x), so the run is truncated at 10.5 s.
  - Body-fitted FSI2 diverges at ≈6.5 s in some variants per PR #491 (#489, since addressed by #492 predictor default; verify).
  - Recorded on v2412 (M1 Ultra, shared machine).
- **Unique contribution:** the community-standard benchmark; three regimes (steady, low-mass periodic, large-amplitude periodic) with a single geometry; the only tolerance-sensitivity evidence.
- **Would strengthen it:** the 4x level for FSI3 (essential for the main paper) and FSI1 (desirable); FSI2 4x after a mesh-motion fix; an independent time-step study at 2x for FSI3.

### 5. 3dTube

- **Files:** `3dTube/{README.md, verification/README.md, verification/reference/{3dTube_verification_references.json, Tukovic2018_fig25a_pointA.csv, Lozovskiy2019_fig3_pointA.csv, Eken2016_fig5.6_pointA.csv}}`; PR #501.
- **Dimension / regime:** **genuinely 3-D** (quarter tube); transient pressure pulse (3 ms), run to 20 ms including the outlet reflection.
- **Fluid:** `pimpleFluid`, ρ = 1000, ν = 3e-6.
- **Solid:** Hookean small-strain linear elastic (literature mostly StVK; peak hoop strain ≈3%). Small-deformation formulation.
- **Feature tested:** 3-D pressure-wave propagation in a compliant vessel with strong added mass (ρ_s/ρ_f = 1.2); the "hemodynamic" partitioned-coupling benchmark of Fernández & Moubachir / Gerbeau & Vidrascu.
- **Reference (cat. C-weak):**
  - No tabulated published values.
  - Tuković et al. (2018) Fig. 25 history (vector-extracted; FV, backward, same code lineage).
  - Lozovskiy et al. (2019) and Eken (2016) histories (raster-digitised ±3e-6 m; independent monolithic FE but implicit Euler Δt = 1e-4, which damps the peak by ≈23%).
  - Analytical wave speeds: Moens-Korteweg 5.2–5.7 m/s (upper bounds); thick-wall with axial-stress correction 4.81 m/s (Tuković 2018 Eq. 36), approximate.
- **QoIs:** u_r,max(A) and its time; u_z,min(A) incident trough (t < 6.5 ms); u_r,min after reflection; front arrival time; pulse-wave speed c_p.
- **Mesh study:** levels 1 (16k fluid / 6.4k solid) and 2 (128k / 51k), Δt halved with mesh. Changes: u_r,max 0.05%, u_z,min 1.85%, t_arr 0.34%, c_p 1.58%. Level 3 (≈1.4M cells) not run.
- **Time-step study:** Δt 1e-4 → 1.25e-5 on level 1. Successive u_r,max changes 1.74, 0.83, 0.016%; t_arr 2.23, 0.44, 0.14%. The tutorial Δt is converged to ≈0.1%.
- **Literature-discretisation study:** Euler, Δt = 1e-4 reproduces the FE peaks (+2.5% vs Lozovskiy, −1.4% vs Eken). This **explains the ≈25% split between FE and FV literature** as a time-discretisation effect, a genuine "revisited" result.
- **Coupling study:** Robin vs IQN-ILS histories differ ≤0.27% of peak, QoIs ≤0.08%. Robin 4.74 vs IQN-ILS 15.53 iterations/step, 4.1× faster. Predictor on/off changes ≤0.001%.
- **Best values (level 2):** u_r,max = 0.15976 mm (1.8% above Tuković 0.15694 mm); c_p = 4.695 m/s (2.4% below 4.81 m/s estimate).
- **Cost:** L1 425 s (4 ranks); L2 7197 s; L3 estimated >10 h on 6 ranks.
- **Limitations:**
  - Two mesh levels only (no order).
  - u_z,min still changing 1.85%.
  - References read from figures.
  - Hookean vs StVK in most literature.
  - Pressure ramp-off 3.0–3.1 ms vs instantaneous.
  - Robin iteration count rises to 7.4 at the smallest Δt (uninvestigated).
- **Unique contribution:** the only genuinely 3-D transient case with a completed mesh, time and coupling study; resolves a literature discrepancy.
- **Would strengthen it:** level 3 (essential for main-paper use); an independent converged reference (COMSOL requested, PR #501 comment).

### 6. mokLidDrivenCavity

- **Files:** `mokLidDrivenCavity/{README.md, verification/README.md, verification/reference/*.csv, mokLidDrivenCavity_verification_references.json}`; PR #497.
- **Dimension / regime:** 2-D; periodic forced (lid velocity 1 − cos(2πt/5)), window t = 20–70 s.
- **Fluid:** `pimpleFluid`, ρ = 1, ν = 0.01.
- **Solid:** membrane h = 0.002 m, E = 250 Pa; high-order TL (default) or standard solid. Large displacement (≈0.29 m on a 1 m cavity), tension-dominated.
- **Feature tested:** classic partitioned-coupling benchmark (Wall 1999; Mok 2001); tension-dominated membrane; near-enclosed flow.
- **Reference (cat. C-weak, ambiguous):**
  - Published curves fall into **two inconsistent groups**. Group A (Mok, Wall, Gerbeau-Vidrascu) has mesh-dependent, unstated openings and is internally inconsistent (peaks 0.177–0.215 m).
  - Group B (Valdés 2007, Kratos example, Tiba et al. 2026) is fully specified. **Tiba et al. 2026 is an arXiv preprint (arXiv:2609.16876)** giving scalars only.
  - No published source reports a mesh or time-step study.
  - Within group B: peak spread 5.7%, amplitude spread 20%.
- **QoIs:** membrane midpoint displacement: mean peak, trough, time-mean over periodic window; history distance to the Valdés–Kratos band.
- **Mesh study:** 16², 32², 64², 96² with Δt scaled. Successive changes 1.56, 1.26, 0.78% of peak; peak order ≈1.2 (→ ≈0.285 m); trough order < 0.5 (uncertain 1–2%). 128² diverged (cause not isolated).
- **Time-step study:** Δt 0.1/0.05/0.025 on 64²: changes 0.17, 0.04% of peak.
- **Solid study:** standard (2 or 8 cells through thickness) vs high-order agree to 0.04%. High-order gives no benefit (tension-dominated) and costs 1.5×.
- **Coupling:** IQN-ILS 2.4–4.8 iterations/step; Aitken "same response" (no number recorded) at 9.9 iterations/step; literature iteration counts compared indicatively.
- **Best values (96²):** peak 0.2861 m (+3.32% above group-B envelope max), trough 0.2044 m (+0.56%), mean 0.2445 m (inside).
- **Independent-code evidence:** Valdés (FE), Kratos (FE), Kassiotis et al. (OpenFOAM-based), all single mesh/Δt. A commercial-code solution has been requested (not yet available).
- **Cost:** 96² mesh 3729 s on 2 ranks; quick run ≈7 min.
- **Limitations:**
  - Trough not mesh-converged to better than ≈1–2%.
  - Finest level is 96² (128² diverges).
  - Only one vector-extracted group-B history.
  - foam-extend differs by 3% in midpoint displacement at t = 5 s.
  - OpenFOAM-9 regression values differ ≈1% between machines (PR #497 comment).
- **Unique contribution:** an "anatomy of a benchmark" case: shows the literature definitions disagree and that solids4foam is converged where references are not. Also the tension-dominated counterpoint to collapsibleChannel's bending wall.
- **Would strengthen it:** independent converged reference (COMSOL/other); diagnosis of the 128² divergence; Richardson estimate for peak and trough.

### 7. beamInCrossFlow (original small-deformation and modified large-deformation forms)

- **Files:** `beamInCrossFlow/{README.md, verification/README.md, verification/reference/{beamInCrossFlow_verification_references.json, original_mesh_t8_backward_vs_references.png, modified_mesh_t8_backward_vs_references.png}}`; PRs #412, #518, #450.
- **Dimension / regime:** **genuinely 3-D** (half-domain, z = 0 symmetry); steady state reached by transient run to t = 8 s (backward).
- **Fluid:** `pimpleFluid`, ρ = 1000, ν = 1e-3, parabolic inflow.
- **Solid:** St Venant-Kirchhoff, nonlinear geometry, both forms.
  - Original: E = 1.4 MPa, U_max = 0.2 m/s (Re = 40), small strain.
  - Modified: E = 10 kPa, U_max = 0.3 m/s, 1 s ramp, **large deformation**.
- **Feature tested:** 3-D external flow around a thick plate; resolution of fluid forces; (modified) large solid displacement plus 3-D mesh motion.
- **Reference:**
  - Original (cat. C): Richter (2015) u_x(A) = 5.95e-5 m, F_x = 1.33 N (GASCOIGNE 3D monolithic FE, **independent**). Which Richter refinement level and whether these are extrapolated values still needs checking against the paper. Tuković et al. (2018) OpenFOAM values u_x = 5.93e-5, u_y = 2.40e-5, F_x = 1.31, F_y = 0.1055 are diagnostics (same lineage).
  - Modified (cat. **C-weak**): Tuković et al. (2018) Fig. 28 only, u = (0.01463, 0.005, −0.000447) m. Same code lineage, single unstructured mesh (273,539 fluid / 6,661 solid cells), **not independent**. Gillebaart et al. (2016) also study a 3-D elastic beam in a channel; whether they tabulate values for this form is unchecked.
- **QoIs:** u_x, u_y, u_z at point A; total interface force F_x, F_y.
- **Mesh study:**
  - Original: 1x/2x/4x/8x uniform (≈1x to 512x cells; 8x ≈7.6M cells), Δt scaled.
  - Modified: 1x/2x/3x/4x/8x.
  - Only finest values are tabulated (PR #412). Original 8x: u_x(A) = 5.95387e-5 m (+0.07% vs Richter), F_x = 1.310774 N (−1.4% vs Richter, +0.06% vs Tuković). Modified 8x: u_x = 0.0143007 (−2.3%), u_y = 0.00496486 (−0.7%), 2u_z = −0.00045743 (+2.3%) *(errors derived here)*.
  - *(read from plot)*: original u_x ≈ 3.93e-5, 5.04e-5, 5.69e-5, 5.95e-5 m (Δx = 0.025 → 0.003125), so the last change ≈4.6%. Modified u_x ≈ 0.00978, 0.01214, 0.01315 (3x), 0.01365, 0.01430, so the last change ≈4.7%. F_y (original) non-monotone.
  - **Not in the asymptotic range;** the finest-level agreement with the references cannot yet be separated from coincidence.
- **Time-step study:** none (Δt scaled with mesh; `--time-scheme Euler` diagnostic only).
- **Coupling study:** Robin vs IQN-ILS on base mesh within 0.4% for both forms (PR #518). PR #518 also has an IQN-ILS option sweep (`relMinSignificant 1e-2` fixes a near-steady stall); two-material variant agrees across couplings to 0.06–0.1%.
- **Cost:** original 8x ≈11 h on 64 cores; modified 8x on 64 ranks (time not recorded).
- **Limitations:**
  - No Richardson / order evidence.
  - Very fine fluid mesh needed near the block (PR #412 comment).
  - Modified-form reference is not independent.
  - PR #518 notes the refined verification meshes were not re-run after the IQN-ILS default change (base mesh changed 0.01%).
- **Unique contribution:** the only 3-D steady case and the only 3-D large-deformation case; both forms share geometry, so they separate geometric nonlinearity from discretisation.
- **Would strengthen it:** a graded (not uniform) fluid-mesh family reaching the asymptotic range (or a 16x level); tabulated per-level values; Richter level check; an independent large-deformation reference (Gillebaart 2016 or COMSOL).

### 8. cavityFlexibleBottom

- **Files:** `cavityFlexibleBottom/{README.md, verification/README.md, verification/reference/{TukovicDisplacements.csv, TukovicForces.csv, *.json}}`; PR #439; **this paper repo:** `simulationData/cavityFlexibleBottom/campaignSummary.tsv`, `sectionTestCases.tex` (current Case 1).
- **Dimension / regime:** 2-D; steady, reached by pseudo-transient (Euler, solid `dampingCoeff`).
- **Fluid:** `pimpleFluid`, Re = 100.
- **Solid:** `nonLinearGeometryTotalLagrangianTotalDisplacement`, StVK, E = 500 Pa; large deformation (u_y ≈ −0.24 m).
- **Feature tested:** internal flow over a flexible cavity floor; steady FSI.
- **Reference (cat. C-weak):** Tuković et al. (2018) mesh-study curves, digitised (u_y −0.2333 → −0.2506 m over Δx 0.1 → 0.0125; F_y per 0.05 m depth). **Same code lineage; no independent code.**
- **QoIs:** steady u_y at (4, −1, 0.5); steady interface F_y.
- **Mesh study (partitioned, Aitken, v2512):** Δx 0.1/0.05/0.025. u_y −0.20246, −0.22883, −0.23563 (errors 13.2, 6.88, 5.59%); net order 1.96 **towards a value ≈5% from the published curve**. F_y within 0.2% on all levels. Level 4 diverges at the tutorial Δt.
- **Pseudo-time-step study (monolithic, paper repo data):** Δt 0.5–40 s on 3 meshes; u_y and F_y insensitive to Δt (e.g. m3: −0.29192 to −0.29193).
- **Coupling study:** none recorded (IQN-ILS variant exists).
- **⚠ Unresolved discrepancy (found in this inventory):**
  - At Δx = 0.025 the monolithic campaign in `simulationData/` gives u_y = −0.29192 m and F_y = −5.25227 N.
  - The partitioned verification gives u_y = −0.23563 m and F_y = −5.25632 N.
  - The **forces agree to 0.08% but the displacements differ by 24%** *(derived here)*, and the published curve (≈ −0.2496 at 0.025) lies between them.
  - Equal loads with different steady deflection point to a difference in solid model, material, monitor point or boundary conditions between the two set-ups, not to the coupling.
  - The tutorial README also quotes an earlier solids4foam coarse-mesh value of −0.219 m vs −0.20246 m now.
- **Cost:** default partitioned sweep ≈10 min serial; monolithic m3 runs 2,000–16,000 s.
- **Unique contribution:** steady internal flow; currently the parent paper's Case 1 and the only case with a monolithic-solver dataset.
- **Recommendation:** exclude from the verification suite until the discrepancy is resolved. Even then the reference is same-lineage and digitised, so at most supplementary.

### 9. blobInTreacle

- **Files:** `blobInTreacle/{README.md, verification/README.md, verification/reference/{LiuInterface_t1.csv, LiuInterface_t10.csv, blobInTreacle_verification_references.json}}`; PR #498. Also the parent paper's current Case 6.
- **Dimension / regime:** 2-D; transient ramp (to t = 2 s) then steady.
- **Fluid:** `pimpleFluid` with `consistent yes` (SIMPLEC), ρ = 1, ν = 1 (highly viscous).
- **Solid:** St Venant-Kirchhoff (preprint describes linear), large deformation (u_x ≈ 0.12 m on r = 0.1 m half-cylinder).
- **Feature tested:** temporal order of accuracy of the coupled scheme (the purpose of the Liu et al. 2014 case); viscous, strongly damped.
- **Reference (cat. C-weak):**
  - Liu et al. (2014, JCP) report only self-convergence norms, no solution.
  - The preprint (arXiv:1401.0082) figures were vector-extracted. At t = 1 s: P5/P4/P5 elements, first order, Δt = 0.025, ≈4% time error of its own. Steady state: coarse mesh (1,239 fluid / 370 solid triangles). Top-point steady displacement (0.1297, −0.0117) m.
- **QoIs:** top-point u_x, u_y; interface shape at t = 1 s and steady.
- **Time-step study:** Δt 0.25 → 0.03125 at t = 1 s. Orders 2.39, 2.54 (u_x) and 2.39, 2.41 (interface). Below Δt = 0.03125 changes reach coupling/solver tolerance (≈1e-5 m).
- **Mesh study:** x0.5–x4. Steady u_x order 2.40 → extrapolated 0.12026 m; finest 0.1% from it.
- **Best values:** steady u_x = 0.120152 m (−7.4% vs published); u_y −0.013996 m (20% larger magnitude than published); t = 1 s shape within 1.6% (inside the reference's 4%).
- **Independent-code evidence:** Liu (independent FE), but coarse and figure-based. COMSOL being arranged.
- **Cost:** whole study ≈33 min (x2, x4 on 4 ranks).
- **Limitations:**
  - Steady-state offset (6–7%) attributed to the coarse reference but unconfirmed.
  - Time study only at t = 1 s during the ramp, with large Δt (the tolerance floor limits smaller Δt).
  - Finest mesh needs Δt = 0.005 (IQN-ILS stagnation).
  - Not cross-checked against the quasi-monolithic solver after the StVK change (PR #498).
- **Unique contribution:** clean demonstration of second-order temporal accuracy in a viscous, large-deformation FSI problem.
- **Would strengthen it:** independent converged reference; time study with tightened coupling/solver tolerances to extend below Δt = 0.03125.

---

## Cases reviewed and excluded

| Case | Location / PR | Reason for exclusion |
|---|---|---|
| immersedHronTurekFsi2 | `tutorials/fluids/immersedBoundary/`, PR #491 | Immersed-boundary method (not body-fitted ALE); compared with the 2006 summary values; tip amplitude 7%/5% low, force amplitudes not mesh-converged (6–9% high on finer mesh); 7.8 h on 4 cores. Separate paper. |
| membraneRoof | PR #499 | Self-declared "qualitative tutorial, not a verification"; single-run reference (von Scheven 2009), reference oscillation not reproduced, discrepancy unresolved. |
| flexibleDamBreak | PR #308 | Free-surface two-phase; the two references disagree after t ≈ 0.6 s; no mesh/time-step independence study. |
| fillingElasticContainer | README | Two-phase; README states only qualitative agreement with Cerquaglia et al. (2019) is expected. |
| perpendicularFlap, oneWayCavity | README | Demonstration/preCICE tutorials; no reference data. |
| cerebralAneurysm | README | Patient-specific demonstration; regression values only. |
| Parent-paper Cases 7, 8 (block with flexible tail; 3-D block with flapping plate) | `sectionTestCases.tex` | No solids4foam tutorial or verification evidence found. |

---

## Cross-cutting observations

1. **Two exact references now exist** (ringAddedMass, womersleyTube). Each is
   computed by a stored, independently checkable script with collocation
   convergence to round-off. This is the backbone of a verification paper.
2. **Only two cases have an independent-code reference of the same problem
   with its own convergence evidence:** collapsibleChannel (oomph-lib) and
   HronTurek FSI1 (Featflow, 6 levels). FSI2/FSI3 Featflow references carry
   ≈1% uncertainty.
3. **Same-lineage references:** Tuković et al. (2018) is the reference for
   cavityFlexibleBottom and modified beamInCrossFlow, and the principal
   converged reference for 3dTube. It is an earlier version of the same code
   family and must not be presented as independent.
4. **Coupling-error evidence is uneven.** Every case except cavityFlexibleBottom
   and blobInTreacle compares two converged couplings. Only HronTurek FSI2
   (outerCorrTolerance 1e-5/1e-6/1e-7) and FSI1 (fluid solver tolerance) vary a
   tolerance systematically. collapsibleChannel shows a Robin–IQN-ILS gap of
   0.1–0.2% that does not shrink with refinement.
5. **Version heterogeneity.** Recorded on v2412 (ring, Womersley,
   collapsible, HronTurek, 3dTube, Mok) and v2512 (blob, cavity, beam
   #518), on macOS and Linux. Womersley Robin moves 0.26% of amplitude
   v2412 → v2512; ring regression half-period moves 0.05%. One consistent
   version and commit is needed for the paper.
6. **Sub-nominal temporal orders.** Womersley (1.3–1.8) is the only case where
   the observed time order is clearly below 2 with an exact reference.
   ring, collapsible and blob all show ≈1.8–2.5.

---

## Proposed Section 4 structure (outline only, not written)

A commented LaTeX version is in [`section4_outline.tex`](section4_outline.tex).

**4 Verification benchmarks**

- 4.0 *Verification protocol* (short): QoIs, error measures, observed order
  (successive differences vs. against reference), Richardson extrapolation,
  periodic-window statistics, coupling-tolerance policy, software version and
  reproducibility (`Allverify`, fingerprints, Zenodo archive).

**4.1 Analytical and semi-analytical reference problems**

| Case | Recommendation | Role |
|---|---|---|
| 4.1.1 ringAddedMass | **Main paper** | Exact added mass; μ sweep 0.1–10; separates time, solid and coupling errors; IQN-ILS vs Robin. |
| 4.1.2 womersleyTube | **Main paper** (conditional on diagnosing the time order, else report it as an open finding) | Exact viscous wave propagation; axisymmetric; amplitude and phase QoIs. |

**4.2 Independently converged numerical reference problems**

| Case | Recommendation | Role |
|---|---|---|
| 4.2.1 collapsibleChannel | **Main paper** | oomph-lib monolithic reference regenerated from stored source; large-deformation bending wall; partitioned failure map vs. wall density. |
| 4.2.2 HronTurek FSI1 | **Main paper** (short; could instead open 4.3.1) | Steady benchmark whose Featflow reference is converged to ≈0.2% over six levels. |

**4.3 Classical FSI benchmarks revisited**

| Case | Recommendation | Role |
|---|---|---|
| 4.3.1 HronTurek FSI3 + FSI2 | **Main paper** (one table; per-quantity details in supplement) | Community benchmark; table-vs-summary ambiguity; tolerance sensitivity of lift; FSI2 mesh-motion limit. |
| 4.3.2 3dTube | **Main paper** | Only 3-D transient study; explains the FE/FV literature split as time discretisation. |
| 4.3.3 mokLidDrivenCavity | **Main paper** | Literature definitions inconsistent; solids4foam converged where references are not; tension-dominated counterpoint. |
| 4.3.4 beamInCrossFlow (original + modified) | **Supplementary** (promote modified form to main only if asymptotic convergence and an independent reference are obtained) | 3-D steady; geometric nonlinearity. |
| 4.3.5 blobInTreacle | **Supplementary** | Second-order temporal accuracy; weak reference. |
| cavityFlexibleBottom | **Exclude** (revisit after the discrepancy is resolved; then at most supplementary) | Same-lineage digitised reference; unexplained monolithic/partitioned difference. |
| immersedHronTurekFsi2, membraneRoof, flexibleDamBreak, fillingElasticContainer, others | **Exclude** | See table above. |

A closing **4.4 Summary table** would list, per case: reference category,
observed orders, finest error, coupling agreement and cost. This is the
"suite card" readers will use.

---

## Missing evidence worth computing

### Essential before publication

1. **Resolve the cavityFlexibleBottom monolithic vs partitioned u_y discrepancy**
   (24% at equal force). Start by diffing the two case set-ups (solid model,
   E, ν, monitor point z, boundary conditions, `dampingCoeff`); no new runs
   needed initially. This affects the parent paper's Case 1 as well.
2. **Diagnose the womersleyTube temporal order (1.3–1.8).** Time study on mesh
   4 and/or at n = 400; isolate the interface-velocity condition and the
   start-up transient. If unresolved, report it explicitly.
3. **Single-version rerun of every main-paper study** (v2512 by default, fixed
   commit, Linux), with tabulated per-level values written to CSV. Archive the
   case set-ups and drivers (Zenodo).
4. **HronTurek FSI3 4x level** (Δt = 0.00025, 16 ranks). Gives a third level for
   reference-free orders and tests whether the 13% lift-amplitude error closes.
5. **3dTube level 3** (≈1.4M cells, >10 h on 6 ranks). A third level is needed
   to claim convergence (u_z,min still changes 1.85%).
6. **A tolerance (iteration-error) study on at least two main-paper cases**
   (e.g. ringAddedMass strong level and collapsibleChannel). Vary
   `outerCorrTolerance` over two decades to show the iteration error sits
   below the discretisation error. Only FSI2 has this today.
7. **Clarify reference provenance** (documentation, no runs). Which Featflow
   tables vs. 2006 summary values are quoted; which Richter (2015) level gives
   u_x = 5.95e-5 m, F_x = 1.33 N; mark Tiba et al. (2026) as a preprint.

### Desirable

1. Independent commercial-code references (COMSOL requested in PRs #497, #498,
   #501, #507). Priority: mokLidDrivenCavity (group B), 3dTube, blobInTreacle.
2. collapsibleChannel time-step study on fluid level 3 or 4, and a rerun in
   parallel once the high-order Jacobian is parallelised (including the
   4-rank run that gave a wrong displacement).
3. HronTurek FSI1 4x level (>14 h on 3 ranks) and FSI2 4x after improving mesh
   motion.
4. Quasi-monolithic solver (open PR #233) on ringAddedMass (strong),
   collapsibleChannel (including massless ρ_s = 0) and womersleyTube: a
   same-code but different-algorithm cross-check that also shows robustness.
5. ringAddedMass: explain the mesh-dependent IQN-ILS amplitude loss that Robin
   does not show.
6. mokLidDrivenCavity: diagnose the 128² divergence; Richardson estimate of the
   trough.
7. beamInCrossFlow: graded fluid-mesh family reaching the asymptotic range;
   tabulate per-level values; check Gillebaart et al. (2016) for an
   independent modified-form value.

### Optional

1. womersleyTube tethered variant (after #511) and high-order wedge solid
   (after #512).
2. ringAddedMass n = 3 mode and thinner ring (h/a sweep).
3. blobInTreacle time study with tighter solver/coupling tolerances, below
   Δt = 0.03125.
4. cavityFlexibleBottom level 4 with smaller Δt (if the case is retained).
5. Cross-fork reproducibility (OpenFOAM.org, foam-extend) for the cases that
   support it.
6. Robin-Neumann cost on HronTurek (#493) as a coupling-performance appendix.
