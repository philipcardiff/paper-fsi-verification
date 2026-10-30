# PHASE C0: INTERFACE DIAGNOSTIC (free-slip, impermeable ring-fluid interface)

Written 2026-10-29. **Sections 1 to 3 were written and committed to disk before the C0 simulation was submitted**; results are in section 4 onward and were added afterwards.

## 1. Question and design

D5 (rigid outer wall free-slip, interface no-slip `w = 0`) gave a wet/dry ratio error of −1.11 %. C0 changes **only** the ring-fluid interface tangential condition: the `w = 0` constraint is removed from the interface nodes so FEBio's natural conditions apply, while impermeability must be retained. Everything else is identical to D5: same mesh (ring 4×48, fluid 12×48 hex8), K = 1e8, ρ_f = 140, μ = 1.4e-5 Pa·s, E_mesh = 1 Pa, `rhoi = 1`, Δt = T_wet/100, same tolerances, same corrected initial velocity (the initial `w` at interface nodes stays 0 in both, so the absolute initial fluid velocity is identical), same single core, same pristine v4.13 binary, same analysis.

## 2. Is w·n = 0 a correct natural condition on this interface? (source and manual inspection, before running)

| Check | Evidence | Result |
|---|---|---|
| Is the normal relative flux the natural boundary quantity of the mixture domain? | Theory Manual §3.6, eq. 3.6.3: δW_ext contains `+∮ δJᶠ w_n da` and `δw · t^τ da`; the internal virtual work has the integrated-by-parts term `w · grad δJᶠ` (eq. 3.6.2). Source: `FEFluidFSIDomain3D::ElementInternalForce`, `fJ = (H·(… + grad e·w/J_f) + ∇N·w)·detJ`. | Yes: w_n appears only as a boundary integral multiplying the virtual dilatation. |
| What does omitting it mean? | With no prescribed w_n (`fluid normal velocity` load) and no dilatation BC at those nodes, the boundary term is zero: **w_n = 0 weakly**, and t^τ = 0 (frictionless tangentially). User Manual §4.20: "On a boundary where no conditions are prescribed explicitly, w_n = 0 and t^τ = 0, which represents a frictionless wall." | Applicable. |
| Is w measured relative to the moving ring at this interface? | The mixture domain's boundary is defined on the solid; the shared interface nodes carry the ring's mesh displacement and velocity (conformal mesh, 98 shared nodes, audited). Absolute fluid velocity v_f = v_s + w, so w_n = 0 means v_f·n = v_s·n, the kinematic (no-penetration) condition relative to the moving ring. | Yes. |
| Any special handling of interface nodes in the fluid solver? | `FESurface::m_bitfc` ("interface") is used only in `FEElasticSolidDomain`/`FEBiphasicSoluteSolidDomain` LM substitution for **solid-shell** interfaces, not for fluid DOFs. The `fluid-FSI traction` load only adds p·n and viscous traction to the displacement rows. | No hidden constraint; nothing re-imposes w at these nodes. |
| Essential/natural conflict? | Dilatation is not prescribed at any node (checked: no `zero/prescribed fluid dilatation` BCs). `w_z` is fixed on all fluid nodes, symmetry-plane normal components stay fixed, so no node has both an essential and a natural condition on the same component. | None. |
| Caveat that cannot be removed by reading | The condition is **weak** (satisfied in the Galerkin sense per nodal dilatation test function), so per-node leakage w_n ≠ 0 is possible. It must be measured, and it is (section 3). | Measured, not assumed. |

Input verification (`runs/phaseC/C0_wet_f1_wallslip_ifslip/`): regenerating D5 with the modified generator reproduces the stored D5 `.feb` **byte for byte**; C0 differs from D5 in exactly two node-set lines (`FixWX`, `FixWY`; 4 diff lines). In C0: `w_z` fixed on all 1274 fluid nodes (unchanged), `u_z` on all 1666 nodes (unchanged), `FixWX`/`FixWY` contain only the 26 + 26 symmetry-plane nodes (on x = 0 and y = 0 respectively), interface nodes keep only their symmetry-plane components (2 each), the `fluid-FSI traction` load is present, no dilatation BCs, node coordinates and all BC definitions identical.

## 3. Pre-registered consistency criteria (fixed before the run)

Frequency agreement is deliberately **not** a criterion for selecting the interface formulation. The selection rule: adopt the slip interface as the analytical-limit formulation only if it passes P1–P6; otherwise keep the no-slip interface. If both pass, the choice is made on physical justification (the reference's inviscid limit, in which tangential ring motion does not couple to the fluid), not on the smaller frequency error, and both are reported.

| # | Criterion | Threshold |
|---|---|---|
| P1 | **Leakage** relative to the ring's normal motion: rms over interface nodes and time of w_n divided by rms of the ring normal velocity at the same nodes; max over nodes and time of ∣w_n∣/V0; net leaked flux ∮ w_n dA (n = 0 content), peak, relative to peak ∮∣v_s,n∣ dA | ≤ 2 %; ≤ 10 %; ≤ 1 % |
| P2 | n = 2 pressure at the added-mass relation: layer-1 amplitude vs m_a·ω²·D, layer-0 amplitude; phase vs D | within 5 %; within 8 %; 178° ± 3° |
| P3 | Spurious pressure: n = 0 amplitude at ω; n = 2 residual after harmonic fit | ≤ 10 % of the analytic n = 2 amplitude (0.03 Pa); ≤ 10 % of the signal |
| P4 | Energy: total kinetic energy at same-sign displacement zero-crossings (integer periods) | within 1 % of initial |
| P5 | Solver quality: every step converged, no Jacobian/NaN/inversion message, volume: max ∣mean(J−1)∣ | 400/400; none; ≤ 1e-7 |
| P6 | Integrated n = 2 radial loading from the fluid vs analytic | within 5 % |

Metrics are computed from the existing logs by `tools/run_report.py`, `tools/pressure_modes.py` and the new `tools/interface_leakage.py`. For the D5 side, w_n = 0 is imposed strongly at interface nodes, so its leakage is identically zero by construction; the meaningful D5-vs-C0 comparisons are the frequency, ratio, pressure harmonics, force and energy.

STOP rule (from the brief): if the impermeable, frictionless interface turns out not to be implemented correctly (leakage fails P1 or the flux is not weakly zero), the model is unphysical and no further Phase C runs use it.

---

## 4. Results (written after the run; run `runs/phaseC/C0_wet_f1_wallslip_ifslip`, Slurm 13619, main, 1 core, 255.6 s, 400/400 steps, normal termination)

Pristine FEBio 4.13.0.32ae206ff (sha256 prefix `f5ae1ed6f347`), input hash in `results/diagnostic_matrix.csv`. D5 is the comparator (same input apart from the two node sets).

| quantity | D5 (interface no-slip) | **C0 (interface natural slip, impermeable weakly)** | reference / note |
|---|---|---|---|
| wet ω (rad/s), damped-harmonic fit 0–4T | 1.830673 | **1.854296** (zero-crossing 1.854295; window spread 2e-5 %) | exact wet 1.804513 |
| error vs exact wet ω | +1.45 % | **+2.76 %** | includes the dry ring's +2.58 % (hex8, B1) |
| wet/dry ratio (matched dry B1) | 0.698447 | **0.707459** | exact 0.706256 |
| ratio error | −1.11 % | **+0.170 %** | screening 0.3 %, not proof |
| effective added-mass ratio μ_eff | 1.050 | **0.998** | exact 1.008 |
| n = 2 pressure amplitude (layer 0 / layer 1, Pa) | 0.3356 / 0.3225 | **0.3254 / 0.3114** | analytic profile (below) |
| n = 2 pressure phase vs ring displacement | 178.3° | **178.1°** | 180° ideal; measured with a 5-parameter harmonic fit |
| measured ÷ exact radial pressure profile (layer 0 / layer 1) | 1.044 / 1.043 | **0.999 / 0.999** | profile factors 0.979 / 0.941 from the exact potential |
| integrated n = 2 radial force vs analytic (from layer-0 pressure) | +2.3 % | **−2.1 %** | the exact layer-0 profile factor is −2.1 % |
| n = 0 pressure at ω | 0.0000 Pa | **0.0000 Pa** | 2ω term 0.025 Pa is the physical second-order effect |
| n = 2 pressure residual after fit | 0.0108 Pa | **0.0097 Pa** (3 % of signal) | |
| fluid volume: max ∣mean(J−1)∣ | 2e-8 | **5.1e-10** | |
| total KE at integer periods / initial | 0.9986 to 0.9995 | **0.9995 to 0.9999** | half-period values 1.003 (corrected-IC cross term) |
| amplitude loss per period | 2.7e-5 | **3.1e-6** | |
| ring velocity vs d(displacement)/dt | 0.37 % | **0.21 %** | |
| solver | 400/400, 0 Jacobian/NaN | **400/400, 0 Jacobian/NaN, 0 "max iterations" events** | |

### Leakage through the interface (`tools/interface_leakage.py`; D5 gives exactly 0 by construction)

| metric | C0 | pre-registered limit |
|---|---|---|
| rms w_n ÷ rms ring normal velocity at the 98 interface nodes | **1.04 %** | ≤ 2 % |
| max ∣w_n∣ ÷ V0 | **1.25 %** | ≤ 10 % |
| net (n = 0) leaked flux, peak ÷ peak ring ∣flux∣ | **0.43 %** | ≤ 1 % |
| n = 2 leaked flux ÷ n = 2 ring flux (peak) | 1.2 % | (reported) |
| accumulated leaked volume, peak | 1.6e-6 m³ (≈ 2 % of the 8e-5 m³ swept per half cycle; 2e-5 of the fluid volume) | (reported) |

The impermeability is therefore satisfied weakly with ~1 % residual, as expected for a Galerkin natural condition. Whether it improves with refinement is not yet known and is checked at M2.

## 5. Pre-registered criteria, outcome (stated literally, with one flagged correction)

| # | criterion | outcome |
|---|---|---|
| P1 | leakage ≤ 2 % rms, ≤ 10 % max, net ≤ 1 % | **PASS** (1.04 %, 1.25 %, 0.43 %) |
| P2 (literal) | layer-1 n = 2 amplitude within 5 % of m_a ω² D; layer-0 within 8 %; phase 178 ± 3° | **FAIL at layer 1 as written** (0.941, −5.9 %); layer 0 0.978 passes; phase 178.1° passes |
| **P2′ (profile-corrected, flagged)** | the same amplitudes compared with the *exact radial profile* of the n = 2 pressure at the element-centre radii (0.979, 0.941) | **PASS** (0.999, 0.999; D5 fails it by +4.4 %) |
| P3 | n = 0 at ω ≤ 0.03 Pa; n = 2 residual ≤ 10 % | **PASS** (0.0000 Pa; 3 %) |
| P4 | KE at integer periods within 1 % | **PASS** (0.05 %) |
| P5 | converged, no messages, volume ≤ 1e-7 | **PASS** (5e-10) |
| P6 | integrated n = 2 force within 5 % of analytic | **PASS** (−2.1 %, equal to the exact profile factor) |

**Disclosure about P2.** The literal P2 compared layer pressures with the *interface* value and so ignored the physical radial decay of the n = 2 pressure. That decay (2 % at layer 0 and 6.5 % at layer 1) was already computed in `PHASE_B_COUPLING_AUDIT.md` item G before this run, so the correct comparison was known; my pre-registered threshold was mis-specified. I record the literal result (FAIL at layer 1) and the corrected result (PASS) and do not treat the correction as a relaxation: P2′ removes an allowance rather than adding one, and the corrected comparison separates the two runs sharply (C0 0.1 % from the exact profile; D5 4.4 % high). No other criterion was changed.

## 6. Selection

**Selected: the slip (natural, frictionless, impermeable-in-the-weak-sense) interface, together with the free-slip outer wall**, for the Phase C baseline, labelled **"inviscid-limit formulation"**. Justification, independent of frequency:

1. It is the consistent analytical limit: the reference states the unresolved viscous layer makes the computed flow effectively inviscid and that tangential ring motion does not couple to an inviscid fluid; strong nodal no-slip over a cell ~100× the Stokes layer (δ ≈ 3.3e-4 m vs 3.75e-2 m) is a discretisation artefact whose signature (4 % pressure overshoot, spurious n = 0 pressure, mesh-dependent added mass) was observed in R1/D2/D5.
2. FEBio's natural condition is exactly the intended kinematic condition (section 2) and leakage measured at the registered level.
3. The n = 2 pressure follows the exact radial potential profile to 0.1 % at two layers, and the integrated loading follows the exact profile to the same level.
4. Energy and volume are conserved better than in D5.

The smaller frequency error (+0.17 % vs −1.11 %) was **not** the basis of the choice; D5's ratio error would have been "more conservative". Caution: both runs share the same coarse mesh and a coarse dry ring (+2.58 %), so a ratio near zero may include cancelling errors; the M2 results are the test.

**Labelling requirement.** Every Phase C table must say: *inviscid-limit formulation (free-slip rigid wall, frictionless impermeable ring interface, K = 1e8 Pa), not a replica of the finite-viscosity no-slip benchmark.* The no-slip results (B2, R1, D2, D4, D5) are preserved and unchanged.

## 7. Stop conditions

No unphysical leakage, invalid traction, unstable pressure or unexplained behaviour; nothing blocks Stage C1.
