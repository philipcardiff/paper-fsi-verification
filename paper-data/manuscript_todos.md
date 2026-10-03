# Manuscript TODOs: FSI verification-benchmark paper

Created during the first drafting pass of Sections 4–7 (2026-10-03).
The evidence base is [`benchmark_inventory.md`](benchmark_inventory.md) and
[`benchmarks.yaml`](benchmarks.yaml) (solids4foam `development` at `7c928c55`).
It was re-checked against `development` at `c706be68`, where there are no FSI
tutorial changes. Every item below is exposed by a statement or an
`\todo{}` in the LaTeX.

Ground rules carried into the draft:

- **Monolithic results are not used.** The quasi-monolithic solver (open PR
  #233) and all of its results, including `simulationData/cavityFlexibleBottom/`
  and the parent paper's Case 1 figures, are not used as evidence.
- **No broad reruns for version consistency.** A blanket single-version rerun
  is not recommended. Reruns are listed only where a version difference is
  scientifically material (womersleyTube, item E2).

"COMSOL value" says whether an independent COMSOL (or other independent-code)
solution would add something that the solids4foam computation alone cannot.

---

## 1. Essential to support an existing main-text claim

| # | Case | Exact computation needed | Reason | COMSOL value |
|---|---|---|---|---|
| E1 | womersleyTube | (a) Time study at 50/100/200/400 steps per period on mesh factor 4. (b) The same on mesh 2 with the interface-velocity condition switched to Euler, and with a start from rest plus ramp versus the exact initial state. (c) Log the start-up transient amplitude against Δt. | §5.1.2 reports a sub-nominal time order (1.29–1.84) and a start-up transient that grows as Δt shrinks. The text must either keep this as an open finding (current draft) or report a diagnosis. Do **not** tune the result away. | None. The reference is exact; this is a property of the coupled scheme. |
| E2 | womersleyTube | Rerun the Robin–Neumann mesh-4/n200 and mesh-2/n100 cases on the OpenFOAM version used for the final tables, and reconcile them with the 0.26%-of-amplitude shift between v2412 and v2512 (PR #513 comment). | The shift is larger than the finest-mesh errors quoted in Table `tab:womersley` (0.07–0.12%), so the version does matter for this case. | None. |
| E3 | Hron–Turek FSI3 | Run the 4x level (Δt = 0.00025 s, about 16 ranks), IQN-ILS, with the same QoI window. Report three-level successive-difference orders for all seven QoIs. | §5.3.1 can currently say only "two levels, not demonstrably asymptotic", with apparent orders of 1.1–3.8. The 4x level also tests whether the 13.3% lift-amplitude error closes. | Low. Featflow is already independent. |
| E4 | 3dTube | Run mesh level 3 (about 1.4M cells; estimated >10 h on 6 ranks), BDF2, Δt halved. | Two levels cannot support a convergence claim, and u_z,min still changes by 1.85% (§5.3.2). | **High.** An independent converged BDF2-type solution would replace the same-lineage Tuković (2018) reference. It would also test the time-discretisation explanation of the 21–24% FE/FV literature split from the FE side. |
| E5 | ringAddedMass, collapsibleChannel | Coupling-tolerance sweep: vary `outerCorrTolerance` over at least two decades (e.g. 1e-4 → 1e-7) at fixed mesh and Δt. Use the ring strong level (mesh 2, 100 steps/period) and the collapsible fluid level 3, for both IQN-ILS and Robin–Neumann. | §4.2 and §6.3 rely on "iterative error ≪ discretisation error". Only FSI2 tests this directly today, and the collapsible Robin/IQN-ILS gap does not shrink with refinement. | None. |
| E6 | Hron–Turek FSI1 | Run the 4x level (estimated >14 h on 3 ranks). | §5.2.2 needs a three-level, reference-free order for u_y (4.64% error on 2x; apparent order 1.5–1.7). Without it, the draft must stay at "converging towards Featflow". | None (Featflow 6 levels). |
| E7 | collapsibleChannel | Time-step study (0.0125 → 0.0015625 s) on fluid level 3 or 4. | The temporal error on the finer meshes is not measured. §5.2.1 attributes the 0.5% floor to the reference and model form, so the time error there must be shown to be small. | Low. |
| E8 | all main-text cases | Record and state the exact solids4foam commit, OpenFOAM version and platform for every number in Sections 5–6 (documentation only; no runs). | Required for reproducibility (§4.3, §7.3). | – |

## 2. High-value strengthening

| # | Case | Exact computation needed | Reason | COMSOL value |
|---|---|---|---|---|
| H1 | mokLidDrivenCavity | Independent converged solution of the **group-B** definition, with its own mesh study (at least 3 levels) and Δt study. Report mean peak, trough and mean over t = 20–70 s. | It would turn a C-weak, self-convergence-only case into category B, and settle whether the +3.3% peak offset from the group-B envelope belongs to the literature or to solids4foam. It would also justify core-suite membership (§7.1). | **High.** No published source has a convergence study. |
| H2 | mokLidDrivenCavity | Diagnose the 128² divergence; run 128² (or a 96²→144² sequence with ratio 1.5); compute Richardson estimates for peak and trough; quantify the IQN-ILS vs Aitken difference (currently "same response", no number). | Peak order ≈1.2 and trough order <0.5 are poorly established. | – |
| H3 | beamInCrossFlow | Graded fluid-mesh family with boundary-layer resolution near the plate (at least 3 levels, constant ratio), at fixed Δt, for both forms. Tabulate u_x, u_y, u_z(A), F_x and F_y per level in CSV. Add a separate Δt check on one mesh. | §5.4.1 is provisional because uniform refinement is non-asymptotic (last change about 4.6–4.7%; F_y non-monotone). This would also test the hypothesis that displacement and force converge differently because of near-wall resolution. | **High** for the modified (large-deformation) form, whose only reference is same-lineage. Moderate for the original form. |
| H4 | blobInTreacle | Independent converged solution at t = 1 s and at steady state; compare top-point u_x, u_y and interface shape. | The steady offset (−7.4% u_x, +20% u_y magnitude against the preprint) is currently attributed to the coarse reference without evidence. | **High.** COMSOL is already being arranged. |
| H5 | blobInTreacle | Time study below Δt = 0.03125 s with coupling and solver tolerances tightened by about two decades. | The present orders (2.39–2.54) are super-nominal and limited by a ~1e-5 m tolerance floor. The draft claims only "at least second order over this range". | None. |
| H6 | cavityFlexibleBottom | (a) Explain the ~5% steady u_y offset at equal F_y relative to Tuković (2018). Diff the solid model, material, monitor point, BCs and `dampingCoeff` against the 2018 set-up, and the change of the coarse-mesh value from −0.219 m (earlier solids4foam) to −0.20246 m (now). (b) Partitioned level 4 (Δx = 0.0125) with reduced pseudo-Δt. (c) Partitioned pseudo-Δt check and coupling (Aitken vs IQN-ILS) comparison. | Promotion criterion in §5.3.4. Monolithic data must not be used for any of this. | **High.** An independent converged solution of the identical set-up is the only way to decide between the published curve and the present converged value. |
| H7 | Hron–Turek FSI3 | Independent time-step study at 2x (Δt = 0.001, 0.0005, 0.00025 s). | No separate Δt study exists for FSI2/3; Δt is scaled with h. | None. |
| H8 | ringAddedMass | Investigate the mesh-dependent IQN-ILS amplitude loss (up to 2.1%/period, halving with refinement, absent with Robin–Neumann). For example, vary IQN-ILS reuse and filter settings, and the tolerance (links to E5). | Reported as unresolved in §5.1.1 and §6.3. | None. |
| H9 | Hessenthaler | Set up the case and the verification driver; define QoIs; run mesh (≥3), Δt (≥3) and coupling-tolerance studies. | §5.4.2 is a placeholder; the cross-benchmark conclusions on cost and dimensionality (§6.5) depend on it. | Moderate. An independent numerical solution would give the capstone case a category-B reference. The experimental (MRI) data are validation only. |

## 3. Optional

| # | Case | Exact computation needed | Reason | COMSOL value |
|---|---|---|---|---|
| O1 | Hron–Turek FSI2 | 4x level after improving mesh motion (non-orthogonality reaches 58° by 10.5 s; IQN-ILS fails at 11.4 s on 2x). | Lift amplitude does not converge (7.9% → 7.8%). | None. |
| O2 | womersleyTube | Tethered-wall variant (after issue #511) and high-order solid on wedges (after #512). | Allows a comparison with Womersley's classical constrained tube. | None. |
| O3 | ringAddedMass | n = 3 mode; h/a sweep. | Shows the O((h/a)²) thin-ring trend. | None. |
| O4 | collapsibleChannel | Parallel rerun once the high-order Jacobian is parallelised, including a recheck of the 4-rank linear-reconstruction run that gave a wrong displacement. | Robustness; the code correctness issue is noted as a limitation in §5.2.1. | None. |
| O5 | 3dTube | Investigate why Robin–Neumann iterations rise to 7.4 at the smallest Δt. | Noted in §5.3.2. | None. |
| O6 | all | Cross-fork reproducibility (OpenFOAM.org, foam-extend) where supported. | Supporting material only. | None. |

## 4. Citation / literature checks

Missing from `bibliography.bib` (marked `\citetodo{}` in the text):

- Roache (1998) *Verification and Validation in Computational Science and Engineering*; Roache (1994) GCI; Oberkampf & Roy (2010); Richardson (1911).
- Womersley (1955, Phil. Mag. 46:199–221); Womersley (1957, Phys. Med. Biol. 2:178–187); Filonova et al. (2020, IJNMBE 36:e3266).
- Brennen (1982), NCEL report CR 82.010 (added mass).
- Jensen & Heil (2003, JFM) collapsible channel; Heil & Hazel (2006) oomph-lib.
- Wall (1999, thesis); Mok (2001, thesis); Gerbeau & Vidrascu (2003); Valdés (2007, thesis); Kratos example (Zorrilla); **Tiba et al. (2026), arXiv:2609.16876. Keep it marked as a preprint.**
- Formaggia et al. (2001); Fernández & Moubachir (2005) for the 3dTube origin.
- Liu et al. arXiv:1401.0082 (preprint version of Liu2014 with the figures used).
- Hessenthaler et al. (2017, IJNMBE) and the companion numerical paper(s).
- Featflow FSI benchmark tables: identify the citable source of the tables, as distinct from the Turek & Hron (2006) summary values.

Checks on existing statements:

- **Hron–Turek:** verify every entry of Table `tab:hronturek_ref` against the original 2006 chapter and the Featflow tables, not against the solids4foam README. Identify which later papers quote which set.
- **Richter (2015):** confirm which refinement level gives u_x(A) = 5.95e-5 m and F_x = 1.33 N, and whether these are extrapolated values.
- **Gillebaart et al. (2016):** check whether they tabulate values for the modified beamInCrossFlow configuration.
- **3dTube:** confirm the Lozovskiy (2019) and Eken (2016) time-integration settings (implicit Euler, Δt = 1e-4 s) from the papers themselves.
- **Mok:** confirm the group-A/group-B classification and the Kassiotis (2011) placement against the original papers.
- Check that the `Juretic2010` citation (cavity origin) is appropriate.
- §4 intro: add 1–2 FSI-specific verification references, e.g. method-of-manufactured-solutions work for FSI, if the authors want it.

## 5. Figure / table work

| # | Item | Status |
|---|---|---|
| F1 | Geometry sketches for ringAddedMass, womersleyTube, collapsibleChannel, mokLidDrivenCavity (only HronTurek3, 3dTube, beamInCrossFlow, cavityFlexibleBottom and blobInTreacle exist in `figures/`). | missing |
| F2 | Replace `figures/verification/collapsibleChannel_mesh_study_wallMid.png` (copied from solids4foam, raw run labels) with a publication figure: clean legend, insets at first trough and second peak, error-versus-h panel. | interim |
| F3 | Replace `figures/verification/mokLidDrivenCavity_references_and_pilots.png` with a figure that adds the converged 96² history and the group-B envelope, and separates the traction-free pilot. | interim |
| F4 | Convergence plots (error or successive difference vs h and vs Δt, log–log, with slope-2 guides) for ringAddedMass, womersleyTube (amplitude vs phase/attenuation on one plot) and collapsibleChannel. These would be the key figures for §6.1–6.2. | missing |
| F5 | Hron–Turek: once 4x exists, a QoI-vs-h figure (drag mean, u_y amplitude, lift amplitude) showing the slow lift convergence. Add a tolerance-sweep panel (FSI2). | missing |
| F6 | 3dTube: u_r(A) history panel with BDF2 (converged), Euler Δt = 1e-4 (present code) and the Lozovskiy, Eken and Tuković curves. This is the visual evidence for the literature-discrepancy result. | missing |
| F7 | beamInCrossFlow: the stored PNGs (`original_mesh_t8_backward_vs_references.png`, `modified_...`) are diagnostic only. Regenerate after the graded-mesh study (H3). | deferred |
| F8 | Regenerate Tables `tab:hronturek` (1x errors and apparent orders were derived in the inventory) and `tab:3dtube` directly from driver CSV output, so that every table entry is reproducible from a file. | to do |
| F9 | A ringAddedMass table of the full μ range (weak and moderate rows) and a parameter table. The draft shows dry and strong rows only. | to do |
| F10 | Suite-overview Table `tab:suite_overview`: add cost and finest-error columns once E3/E4/E6 are done ("suite card"). | to do |

## 6. Manuscript-structure items (not computations)

- The title, abstract, Introduction and Sections 2–3 remain those of the parent monolithic paper. Section 3 must be revised to describe the partitioned procedures (IQN-ILS, Robin–Neumann, Aitken), the fluid and solid solvers, and the time schemes actually used.
- Section 8 (Discussion) is a skeleton; Section 9 (Conclusions) is unchanged (WIP).
- `sectionTestCases.tex` (parent-paper monolithic results) is preserved on disk but no longer `\input`.
- The Data Availability statement still points to the `feature-petsc-snes` branch and `solid-benchmarks`. Update it to the `development` branch, the commit and the Zenodo DOI.
