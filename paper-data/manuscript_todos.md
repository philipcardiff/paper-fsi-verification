# Manuscript TODOs: FSI verification-benchmark paper

Created during the first drafting pass of Sections 4–7 (2026-10-03), and
re-prioritised after the compression and methodological-correction pass on the
same day. The evidence base is [`benchmark_inventory.md`](benchmark_inventory.md)
and [`benchmarks.yaml`](benchmarks.yaml) (solids4foam `development` at
`7c928c55`; no FSI tutorial changes up to `c706be68`). Each item corresponds to
a statement or a `\todo{}` in the LaTeX. Section references use LaTeX labels
because the section numbers may change.

Ground rules:

- **Monolithic results are not used.** The quasi-monolithic solver (open PR
  #233) and all of its results, including `simulationData/cavityFlexibleBottom/`
  and the parent paper's Case 1 figures, are not evidence for anything.
- **Proportionate verification.** For each benchmark and QoI, show that the
  error sources that materially affect that QoI are controlled. Do not fill a
  full mesh × Δt × coupling × tolerance matrix for every case.
- **No broad reruns for version consistency.** Rerun on a different version
  only where the version difference is scientifically material
  (womersleyTube, item 1).

---

## 1. Highest priority (needed to support or settle current main-text claims)

| # | Case | Computation | Why it matters |
|---|---|---|---|
| 1 | womersleyTube | **Finer temporal investigation.** Run 50/100/200/400 steps per period on mesh factor 4. On mesh 2, switch the interface-velocity condition to Euler, and compare an exact-initial-state start with a start from rest plus a ramp. Log the start-up transient amplitude against Δt. Also reconcile the 0.26%-of-amplitude Robin–Neumann shift between v2412 and v2512 (PR #513 comment). That shift exceeds the finest-mesh errors in `tab:womersley`. | The sub-nominal time order (1.29–1.84) is currently reported as an unresolved observation (`sec:womersley`, `sec:cross_convergence`). Report a diagnosis if one is found; otherwise keep it open. **Do not tune the result away.** |
| 2 | Hron–Turek FSI3 | Run the **4x level** (Δt = 0.00025 s, about 16 ranks), IQN-ILS, with the same QoI window. Report three-level successive-difference orders for all seven QoIs. | Without it, FSI3 stays at "two levels, not demonstrably asymptotic" (apparent orders 1.1–3.8). It also tests whether the 13.3% lift-amplitude error closes. Required for any strong spatial-convergence claim. |
| 3 | 3dTube | Run **mesh level 3** (about 1.4M cells; >10 h on 6 ranks), BDF2, with Δt halved. | Two levels cannot support a convergence claim, and u_z,min still changes by 1.85% (`sec:3dtube`). |
| 4 | ringAddedMass and/or collapsibleChannel | **One or two representative coupling-tolerance sweeps.** Vary the coupling tolerance over at least two decades (e.g. 1e-4 → 1e-7) at fixed mesh and Δt. Use the ring strong level (mesh 2, 100 steps/period) and collapsible fluid level 3, with IQN-ILS and Robin–Neumann. | The claim "iterative error ≪ discretisation error" currently rests on FSI2 alone (`sec:iterative_errors`). The collapsible Robin/IQN-ILS gap does not shrink with refinement. |
| 5 | beamInCrossFlow | **Redesigned graded mesh study.** Use a fluid-mesh family with near-body/boundary-layer grading, at least 3 levels at a constant ratio, at fixed Δt, for both forms. Tabulate u_x, u_y, u_z(A), F_x and F_y per level in CSV. Add one separate Δt check. | The existing uniform family is not asymptotic: the last change is about 4.6–4.7% and F_y is non-monotone. This tests the displacement-versus-force hypothesis. It is a condition for keeping this provisional case (`sec:beam`). |
| 6 | all main-text cases | Record the exact solids4foam commit, OpenFOAM version and platform for every number in Sections 5–6. This is documentation only, with no runs. | Reproducibility (`sec:qoi`, `sec:infrastructure`). |

## 2. Independent COMSOL references (in priority order)

An independent converged solution, with its own mesh and Δt study, adds what
solids4foam alone cannot. It is **not** needed for every benchmark.

1. **3dTube**: BDF2-type time integration at a converged Δt. It would replace the same-lineage Tuković (2018) reference and test the time-discretisation explanation of the FE/FV literature difference from the FE side.
2. **beamInCrossFlow**: at least the modified (large-deformation) form, whose only reference is same-lineage. The original form is of moderate value.
3. **mokLidDrivenCavity**: the group-B definition, with peak, trough and mean over t = 20–70 s. It would make this a category-B case and decide core-suite membership (`sec:core_suite`).
4. **Hessenthaler**: may become high priority once the solids4foam implementation is mature, since it would give the capstone case a category-B reference.

Lower priority: blobInTreacle (already being arranged) and cavityFlexibleBottom
(part of its rehabilitation route; see §3).

## 3. Hessenthaler verification campaign (when the implementation is ready)

Set up the case and the verification driver, and define QoIs. Run the mesh
(≥3 levels), Δt (≥3 levels) and coupling-tolerance studies. Keep the
experimental MRI data separate as validation context only. `sec:hessenthaler`
must contain no numbers until this exists.

## 4. Strengthening (valuable, not essential)

| # | Case | Computation | Note |
|---|---|---|---|
| S1 | Hron–Turek FSI1 | 4x level (>14 h on 3 ranks). | Gives a reference-free order for u_y. Until then, write "converging towards Featflow". |
| S2 | collapsibleChannel | Δt study on fluid level 3 or 4. | The finer-mesh temporal error is unmeasured. |
| S3 | Hron–Turek FSI3 | Independent Δt study on 2x. | There is no separate Δt study for FSI2/3. |
| S4 | mokLidDrivenCavity | Diagnose the 128² divergence. Compute Richardson estimates for peak and trough, and quantify the IQN-ILS vs Aitken difference. | The trough order (<0.5) is poorly established. |
| S5 | cavityFlexibleBottom | Partitioned level 4 at reduced pseudo-Δt, a pseudo-Δt check and a coupling comparison. Explain the ~5% u_y offset at equal F_y by diffing the set-up against Tuković (2018). | This is the rehabilitation route, together with a COMSOL reference. Partitioned evidence only. |
| S6 | blobInTreacle | Time study below Δt = 0.03125 s with tolerances tightened by about two decades. | The present orders (2.4–2.5) are super-nominal and limited by the tolerance floor. |
| S7 | ringAddedMass | Investigate the mesh-dependent IQN-ILS amplitude loss. | Links to item 4. |

## 5. Explicitly **not** essential for this paper

- Resolving any discrepancy involving the experimental monolithic/quasi-monolithic solver.
- Obtaining COMSOL results for every benchmark.
- Rerunning every benchmark on one OpenFOAM version for cosmetic consistency.
  A single-version archival rerun can be a later reproducibility task, once the
  final benchmark suite is frozen.
- Optional extensions: FSI2 4x after a mesh-motion fix; a tethered womersleyTube
  variant (#511) and a high-order solid on wedges (#512); a ring n = 3 / h/a
  sweep; a parallel collapsibleChannel rerun after the Jacobian fix, including a
  recheck of the 4-rank linear-reconstruction run; 3dTube Robin–Neumann
  iteration growth at small Δt; cross-fork reproducibility.

## 6. Citation / literature checks

Missing from `bibliography.bib` (marked `\citetodo{}`):

- Roache (1994, 1998), Oberkampf & Roy (2010), Richardson (1911).
- Womersley (1955, 1957); Filonova et al. (2020, IJNMBE 36:e3266).
- Brennen (1982), NCEL CR 82.010.
- Jensen & Heil (2003, JFM); Heil & Hazel (2006), oomph-lib.
- Wall (1999), Mok (2001), Gerbeau & Vidrascu (2003), Valdés (2007), the Kratos example, and **Tiba et al. (2026), arXiv:2609.16876, which must stay marked as a preprint**.
- Formaggia et al. (2001); Fernández & Moubachir (2005).
- Liu et al., arXiv:1401.0082.
- Hessenthaler et al. (2017, IJNMBE) and the companion numerical paper(s).
- The citable source of the Featflow FSI tables, as distinct from the Turek & Hron (2006) summary values.

Checks:

- Verify every entry of `tab:hronturek_ref` against the 2006 chapter and the Featflow tables (not the solids4foam README). Identify which later papers quote which set.
- Richter (2015): which level gives u_x(A) = 5.95e-5 m and F_x = 1.33 N, and whether the values are extrapolated.
- Gillebaart et al. (2016): whether they tabulate the modified beamInCrossFlow configuration.
- 3dTube: confirm the Lozovskiy (2019) and Eken (2016) time-integration settings from the papers.
- Mok: confirm the group-A/B classification and the placement of Kassiotis (2011).
- Check that `Juretic2010` is the right origin citation for the cavity case.
- Optional: 1–2 FSI-specific verification references for §4 (e.g. MMS for FSI).

## 7. Figure / table work

| # | Item | Status |
|---|---|---|
| F1 | Geometry sketches for ringAddedMass, womersleyTube, collapsibleChannel and mokLidDrivenCavity. | missing |
| F2 | Publication version of `fig:collapsible_mesh` (clean legend, insets, error-vs-h panel). | interim |
| F3 | Publication version of `fig:mok_refs` (add the 96² history and the group-B envelope; separate the traction-free pilot). | interim |
| F4 | Log–log convergence plots (ringAddedMass; womersleyTube amplitude vs phase/attenuation; collapsibleChannel) for `sec:cross_convergence` / `sec:qoi_dependence`. | missing |
| F5 | Hron–Turek QoI-vs-h figure (after item 2), with an FSI2 tolerance panel. | missing |
| F6 | 3dTube u_r(A) panel: BDF2 vs implicit Euler (Δt = 1e-4) vs the literature curves. | missing |
| F7 | beamInCrossFlow figures: regenerate after item 5. | deferred |
| F8 | Generate `tab:hronturek` and `tab:3dtube` directly from driver CSV output. | to do |
| F9 | ringAddedMass parameter table and full-μ results. | to do |
| F10 | Add cost and finest-error columns to `tab:suite_overview` / `tab:suite_tiers` once items 2–3 are done. | to do |

## 8. Manuscript-structure items (not computations)

- The title, abstract, Introduction and Sections 2–3 are still from the parent monolithic paper. Section 3 must describe the partitioned procedures, solvers and time schemes actually used.
- Section 8 (Discussion) is a skeleton; Section 9 (Conclusions) is unchanged.
- `sectionTestCases.tex` is kept on disk but is not `\input`.
- Data Availability: update to the `development` branch, the commit and the Zenodo DOI.
- Decide whether to report GCIs or only Richardson estimates (`sec:convergence_method`).
