# Cleanup of Sections 5–8, Conclusions and appendices (Philip's principles)

Branch `text/sections-5-7-cleanup`. This applies the principles of `philip_notes.md` and `philip_notes_sections1-4_response.md` to Section 5 (benchmark suite), Section 6 (cross-benchmark analysis), Section 7 (recommended suite), Section 8 (Discussion), the Conclusions, Appendix A (beam) and Appendix D (lessons). In `article.tex` the order is 5 Benchmark suite, 6 Cross-benchmark, 7 Recommended suite, 8 Discussion, 9 Conclusions.

The list below gives, per subsection, what was **removed**, **moved** (and where) or **rewritten**. Nothing listed as moved was dropped: the text now sits in Appendix D (`appendixLessons.tex`) unless stated otherwise.

Build: `latexmk -pdf article.tex` gives no errors and no undefined references. **93 pages before, 89 after.** The warnings are the duplicate hyperref destinations and the three overfull boxes that were already there.

## Section 5 (`sectionBenchmarkSuite.tex`)

### Table 2 (suite overview)
- Rewritten: the FSI3 row now gives the final sequence only (4 levels; 2 time steps on 2x′), with a footnote `\hl{[FSI3: ring+box, to be replaced]}`. Removed: "3; 4 (corrected)".

### 5.1.1 Ring added-mass problem
- Rewritten: the FEBio "Resolution of the viscous layers" paragraph is now two sentences. It keeps the finding that FEBio's nodal no-slip needs the Stokes layer resolved (−6.2% on 12 uniform radial cells; +0.001% with a first cell of δ/4).
- Removed: "Source inspection found no defect in the enforcement of the no-slip condition or in the traction transfer…" (this was the contributor's own debugging), and "reduced the net energy creation of the unresolved cases". The supplement keeps the energy detail.

### 5.1.2 Womersley tube
- Rewritten:
  - The *Time* item now states the final flux-consistent end treatment, its formal O(Δt h) property and the fixed-mesh limit (3200/6400 steps).
  - It states in one sentence that the usual mixed end condition gives orders 1.29–1.84 and 0.75–1.79.
  - The *Lesson* no longer tells the story ("first showed … traced …").
- Moved to App. D, new paragraph "Flux consistency at truncated boundaries":
  - the end-flux formula;
  - the diagnostics that ruled out the coupling (IQN-ILS/RN, fluid-only, mesh, start-up, window, ddtCorr, tolerances);
  - the before/after time errors (1.6, 0.9, 3.8e-3 → 1.3e-4, 3.3e-5, 1.6e-4);
  - the exact-flux variant (2.00–2.08, not robust);
  - the ddtCorr floor (~1e-6);
  - the zero-boundary-value start-up transient (3.4%) and the 0.26% v2412/v2512 difference it explained. The release names are now written as "two OpenFOAM releases".
- Removed from the main text: "All values are from the verification driver of the case with this configuration", and the build in the Table 6 caption (`fluxConsistentPatches`, v2512, Linux).

### 5.1.3 Manufactured solution
- Rewritten:
  - "Partitioned fluid" now states the four findings directly: the errors belong to the fluid, the outlet causes the corner error, ddtCorr makes the error grow as Δt falls, and the interior is second order with an outlet layer. All numbers are kept.
  - "corrections under review" is now "three corrections, not yet in a released version" (paragraph, captions, Table 19 row).
- Moved to App. D ("Errors that only an exact finite-deformation solution exposes"):
  - the list of ruled-out candidates (boundary non-orthogonality, gradient scheme, interface p and U conditions, skew-corrected flux, exact bottom-wall gradient);
  - the open-cause candidates (single outer corrector, outlet/inlet mass balance, mesh flux);
  - the near-uniform offset at Δt = 6.25e-4 (mean 3.4e-2);
  - the quasi-monolithic attempts to add k·∇p to the boundary pressure only or to the traction only.
- Consolidated: the "Non-orthogonal correction of the moving-wall velocity" paragraph now points to Sec. 3.1, which holds the measured FSI3 effect. Its old statement that the effect "has not been quantified" contradicted Sec. 3.1. The MMS-specific number (QM velocity condition changes the errors by <1%) is kept.

### 5.2.1 Collapsible channel
- Rewritten: two sentences that both said "cannot be attributed to either code" are now one sentence.
- The "remains to be re-checked" 4-rank run is now `\hl{[pending: re-check of the parallel run]}`.

### 5.3.1 Turek–Hron family (main target)
- Rewritten:
  - The FSI3 todo no longer says "level 3 running / expected 24–25 Nov".
  - The Featflow-uncertainty part of *Reference* is now four sentences instead of six. The numbers are unchanged.
- **Removed from the main text:**
  - the baseline-set-up Table (`tab:hronturek`, 1x/2x/4x with velocity-Laplacian motion and zero-gradient wall pressure);
  - the items "FSI3 displacements" (baseline), "FSI3 frequency and drag mean" (baseline) and "FSI3 time step and coupling" (baseline);
  - the whole "Diagnosis of the FSI3 displacement order" (solid alone, CFD3, code defects, prescribed motion, surface decomposition, corrected fluid under prescribed motion on both families).
- **Moved to App. D ("Path-dependent mesh motion")**:
  - the baseline displacement sequence (29.85/34.18/36.23 mm, order 0.9–1.1);
  - the 2x crossing of Featflow (+3.5%, −7.1% at 4x; two-level orders 2.6–3.8);
  - the solid-alone orders (1.4–1.6 → 1.8–1.9) and the solid refinement inside FSI3 (+0.33%; lift/drag −12.3/−5.5%);
  - the CFD3 agreement (0.16%, 0.05%);
  - "defects have no effect in 2-D";
  - Q_in under prescribed motion (50.3/72.2/81.0, order 1.3; treatments 1.1–1.3; corners 22%);
  - the drift (+0.3 N/m/s) and the displacement-Laplacian fix (±0.03).
- Not carried over:
  - the CFD3 detail beyond drag mean and lift amplitude (drag amplitude +2.0%, frequency +1.2%; rigid-flag drag/lift sequences 446.01… and 296.0…);
  - the Q_in values of the corrected treatments on the baseline family (43.66, 76.17, 81.22, 81.96) and on the ring+box family (67.18 … 77.16);
  - the 4–5° rotation per level;
  - the 6% family difference;
  - "net fluid work per cycle".

  The main text keeps the one consequence that matters, the −2.7% change of the in-phase load from 4x′ to 8x′, which was already quoted in Sec. 6.2. Everything else is in `paper-data/hronturek_fsi3_4x_findings.md`.
- **Wall-condition history:**
  - Removed: "first computed with an earlier version of the moving-wall velocity condition … subsequently corrected", and every "earlier wall condition" caveat in the text.
  - Each instance is replaced by one statement: the 8x′ level and the time/coupling column were computed with the "uncorrected wall-normal velocity derivative". It is stated once in the table caption and footnote and once in the FSI3 displacement item, together with its effect (1.8%, 1.5%, 0.86% on 1x′–4x′, ≈0.4% expected at 8x′).
  - The 2x′ A/B numbers (−3.1, −1.5, −3.8, −2.2%; other cases <1e-3, ~1e-4, ≤0.4%) moved to App. D ("Wall pressure and wall velocity conditions"), which before gave only "1.5–3.8%".
- Removed: the earlier-wall-condition 2x′→4x′ changes (−0.14%, +0.4%), the "apparent plateau" story, and "did not survive the correction of that condition".
- Kept, as a short item on the "uniform family": the coupling studies of the baseline family (lift noise 7.3 vs 1.7 N/m rms; amplitudes +15.0/+18.1%; 1e-6 noise eightfold lower; IQN-ILS stall 5–9e-6; RN vs IQN-ILS 5.1/1.9%, 3.5/2.8%). Sec. 6.3 relies on them, and they do not depend on the mesh-motion set-up.
- FSI2 items: four items are now one, with the same numbers.
- FSI1 with the high-order solid: tightened; "baseline set-up" wording removed.
- Moved to App. D (Robustness): "Limit-cycle restarts of the solid diverge…".
- Todos:
  - The "8x′ / high-order pair / diagnosis" decision is now limited to the 8x′ level, the cross run and the second-solid-code comparison. The diagnosis is resolved: it moved to App. D. The 2x′/4x′ high-order-solid pair appeared only in the todo and is not mentioned anywhere now.
  - The FSI2/evidence todo no longer gives the "Section 12 retracted…" history.
- Lesson: rewritten. It states the reference findings, the 4x′ differences from Failer and Richter, and that our unresolved sequence carries the error. The "disagrees vs not resolved" argument is consolidated in Sec. 6.1, and the Lesson points there.
- Removed: the build in the Table `tab:fsi1` caption (v2512, MeluXina).

### 5.3.2 3-D pressure-pulse tube
- Rewritten:
  - The two defects (a)/(b) are now one *Platforms* item. It keeps them briefly as evidence (the 1.5–3% u_z,min platform difference; weights switching on between levels 1 and 2) and points to Sec. 6.4 and App. D for the mechanisms.
  - The meshes used for the time study and the replication are now described as "a variant of the meshes whose wall arcs are defined to fewer digits", not "earlier meshes".
  - The replication result is kept in one sentence.
- Removed from the main text:
  - the inner-product expression `nf() & (-devReff()…)`;
  - the ×2.5 first-step shear error;
  - "about 15 further sites" (this stays in App. D);
  - the per-replicate builds (GCC 11.4/24 ranks, v2412 -O2/128 ranks, GCC 13.3/32 ranks) and c_p values 4.598/4.594/4.597 (now in the provenance appendix and paper-data);
  - the caption builds.
- Defined: MVQN.
- Todo: rewritten as pending results without "earlier meshes" or branch names. The "backward-scheme instability" is kept, as "BDF2 instability".
- Lesson: rewritten without "after correction".

### 5.3.3 Mok cavity
- Rewritten: the 128² diagnosis ("With the baseline settings diverges … The cause is…") is now a statement of the set-up used for 128² and why, with a pointer to App. D, which already had it.
- Removed: the build in the coupling-test item, and "PR #700" in the todo.
- Defined: variational multiscale stabilisation (replaces "VMS-stabilised").

### 5.4 Provisional cases
- Removed: the stale sentence "for the Hessenthaler flap the finest-mesh solution … is still running". The rF2 results are in Table 17.
- **Flexible-bottom cavity:**
  - The todo is rewritten as pending items. It now includes how the COMSOL F_y is integrated, the open question in `cavity_plate_material_law_findings.md`.
  - Removed: "which an earlier version of the case set-up recorded".
  - The iteration and tolerance change for the finest level is now stated as the set-up.
  - The skew-correction GCI finding is stated once, in *Comparison*, with a pointer to App. D.
  - Removed: the COMSOL "solver-artefact plateau of an earlier sequence".
  - The uniform-inflow statement is made once. It is now tentative ("may account for"), with `\hl{[pending: COMSOL sequence with the parabolic inflow]}`.
- **Beam:** "its four-level sequence" is now "our four-level sequence", which fixes an ambiguity.
- **Hessenthaler:**
  - The μ* history ("chosen before the fourth level was added") is now "the GCI contains 64.2 kPa".
  - The convection-scheme instability story moved to App. D (Robustness). The main text keeps the outcome: linear-upwind diverges at the tip corners, upwind gives −17%, and the blend stays bounded.
  - The feature-edge/chamfer detail moved to App. D (Mesh families).
  - The r1-floor explanation is shortened, with a pointer to App. D.
  - Defined: MRI.
- **Blob:** "baseline set-up of the case" is now "the solids4foam tutorial"; "It now has" is now "It has".

## Section 6 (`sectionCrossBenchmark.tex`)

- Todo: the dated history ("first written on 2026-10-10", "entered … on 24 November", "level 3 running") is removed. It keeps the pending Gmsh family, the MMS author decision and the SU2 merge.
- Opening: the list of defects is replaced by a pointer to 6.4, so the list is not repeated.
- **6.1:**
  - Removed: the Womersley "first showed 1.3–1.8" story and "(before, the two were similar…)".
  - Table 19: "before end-flux fix" is now "with a standard mixed end condition".
  - The FSI3 row and paragraphs no longer use the baseline sequence.
  - "Two levels" now uses FSI3 (apparent order 3, not confirmed by 8x′), the 3-D tube and FSI1.
  - Super-nominal examples: the baseline lift order 3.1 and the CFD3 force orders 2.8–3.3 are dropped with the diagnosis.
  - The "disagrees vs not resolved" argument is now held here, for both Featflow and Failer and Richter. The FSI3 Lesson points here.
  - Added: one sentence on the manufactured solution as an example of sub-nominal orders.
- **6.2:** the FSI3 row drops the 8x′ "earlier wall condition" numbers. Item (iii) is marked ring+box.
- **6.3:**
  - The FSI3 noise sentence now points to 5.3.1.
  - The Hessenthaler drift paragraph is tightened ("OpenFOAM-v2412" removed) and is now the single place for the legacy-test drift evidence. Sec. 4.2 and App. D point to it, and Sec. 3.5 no longer repeats the sentence.
- **6.4:** the aliasing mechanism ("reuse of temporary storage … about 15 further sites") is replaced by a pointer to App. D. The defects-table row "first divergence after the fix above" is now "remaining cross-platform divergence".
- **6.5 (Table 23):**
  - The beam row "an independent two-level Kratos solution (three levels)" contradicted itself. It is now "three-level", as in `paper-data/kratos_beam_richter/README.md` (levels 0–2).
  - "earlier calculations solved a different beam position…" is now "a same-lineage calculation appears to use a different plate position".
- **6.6:**
  - The FSI3 solid-share sentence now uses the final cross run ("less than the time-step and coupling sensitivities") instead of the baseline 5%. It is marked ring+box.
  - The cost rows for the FSI3 "uniform family" and the 3-D tube "mesh variant" are relabelled.
  - Added: `\hl{[pending: cost on the Gmsh family]}`.
  - Hessenthaler: "still open" is now "not established".

## Section 7 (`sectionRecommendedSuite.tex`)
- The todo no longer says "level 3 running, expected …".
- The FSI2/FSI3 row:
  - "mesh-motion history; set-up-induced sub-nominal order" is now "sensitivity of the in-phase fluid load to mesh motion and wall conditions";
  - "three levels (baseline) and four (corrected)" is now "four levels, not asymptotic" (marker kept).
- Core suite: the FSI3 sentence points to 5.3.1 and App. D.
- Mok todo: the diagnosis remark is removed.
- Extended-suite todo: "COMSOL element order" is now the parabolic-inflow sequence.

## Section 8 (Discussion)
- (ii): "two of the three cases" is now named cases (FSI1, FSI3, 3-D tube).
- (v): "two of the four defects" contradicted Table `tab:defects`, which has seven rows. It is now "two of the defects".
- "Exact references are diagnostic": now covers the manufactured solution as well as the two linear cases.
- "Limits of the suite": "Kratos … beam … two levels each" was wrong. It is now three levels for the beam and two for the tube.
- Mok: "turned out to solve" is now "solves".

## Conclusions (`article.tex`)
- MMS bullet: shortened. Same findings; "corrections under review" is now "three unreleased corrections".
- "Two levels" bullet: rewritten without the baseline set-up, with the ring+box marker kept.
- "Defects" bullet: kept as motivation but shortened to categories, with pointers to Sec. 6.4 and App. D. The Introduction contribution is unchanged.

## Appendix A (beam)
- Todo: "PR #546" is now "the interpolation-weight correction of Appendix D".
- "Earlier calculations" is now "A same-lineage calculation". The history of the earlier set-up (plate position, sampling face, inflow, transient start-up) is removed. What remains is the explanation of the literature value: its transverse displacement is consistent with a plate 0.05 m downstream.
- "has since confirmed" is now "has confirmed".

## Appendix D (lessons)
- New paragraph "Flux consistency at truncated boundaries" (Womersley, see above).
- Expanded paragraphs:
  - "Path-dependent mesh motion": the FSI3 baseline and the component isolation;
  - "Wall pressure and wall velocity conditions": the 2x′ A/B numbers and the statement that 8x′ lacks the correction;
  - "Errors that only an exact finite-deformation solution exposes": the MMS candidates and the QM k·∇p attempts;
  - "Robustness": the Hessenthaler convection instability and the restart lesson;
  - "Mesh families": the Hessenthaler chamfer.
- The Womersley sentence moved from "Errors that only…" to the new paragraph.
- "Boundary non-orthogonal corrections": the FSI3 0.37/0.18% numbers are replaced by a pointer to Sec. 3.1, which holds them.
- "Coupling convergence tests": the Hessenthaler drift points to 6.3. "baseline FSI3 sequence" is now "uniform FSI3 family".
- "Unstated benchmark definitions": the uniform-inflow 2% statement is made tentative, with a pending marker.

## Sections 3–4 (consistency fixes only)
- Sec. 3.1, 3.2, 4.2: "corrected FSI3 sequence/set-up" is now "FSI3 sequence". "baseline FSI3 sequence" is now "uniform FSI3 family of the coupling studies".
- Sec. 3.1: "coarsest corrected mesh" is now "coarsest mesh".
- Sec. 3.5: removed the duplicated "the default test can nevertheless accept insufficiently converged steps".
- Sec. 4.2: the steadiness-check sentence points to Sec. 6.3, and its nested brackets are removed.

## Inconsistencies found and fixed
1. Sec. 5.4 said the Hessenthaler rF2 acceleration run was "still running". Table 17 has its result.
2. The MMS paragraph said the effect of the missing moving-wall non-orthogonal correction "has not been quantified". Sec. 3.1 reports it (≤0.37%).
3. Table 23, beam row: "two-level Kratos solution (three levels)". It is three levels.
4. Discussion, Limits: the Kratos beam had "two levels". It has three.
5. Discussion (v): "two of the four defects". Table `tab:defects` lists seven.
6. FSI3, comparison with Featflow: "(at 4x′: 3.2%, 6.6% and 7.0%)" compared Featflow with the 4x′ level computed with the uncorrected wall derivative, not with the 4x′ of the table. The table's 4x′ gives 2.3%, 4.9% and 5.8%. The parenthetical is removed rather than replaced. Sec. 6.1 had the same pair (3.2%→2.2%, 6.6%→4.5%); it is removed too.

No other number was changed.

## Pending markers kept
- `\hl{[FSI3: ring+box, to be replaced]}`:
  - abstract and Conclusions;
  - Table 2 footnote;
  - 5.3.1: the evidence paragraph, the table caption, three items and the Lesson;
  - 6.1 (three), 6.2 (table and text), 6.6;
  - Table 19;
  - Sec. 7: the tier table and the extended suite.
- `\hl{[pending merge of the SU2 data]}`: all occurrences kept.
- COMSOL markers:
  - `\hl{[pending: COMSOL sequence in paper-data; how its $F_y$ is integrated]}`;
  - `\hl{[pending: COMSOL sequence with the parabolic inflow]}` (5.4.1, Table 23, App. D).
- New markers:
  - `\hl{[pending: re-check of the parallel run]}` (collapsible channel);
  - `\hl{[pending: cost on the Gmsh family]}` (Table 24).
- Todos kept, with the history stripped:
  - MMS (release, novelty, A/B, open cause);
  - FSI3 (to be replaced; the 8x′ decision);
  - Turek (2006) summary values;
  - FSI2/Featflow CSM/CFD citation;
  - 3-D tube;
  - Mok;
  - cavity;
  - Hessenthaler and blob status;
  - cross-benchmark;
  - cost;
  - suite tiers;
  - beam A/B.

## Unsure, for Philip
- The ring+box Q_in replay values were not carried over: the corrected treatments on two families, the 4–5° rotation and the 6% family difference. Only the −2.7% 4x′→8x′ consequence is kept. If you want the prescribed-motion study as evidence that the remaining FSI3 error is on the fluid side, it could go to the supplement.
- The "uniform family" coupling studies (force noise, RN vs IQN-ILS) stay in 5.3.1, although that family's displacement results are now only in App. D.
- The FEBio "source inspection found no defect" sentence was removed (contributor debugging).
- I kept the Hessenthaler preliminary drift in 6.3 as the evidence for the strict coupling test. It is a finding about a coupling test, but it comes from a preliminary run.

## Corrections after an independent review of the diff
- The FSI3 sequence statements (Sec. 5.3.1, 6.1, Table 19, Conclusions) now say explicitly that the 4x′→8x′ change of 1–2% compares two levels computed with the uncorrected wall-normal velocity derivative. The apparent order of about three comes from the corrected 1x′–4x′ levels.
- The prescribed-motion load is stated as using the coupled 2x trajectory of the uniform family, not 2x′. The 6% time-step share is attached to that load change, as in the source (§20.1), not to the coupled displacement change.
- The 8x′ effect of the wall derivative (about 0.4%) is stated as expected, not measured.
- The cross run (solid and fluid shares) is qualified as computed with the uncorrected derivative, in Sec. 5.3.1 and 6.6.
- App. D: orders of 1.1–1.3 are stated only for the treatments run on three levels. The wall velocity condition and the mesh-motion diffusivity are two-level results.
