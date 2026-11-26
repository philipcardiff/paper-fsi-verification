# Response to Philip's notes 4 (Section 5 restructure; 3-D tube, Mok, cavity; Sections 6–8)

Branch `text/notes4-restructure`. There is one entry per comment in `philip_notes_4.md`, plus the decisions taken beforehand. Unless stated otherwise, the edits are in `sectionBenchmarkSuite.tex`. The PDF goes from 100 to 101 pages (`latexmk -pdf article.tex`). The build has no errors, undefined references or multiply defined labels.

## Restructure of Section 5 (decisions taken beforehand)

- **Order by reference category.** The cases now appear in this order:
  - **5.1 Exact references (A):** ring, Womersley, manufactured solution.
  - **5.2 Independently converged references (B):** collapsible channel, Turek–Hron family, Mok cavity, flexible-bottom cavity.
  - **5.3 Partly quantified references (C):** beam in cross-flow.
  - **5.4 Unquantified references (D):** 3-D tube, Hessenthaler flap, blob in treacle.
- **Turek–Hron FSI2.** FSI2 is category C, but it stays with its family in 5.2, and both 5.2.2 and 5.3 say so. Splitting the family would repeat the geometry, definitions and references.
- **Labels.**
  - Two new labels: `sec:partly_quantified_refs` and `sec:unquantified_refs`.
  - `sec:classical` and `sec:provisional` are removed.
  - Every `\ref` to them is updated. In Section 4, the list of studies that use the all-measure test now says "the other studies of Section 5".
- **Section 5.4 ("Provisional") removed.**
  - The cavity is now a full category-B case.
  - The beam is a compact category-C case in the main text, with the detail kept in Appendix A, which is retitled "The beam in cross-flow: published inflow and sequences".
  - The Hessenthaler flap and the blob in treacle are category-D cases.
  - All "provisional/mature/[PROVISIONAL]" status wording and the `$^\ast$` provisional markers are removed: Section 5, Tables 2, 18 (orders), 19 (QoI dependence), 23 (cost), 24–25 (suite tiers, drivers) and the provenance table, Appendix E, and the Conclusions.
  - Kept: "provisional indicator" in Section 4 and Table E (FSI3 GCI). This is a methodology term for a GCI outside the asymptotic range, not a case status. Rename it if you prefer (for example "unconfirmed indicator").
- **Table 1** now lists the flexible-bottom cavity under B (SU2–CalculiX and FEniCSx) and the beam under C, with the independent Kratos check. D now also covers "independent but not shown to converge" (blob).
- **Counts.**
  - The count of category-B references among the ten non-exact problems is now five, of which two are published values.
  - It is updated in Section 6.5, the Discussion and the Conclusions.
- **Section 5 preamble.**
  - It states the category order.
  - It states that every case has a self-contained Problem paragraph.
  - It names the solver "solids4foam", with the coupling named where it matters.
  - It defines the L1, L2, … level convention once.

## 5.3.2 Three-dimensional pressure-pulse tube

**Present it later (after B and C)?**
Done. It is now the first category-D case (5.4.1).

**Hookean vs St Venant–Kirchhoff.**
- The new Problem paragraph states that published calculations differ in the wall law. It says the three references use a St Venant–Kirchhoff continuum, and that the Hookean variant is the one studied here, by solids4foam and by the independent Kratos solution.
- At ≈3 % hoop strain the two laws differ in stress by a relative amount of the order of the strain. That could exceed the L2→L3 changes (0.3–1.0 %).
- paper-data does not quantify the difference: `3dtube_level3_findings.md` lists it as "not quantified here". It is therefore marked `\hl{[pending: St Venant–Kirchhoff wall run]}`.
- Suggestion: an L1 St Venant–Kirchhoff run is cheap (2–8 min on a workstation, from the cost table) and would settle this.

**AMI.**
"What it tests" now says that the interface could be meshed conformally, but a non-conformal (AMI) transfer is used deliberately to exercise it.

**"hence, this reference is categorised as …".**
Applied to every case Reference paragraph: ring, Womersley, manufactured solution, collapsible, FSI1/FSI2/FSI3, Mok, cavity, beam, tube, Hessenthaler and blob. Each gives the reason.

**Terms before use (L3 etc.).**
- The preamble defines the level convention.
- The tube now has a Numerical model paragraph, with the L1–L3 meshes and steps, before its Results.
- The collapsible Reference paragraph used "L4" before the meshes were described. It now reads "the finest fluid level, L4, below".
- I checked the other cases. Their mesh names appear in the Numerical model or the preamble before use.

**Fig. 8 (now Fig. 10) "smoothed and resampled" — why?**
From the reference CSV headers and the solids4foam commit that recorded them:
- The two BDF1 curves were read from raster figures by colour masking.
- A raster line is several pixels thick, and its centre is quantised to whole pixels, so the trace steps by one pixel (≈ the stated ±3 µm) between columns.
- The trace was therefore averaged over five pixel columns and resampled to the published time step, 0.1 ms.

The raw pixel traces were not kept, so the raw points cannot be plotted and the figure is unchanged.

Instead, the text now:
- explains the procedure;
- bounds the smoothing error by applying the same smoothing to our L3 history. Moving averages over 0.2 and 0.3 ms lower the radial peak by 0.4 and 0.9 µm and raise the axial trough by 0.3 and 0.6 µm. That is a tenth to a third of the ±3 µm reading accuracy, so the smoothing error is smaller than the digitisation error.

The caption now refers to this explanation. The one assumption, that five pixel columns span less than 0.3 ms (i.e. the 20 ms axis is wider than about 340 px), is marked `\hl{[check: pixel width of the two published figures]}`. Re-digitising from the PDFs would close this. The Lozovskiy and Eken figures are not in the literature folder.

**Figure placement.**
- The solution-gallery figure was `[p]`. Because LaTeX keeps floats in order, it held every later figure of Section 5 back to float pages at the end of the section (pp. 47–54, up to 20 pages from their text).
- It is now `[tbp]`, and the float parameters are relaxed in the preamble (`topfraction` 0.9, `textfraction` 0.1, `floatpagefraction` 0.75, up to 3 top / 5 total floats).
- `placeins` with `[section]` was already loaded.
- Every Section 5 figure and table now sits on the first page of its case or within three pages of it.
- Page count: 100 before, 101 after. The restructure added about 1.5 pages (problem definitions, beam case, smoothing paragraph), and the placement change recovered about 1.

**Appendix tables mentioned first in Reference/Results.**
Moved. Each case now gives its own results first and points to its Appendix E table at the end of Results: ring, Womersley, manufactured solution, collapsible, Turek–Hron, Mok, cavity, beam, tube, Hessenthaler and blob.

**Run L4 for the tube?**
Not run.
- `3dtube_level3_findings.md` estimates L4 at ≈8 M fluid cells, ≈16× L3 and ≈10–20 k core-hours.
- It is the only route to a convergence statement on u_r,max, so it is kept in the todo, alongside a third Kratos level.
- If the time-step study now being rerun shows that the time error is as small as before, L4 would be scaled.

**"Measured on a variant of the meshes …".**
- The Time item now reads `\hl{[pending: time-step study on the corrected meshes]}`. The mixed-mesh framing is removed from the main text and from the Table 14 caption ("Wall arcs defined to full precision").
- The earlier mesh variant is mentioned only in Appendix D (the platform replicates), in the cost table and in the provenance footnote.
- The time-step numbers in Section 6.1 and Table 18 are marked with the same pending flag.

**"Platforms. This case exposed two defects …".**
- Moved to Appendix D, into the aliasing paragraph: it now reports the two- and three-platform agreement after both corrections, and that these replicates used the earlier mesh variant.
- The tube's Lessons keep one sentence: "Replicating the case on several platforms exposed two defects in solids4foam that no single-platform refinement study would have shown (Appendix D)."

**"Implicit Euler" → BDF1.**
- First use (Section 3.1): "first-order backward Euler (BDF1)".
- Thereafter BDF1: the tube Reference, figure caption, literature-discrepancy item, Kratos paragraph, Table 22 and the Discussion. BDF2 is unchanged.
- The legend of Fig. 10 now reads "BDF1" as well. It was regenerated with `paper-data/3dtube_level3_corrected_ecore/figure/tube_histories_figure.py` and `paper-data/figstyle`, with the data unchanged.

**"wall model, pressure switch-off".**
The Reference paragraph now states:
- All three references model the wall as a three-dimensional St Venant–Kirchhoff continuum, not as a shell.
- The published definition (Fernández & Moubachir) prescribes the inlet pressure for 3 ms, and the references state no ramp at its switch-off, whereas here it is ramped to zero over 0.1 ms.
- Their meshes differ.

The literature-discrepancy item now says that these remaining differences are not quantified. Source: `3dtube_level3_findings.md` §3 and §8 ("StVK vs Hookean, instantaneous vs ramped switch-off").

## 5.3.3 Mok lid-driven cavity

**"The reference taken here is …".**
Done.

**Problem definition missing.**
- Added a Problem paragraph with the domain, fluid properties, lid velocity 1 − cos(2πt/5), openings (left inflow v̄(8y − 7) over 0.875 ≤ y ≤ 1, right p = 0, fixed in space), no-slip walls, membrane (0.002 m, E = 250 Pa, ν = 0, ρ = 500, St Venant–Kirchhoff, clamped, traction-free underneath), start from rest, end time 70 s, and the QoI window.
- The Kratos paragraph now gives only Kratos particulars.
- I checked every case. Problem definitions were added or completed for:
  - Womersley: dimensional E, ν, ρ, μ, f and pressure amplitude;
  - collapsible: geometry, E, ν, fluid, inlet drive and p_ext(t);
  - Turek–Hron: a new Problem paragraph with ρ_S, μ_S, ν_S, Ū and the inflow ramp for FSI1–3;
  - flexible-bottom cavity: geometry, plate size, boundary conditions;
  - beam: duct, plate, inflow, outlet;
  - 3-D tube: geometry, properties, pulse, ramp, QoIs;
  - blob: geometry, inflow g(t), properties.

  Ring, manufactured solution and Hessenthaler already had complete ones.

**Naming our solver.**
- "solids4foam" is used throughout, with the scheme named where it matters. The Section 5 preamble says so.
- "the code under test" is replaced in every file: by "solids4foam" where it means our code, and by "the code being verified" in the generic statements of Section 1 and Section 4 ("The code being verified is the compiled executable").
- "the present code" is gone as well.

**"Kratos needs no such change".**
Replaced by a qualification: "Kratos solves the velocity and pressure of its stabilised fluid together, without a segregated outer loop, so this change has no counterpart there." This now sits in a new Numerical model paragraph.

**"legacy coupling test" — define without "legacy".**
- Section 3.5 (Convergence test) now names the two tests by what they require:
  - the **either-measure test**, the solids4foam default, accepts when min(r1, r2) is below the tolerance;
  - the **all-measure test** requires every measure, i.e. max(r1, r2).
- "legacy" and "strict" are replaced everywhere in the paper (Sections 3–7, Appendices A, D and E, captions and tables). "Strictly periodic" is untouched.
- The FEniCSx supplements, which are separate documents, still use the old words.

**Strict-test stall author note — resolve.**
Resolved, from `drafts/mok_pr760_rerun.md`. The Coupling-test item now says:
- About 98 % of the steps were accepted on one residual measure only.
- The all-measure test on 32² and 64² changes the QoIs by at most 0.025 %.
- Rerunning 16²–96² with the later build, in which the boundary non-orthogonal correction is restored, changes them by at most 0.041 %, and every mesh change by < 0.05 percentage points. The `\hl{[pending: 128² relaxed run with the later build]}` flag covers the 128² run.
- The all-measure test was not used for the reported values because it stalls on some meshes later in the run, at t ≈ 14–65 s depending on the mesh, either far above the tolerance (16², unrelaxed 128²) or just above it (96²).
- The todo is reduced to the IQN-ILS/Aitken difference and the tutorial decisions.

Pushback on the wording of the decision note:
- The 0.041 % comes from rerunning 16²–96² with the either-measure test on the new build. It is not an all-measure rerun of 16²–96²: with the new build the all-measure test completes only on 32² and 64², where it changes the QoIs by ≤ 0.025 %. The text states it that way.
- The stall times run to 64.6 s (16², no mode reuse), so the text says 14–65 s, not 14–60 s.

## 5.4 Provisional benchmarks (now 5.2.4 cavity, 5.3.1 beam, 5.4.2 Hessenthaler, 5.4.3 blob)

**Tuković's constitutive law.**
- The text now says that a neo-Hookean law was not available in that code at the time (Ž. Tuković, private communication). The published plate was therefore St Venant–Kirchhoff or linear elastic, and its 4.3 % larger displacement is not explained here: `\hl{[pending: linear-elastic check]}`.
- The earlier inference "neo-Hookean 1.047–1.052 against 1.043" is removed, together with its entry in Table 22.

**"so the monitoring point must be stated".**
Rephrased as asked: "The through-thickness position of the monitoring point changes u_y measurably: the mid-thickness value differs … by 0.5 %, and … underside … 1.2 % smaller …, so point A is specified here as the point on the top surface."

**"of the code under test".**
Replaced: Table 13's "Relation: to solids4foam", and everywhere else as above.

**FEniCSx independence.**
- Agreed. FEniCSx and SU2–CalculiX are now both treated as independent solutions: Table 13, the text, Table 1, Section 6.5, the Discussion and Appendix D.
- SU2–CalculiX is also noted as computed blind.
- The quasi-monolithic solution stays same-lineage.
- The two independent limits agree to 0.006 % (u_y) and 0.0003 % (F_y). They form the category-B reference, with uncertainties of ≈0.05 % (u_y) and 0.001 % (F_y).

Consequence for the blob: FEniCSx is also an independent implementation there. The blob stays D, as you decided, because its FEniCSx displacement converges at only first order, and the text now says exactly that. I added "or independent but not shown to converge" to the D definition in Table 1.

**"the independent solutions expose an interpolation error …".**
- The Results paragraph now states only the final skew-corrected results and their agreement with both independent solutions.
- The lesson is in the Lessons paragraph as one sentence pointing to Appendix D ("Interpolation of the continuity flux"), which keeps the detail and is updated for two independent solutions.

**COMSOL "by E. Tandis".**
Now "computed by a co-author". The ANSYS solution is treated the same way. The contributor names remain only in the author-note todo for pending COMSOL/ANSYS items.

**Blob in treacle.**
Kept, with `\hl{[pending: second independent solution]}` under the heading. A full Problem paragraph was added.

**Hessenthaler.**
Rewritten shorter, as an honest mesh-sensitivity and MRI-validation case.
- The Numerical model, Results and validation text are condensed. The three tables are kept.
- What can be claimed: on L3 the tip is 16.73 mm with either pressure condition, the last mesh change is ≤ 1.0 %, and the quasi-monolithic solution is within 0.09 mm.
- What cannot be claimed: an error bound, and a judgement of the 2 % discrepancy from the measurement, which is about twice the last mesh change and has no experimental uncertainty.

**Beam in cross-flow.**
- A compact main-text case (5.3.1) with Problem, Reference (inflow correction and category C), Numerical model, Results, Kratos check and Lessons.
- Appendix A keeps the evidence, the full sequences (Table A1) and the Kratos detail. Its "kept provisional" and "[author decision: category]" sentences are removed.
- The L5 statement is flagged `\hl{[pending: L5 or author decision]}`.

## Sections 6, 7 and 8 — merge?

Not merged in this PR, as instructed. Proposed outline: a single section, "Synthesis and recommendations", replacing Sections 6–8, which now run from p. 53 to p. 65 (about 13 pages).

1. **Convergence against references**: today's 6.1 and 6.2, merged.
   - Table 18 (orders) gains a "slowest QoI" column, and Table 19 (QoI dependence) is dropped.
   - Keep the three observations (sub-nominal orders, two levels are not enough, keep "disagrees" apart from "not resolved").
   - Cut the FSI3 crossing narrative, which is told three times in Sections 5 and 6.
2. **Iterative error and coupling tests**: today's 6.3. Keep Table 20 (tolerance sweeps) and the Hessenthaler drift evidence, shortened to one paragraph, and the "internal agreement is not external verification" point.
3. **Published references and model-form differences**: today's 6.5, plus Discussion "Model-form differences" and "Recommendations to benchmark authors". Keep Table 22 (published issues). Shorten the SU2 artificial-compressibility paragraph to two sentences.
4. **Defects and what a verification claim should contain**: today's 6.4 observations, plus Discussion "What a verification claim should contain" and "Exact references are diagnostic". Table 21 (defects) moves to Appendix D, which already describes each defect.
5. **Cost and a tiered suite**: today's 6.6 plus 7.1–7.2.
   - Keep Table 24 (tiers).
   - Move Table 23 (cost) to an appendix and keep one paragraph on the cost/evidence trade-off.
   - The FEBio cost comparison goes to Appendix B.
6. **Verification and validation; limits of the suite**: today's Discussion paragraphs, unchanged in substance.
7. **Infrastructure** (7.3): reduce to a short paragraph pointing to the Data Availability statement. Table 25 (driver properties) moves to Appendix F (provenance).

Estimated effect:
- About 4–5 pages fewer in the main body: Tables 19, 21, 23 and 25 out of the main text (≈3 pages), and the removal of repeated numbers and narrative (≈1.5–2 pages).
- About 2–3 pages fewer overall, since three of the tables move to appendices.
- The section numbering of the Conclusions would become 7.
