# Response to Philip's notes on Sections 5.3 and 5.3.1 (Turek–Hron FSI1, FSI2, FSI3)

Branch `text/section531-review`. Each entry paraphrases one comment in `philip_notes_5.3_5.3.1.md` and says what was done and where, or why not. Unless stated otherwise, the edits are in `sectionBenchmarkSuite.tex`, Sections 5.3 and 5.3.1. The PDF goes from 100 to 99 pages, and the subsection from about 3 970 to 3 470 words. This is despite a new FSI1 table and the new explanations asked for below.

## "softer" → "more flexible"

Done: "heavier and more flexible flag". This was the only place in the paper where "softer" stood for stiffness.

## Re-defining categories A/B/C; self-descriptive names

- **Re-definitions.** The categories are defined once, in Table 1 (Section 4). In 5.3.1 the first use of each category now has a short gloss: "category B: independent and converged" and "category C: published, with partly quantified error". Later uses give the letter only. The other category uses in the paper are labels, not re-definitions, so they were left alone.
- **Names.** Not renamed yet, as you asked. Proposals, checked against the Table 1 definitions:
  1. **Exact / Independent-converged / Published / Weak.** "Weak" covers what C-weak actually lists: digitised, single-mesh, mutually inconsistent, preprint-only or same-lineage values.
  2. **Analytical / Verified-numerical / Published-numerical / Indicative.**
  3. **Exact / Converged / Partly-quantified / Unquantified.**
- **A caution on "Self-convergence-only" for C-weak.** Section 4 says self-convergence is the primary evidence for C *and* C-weak. A name based on the evidence used would therefore not tell C and C-weak apart. Names that describe the reference itself, as in the sets above, fit the definition "the category describes the reference, not the benchmark problem and not the quality of the present calculations".
- **Category B in sets 2 and 3.** "Verified" or "Converged" alone hides the independence requirement, which is half of the B definition. Set 1 keeps it.

## "the definitions are the same": which definitions?

The benchmark definitions are now stated once, in What it tests. In each period of `u_y`:
- mean = (max+min)/2 and amplitude = (max−min)/2;
- these are averaged over the analysis window;
- the frequency is that of `u_y`.

The Failer & Richter sentence now reads: "Its means and amplitudes are defined as above but averaged over t = 9–10 s rather than over the last five periods before t = 6 s."

## Tone of the opening of 5.3

The Section 5.3 opening is rewritten:
- It starts with the compliment: these are the most widely used FSI benchmarks, and the Featflow multi-level results are among the most carefully documented FSI references available.
- It then states the question: how far such a published multi-level result, whose own discretisation error is only partly quantified, can bound the error of another code.

The substance is unchanged. The later compliment sentence in 5.3.1 is dropped to avoid repetition.

## "≈4 N/m per level"

Now "the lift amplitude rises by 2.7% per level without slowing". The level 2 → 3 → 4 values are 146.00, 149.91 and 153.91 N/m, so each step is +2.7%.

## Featflow per-period lift drift: monotone?

**Not monotone, so the cumulative-error sentence was not added.**

I downloaded the published level-4 history (`ref_fsi3.point`, Δt = 0.00025 s, t = 5–6.44 s) from the Featflow site and computed the per-period lift amplitude. Periods run between successive upward crossings of the `u_y` mean:

157.74, 154.84, 153.09, 153.82, 153.77, 153.72, 153.53 N/m

The amplitude falls by 2.9% over the first two periods. It then stays within 153.1–153.8 N/m (±0.2%), with no trend. The drag amplitude declines more slowly, from 27.81 to 27.45 N/m (1.3%) over the first five periods, and then levels off (27.46, 27.51). It is not part of the text.

The text now says:
- the amplitude falls from 157.7 to 153.1 N/m over the first two of its seven full periods and then stays within 153.1–153.8 N/m;
- the tabulated amplitude therefore depends on the window by up to 2.7%;
- the variation is not monotone, its cause is not known, and it may be a transient that decays within about two periods.

No cause is attributed to Featflow. The Appendix D lesson on path-dependent mesh motion is not invoked: an accumulating error would keep drifting, and this history does not. The per-period values are recorded in a source comment next to the sentence.

In Section 6, "drifts by 2.7%" is changed to "the per-period lift amplitude varies by 2.7%" in two places: the paragraph on the Featflow table and Table 22.

## Which Featflow values we use, and which are more trustworthy

The paragraph is rewritten to answer this directly:
- **Which we use.** We use the detailed level-by-level tables, which carry convergence information and published time histories. We do not use the rounded summary values.
- **Which is more trustworthy.** The detailed tables, for three reasons:
  1. their time steps are smaller;
  2. in the chapter itself, the level-4 FSI3 drag amplitude changes by 13% when Δt is halved (20.00 → 22.66 N/m), so the summary value is not converged in time;
  3. the independent Failer & Richter sequence agrees with the table (27.91 against 27.74 N/m), not with the summary value (22.66 N/m).
- **Why the 18% matters.** A code compared against the summary values would carry an 18% difference in the FSI3 drag amplitude that belongs to the reference. Any comparison must therefore state which set it uses.
- **The confusing "driver reproduces…" sentence.** Now: "Applied to the published Featflow time histories, our post-processing reproduces Featflow's own tabulated statistics to within 0.6% for FSI3 (0.8% for the drag amplitude, 0.3% for the frequencies) and 0.3% for FSI2, which validates our implementation of the definitions above."

## The author note on the Springer chapter

Resolved from the PDF, and the note is removed.

- The summary values are the "ref." rows of the chapter's Table 14 (FSI2) and Table 15 (FSI3).
- Each is **level 4 (15 872 elements) at the smaller of the chapter's two time steps**: Δt = 0.0005 s for FSI3 (pair 0.001/0.0005) and Δt = 0.001 s for FSI2 (pair 0.002/0.001).
- The Featflow pages keep these values only inside HTML comments marked "old values". Their visible tables are later computations, down to Δt = 0.00025 s (FSI3) and 0.0005 s (FSI2).

This is now stated in the text, with the source recorded in a LaTeX comment.

**Who quotes which.** I did not survey the later literature, so the text names no papers. It says only that the chapter's values are the frequently quoted ones, as Section 8 already says ("often the summary values of Turek and Hron"). If you want named examples, that needs a short citation check.

## Table 8: purpose

Kept in the main text, because it supports a needed point. The 18% is cited again in Section 6 (Table 22 and the paragraph after it) and in Appendix D.
- The caption now opens with the point: "the two published Featflow sets differ by up to 18%, so a comparison must state which it uses".
- It gives the level and time step of both sets.
- The rows are reordered FSI2, then FSI3.
- The text cites it at the sentence that makes the point.

## FSI3 numerical model

- **Why the unstructured mesh.** Added: the Gmsh family replaces the block-structured family, whose neighbour cell-size ratios of 2–2.8 through the wake persist under refinement and spoil the observed orders (Appendix D, "Mesh families and geometry").
- **"PIMPLE 3/3/1".** Spelled out: each fluid time step takes three outer (PIMPLE) correctors, each with three pressure correctors and one non-orthogonal corrector. Checked against `system/fluid/fvSolution` of the L1–L3 cases: `nOuterCorrectors 3`, `nCorrectors 3`, `nNonOrthogonalCorrectors 1`.
- **Inconsistency with the uniform family.**
  - What was run: `paper-data/hronturek_fsi3_gmsh/README.md` confirms that the Gmsh family is the FSI3 result. On it, there is no coupling-tolerance sweep, coupling-algorithm comparison or separate time-step study.
  - The tolerance and Robin–Neumann/IQN-ILS studies exist only on the uniform block-structured family, which has a path-dependent mesh motion and a zero-gradient wall pressure.
  - Removing them would break references in Section 4 (Robin–Neumann vs IQN-ILS), Section 6 (tolerance noise; cost table), Table 2 and Appendix D. I therefore kept them, once, as an explicitly labelled paragraph, "FSI3 coupling sensitivity, on a separate family".
  - That paragraph says they were run before the Gmsh family existed and that the family's displacements are superseded. They remain informative for the force amplitudes because the iterative error depends mainly on the coupling tolerance and the time step.
  - The uniform family is now defined once, in Numerical model, as the FSI1/FSI2 family.
- **Recommendation: a cheap coupling check on Gmsh L1.** Worth running, to replace this paragraph. One L1 run (96 712 cells) took 9.7 h on 24 ranks. Two cold runs from t = 0 (restarts are not possible for this solid):
  1. tolerance 10⁻⁷ at the production Δt: about 12–15 h on 24 ranks, given the extra iterations;
  2. optionally, Δt/2 at 10⁻⁶: about 20 h.

  If the QoIs change by much less than the L1→L2 changes, the separate-family paragraph could shrink to one sentence or go to Appendix D. Nothing was run.

## "relative errors" in the FSI2 table

Now "differences", in the caption and the text. The Featflow FSI2 table is category C, with level 3→4 changes of 0.9–2.5%. This matches the earlier oomph-lib choice for the collapsible channel. The FSI1 table also says "difference from Featflow".

## Order FSI1, FSI2, FSI3

Done throughout 5.3.1: What it tests, Numerical model, tables, Results paragraphs, Lessons and Table 8.
- FSI2 keeps its current content.
- It is marked `\hl{[pending: FSI2 on the Gmsh family]}` in Numerical model, the Results paragraph, the FSI2 table caption and Lessons.
- The FSI2 part of the earlier author note ("whether to compute an FSI2 4x level") is dropped, since the rerun supersedes it.

## FSI1: only the high-order results

Done, with one caveat.
- **Main text.** FSI1 now presents the cubic-solid sequence (1x, 2x, 4x) in a new Table 9 (`tab:fsi1_ho`, from `paper-data/fsi1_high_order/`):
  - monotone from below;
  - observed orders 1.3–1.5 (`u_x`, `u_y`, drag) and 2.1 (lift, marked indicative because the 4x lift varies by ±0.2% over the closing 20 s);
  - 4x differences from Featflow of at most 0.48%.
- **The standard solid** gets two sentences: its bending stiffness dominates the error on coarse meshes (`u_y` −14.2% against −1.7% on 1x), and the cubic solid on 4x is closer in `u_y` than the standard solid on 8x.
- **Caveat.** The four-level, strict-test, standard-solid study is the better-controlled sequence:
  - it has four levels against three;
  - it uses one pseudo-time step and the strict test, whereas the high-order runs use the driver's legacy test and pseudo-time steps of 0.025, 0.025 and 0.0125 s;
  - it supplies the FSI1 values in Appendix E and the baseline for the CalculiX comparison.

  I therefore kept it, but moved it out of the main text. Its table is merged with the Appendix E reference-values table, which now carries both labels (`tab:fsi1`, `tab:refvals_fsi1`). The remaining details (pseudo-time-step check, steadiness, IQN-ILS/Robin–Neumann) are two sentences there. The main text keeps one sentence pointing to it.

## Independent numerical comparisons: shorten

Cut from about 520 to about 190 words:
- **FSI1, SU2–CalculiX.** One sentence on the fully independent solution: `u_y` and lift Richardson values within 0.04% and 0.03% of Featflow, and `u_x` and drag too slowly converging to discriminate.
- **FSI1, OpenFOAM–CalculiX.** One sentence on the shared-fluid solution: 0.02% from our standard-solid Richardson values.
- **FSI3, preCICE.** One sentence on the cross-check:
  - agreement on the same coarse mesh to 1.1% in the displacements and 1.0% in the frequency and drag mean;
  - 6% and 15% in the force amplitudes, which one level cannot attribute;
  - not a comparison for the Gmsh sequence.

The FSI3 sentence moved here from the time-step paragraph.

## Other changes made while editing

- **FSI3 comparison with the references.** Shortened: the L3 differences are in Table 11, and the text keeps only their interpretation.
- **FSI3 Lessons.** The two Lessons sentences on FSI3 are merged.
- **Category footnote.** The "(category B)" footnote of the FSI3 table is dropped.
- **Build.** `latexmk -pdf article.tex`: no errors, no undefined references or citations. 99 pages, down from 100.
