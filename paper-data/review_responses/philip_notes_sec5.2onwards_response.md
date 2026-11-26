# Response to Philip's notes on Section 5.2 (collapsible channel) and on reference values

Branch `text/section52-review`. Each entry paraphrases one comment in `philip_notes_sec5.2onwards.md`, says what was done and where, or why not. The notes refer to the draft before the headings were standardised. They are mapped onto the current headings (What it tests / Problem / Reference / Numerical model / Results / Independent numerical comparisons / Lessons). Unless stated otherwise, all edits are in `sectionBenchmarkSuite.tex`, Section 5.2.1.

## Section 5.2.1, collapsible channel

**"channel width" → "channel height".** Done, in both places. While checking, I also fixed "collapses ... *to* ≈21% of the channel width". The peak displacement is 0.216 m in a channel of height 1 m, so the wall collapses *by* ≈21% of the height, not *to* it.

**"Re St = 50": St or the product?** It is the product, and the value is right. oomph-lib's driver (`collapsible_channel_ref.cc`) defines `ReSt = 50` as "Womersley = Reynolds times Strouhal", the coefficient of ∂u/∂t in its non-dimensional equations. The case JSON gives `Re = 50, ReSt = 50`. Time is scaled on the advection time a/U, so St = 1 and Re St = Re. The Problem paragraph now says this explicitly: time is scaled on the advection time, so the Strouhal number is one and the Womersley parameter Re St, the product that multiplies the time derivative, is also 50.

**Jensen & Heil sentence: what they did and what we changed.** Rewritten as separate sentences.
- What Jensen & Heil did: in their problem, as implemented in the oomph-lib tutorial, the wall is a massless Kirchhoff–Love beam with h/a = 1/100 and an axial pre-stress of 10³ times its effective modulus, and the external pressure is applied impulsively.
- What we changed, and why. Each change gets one sentence with its reason:
  - no pre-stress, because a continuum cannot carry the pre-strain;
  - h/a = 1/20 (the value of Heil et al. 2008), with the pressure scaled;
  - clamped ends and a 0.25 s ramp, because a continuum fixed on its end faces is clamped and the impulsive start limits second-order temporal convergence;
  - ρS = ρF, because partitioned coupling fails without it at small steps and on fine meshes. This shifts the response by 0.28% of the peak.

Source: the tutorial README ("Relation to the published case") and the oomph-lib driver.

**"purposely not reproduced".** Done: "The published values are therefore purposely not reproduced; the reference below is computed for the modified problem."

**"the code under test".** Now "the solids4foam finite volume solver used here". Within Section 5.2.1, "not a property of the code under test" is also now "not a property of the solids4foam solution".

**0.06% of what; "open-loop".** Rewritten in Reference:
- "Open-loop" is now described as a replay. The fluid load recorded at every step of one coupled run is applied as a prescribed load, without feedback from the wall motion to the fluid, to the continuum wall alone and to the beam alone.
- The two first troughs differ by about 0.06% of the reference peak displacement (0.057% at the wall midpoint with 320 beam elements).
- The beam trough still deepens at about first order, so this is a lower estimate, not a bound.

Numbers are from the FEniCSx supplement (MF1, Table `tab:cc_off`).

**FEniCSx mentioned with no details.** The Independent numerical comparisons paragraph now says what the implementation is and what it shows:
- a monolithic ALE finite element code (quadratic velocity and displacement, linear pressure), written and run by the authors, sharing no code with solids4foam or oomph-lib, on three meshes and four time steps;
- its finest-mesh first trough is within 0.02% of P of extrapolated oomph-lib and within 0.001% of L4, and its history is within 0.075% of P;
- it shows the same 0.5% floor against the tutorial-resolution reference, consistent with that floor being the reference's own spatial error;
- its successive-Δt trough changes agree with oomph-lib's to 0.003–0.014%;
- it provides the open-loop model-form estimate.

The details stay in the Supplementary Material.

**blockMesh clause out of place.** It is moved to a new, short Numerical model paragraph. That paragraph gives the meshes (blockMesh; fluid L1–L4 with their cell counts and ratios; wall 160×8, high-order solid), the baseline Robin–Neumann coupling at Δt = 0.00625 s, and IQN-ILS as the comparison. Previously the case had no Numerical model heading.

**Fig. 6 caption: "error" or "difference"?** Decision: "difference". oomph-lib extrapolated to zero mesh size is the category-B reference, but it has its own uncertainty and a model-form difference of the same size as the finest-level values, so calling the deviation an error would overstate it. The following now use "difference":
- the caption ("maximum difference from oomph-lib");
- the figure's axis label and legend (regenerated);
- the Results bullets;
- the related sentences in Section 6 (tolerance paragraph) and the Conclusions.

The 0.53% "tutorial-resolution error" keeps the word "error": it is oomph-lib's own spatial error against its extrapolated limit.

**Grid names F1–F4.** Now L1–L4 throughout, including the figure. See "Grid names" below.

**"by the same group".** Now "written and run by the authors of this paper and sharing no code with solids4foam or oomph-lib". The same wording problem existed elsewhere, so it was changed there too:
- the cavity text ("computed by the authors of this paper with a separate implementation");
- the cavity table ("same authors");
- the blob Reference;
- the Supplementary Material appendix.

The audit wording (a same-lineage cross-check, not independent validation) is kept: the text calls it a separate implementation, not an independent one.

**Table of values to copy and compare.** Done for every case: see "Reference values" below.

**"error falls to ≈0.5% the reference uncertainty".** Rewritten, because inserting "of" would give the wrong meaning. The intended meaning is that once the fluid discretisation error is a few tenths of a per cent of the peak, it is of the same size as three other contributions:
- the spatial error of a reference at a published tutorial resolution (0.53%);
- the difference between the two couplings (0.08–0.23%);
- the model-form difference (≈0.06%).

Closer agreement therefore needs all of them reduced.

**Solver-limitation sentence.** Removed (the serial-only high-order solid and the 4-rank linear-reconstruction run, with its pending marker). The scientific limitation stays: "The temporal error on the finer fluid meshes is not measured."

## Grid names (all cases)

**Convention.** Stated once, in the Section 5 introduction: mesh levels are L1, L2, … from the coarsest. There are two exceptions, each with its reason:
1. Uniform square meshes (manufactured solution, Mok cavity) are named by their cells per side, n. Their time step is tied to n (Δt = 0.32/n; the Mok Δt sequence), and the Mok figure's axis is 1/n.
2. The uniform block-structured Turek–Hron family keeps 1x–8x. Its meshes are used for FSI1, FSI2 and the FSI3 coupling studies, and this keeps it distinct from the Gmsh family L1–L3 of FSI3 in the same subsection.

Levels of other codes keep their authors' numbering (Featflow level 4, Failer & Richter level 6, the FEniCSx S1–S3). The exception is where those codes use our meshes: the Kratos 3-D tube and beam runs are named after the meshes they split.

The Hessenthaler rF0–rF2 had no strong reason to differ: the companion paper calls the same meshes levels 1–3. They are now L1–L3, which also matches the companion paper.

**Changed:**
- Ring: table `tab:ring` mesh column 1/2/4 → L1/L2/L3. Its first column is now headed "Added mass", so that "Level" no longer means two things. The Numerical model and the Section 6 tolerance table ("mesh 2" → L2) and cost paragraph ("mesh-2 discretisation") are updated.
- Womersley tube: Numerical model and table caption ("mesh factor 4/2" → L3/L2).
- Collapsible channel: F1–F4 / "fluid level 3" → L1–L4 / L3 in the text, the Fig. 6 caption and the regenerated figure. Also in Section 6 (tolerance paragraph and table) and in the FEniCSx supplement (`supplement_fenicsx_collapsiblechannel.tex`, its solids4foam F1–F4 → L1–L4).
- 3-D tube: table `tab:3dtube` header (Level 1–3 → L1–L3) and all the text, including the Kratos paragraph and the to-do. The figure legend is regenerated ("Level k" → "Lk"). Also changed: the provenance rows, the cost table, the defects paragraph of Section 6, and Appendix D (interpolation degeneracy).
- Flexible-bottom cavity: Results ("levels L1–L4"; "on level~1/2/3" → L1/L2/L3), the Section 6 orders table, and Appendix D (pseudo-time-step paragraph). The cost table's "level 4" → L4.
- Beam in cross-flow: levels 0–3 → L1–L4. This shifts the numbering by one, so "the fifth mesh (level 4)" becomes L5. It is applied in Appendix A (text, table `tab:beam_seq` rows and orders, the Kratos paragraph, costs), the Section 5.4 summary, the Section 6 tables, the Recommended-suite to-do and the Conclusions.
- Hessenthaler: rF0–rF2 → L1–L3 in the text, both tables, the Fig. caption, the Section 6 cost table and the Discussion. The paper-data paths (`hessenthaler_rF/`) are unchanged.
- Blob in treacle: the four mesh levels are named L1–L4, and the baseline is L2.
- Turek–Hron provenance rows: "block-structured family levels 1–3 / level 4" → "1x–4x / 8x". The provenance audit confirms that its level 4 is 8x.
- CRediT entry: "level 3 of the 3-D tube" → "L3".

**Unchanged on purpose:** the manufactured solution (n), Mok (n²), Turek–Hron 1x–8x, Featflow/Failer levels, the FEniCSx cavity and blob supplements (their own mesh names), and the paper-data file names (`F1.dat`, `robin_backward_mesh1.csv`, `rF1_x0.vtk`).

**Figures regenerated:**
- `figures/verification/collapsibleChannel_mesh_study.pdf`, with `paper-data/figure_scripts/collapsible_mesh_figure.py`: legend and annotations L1–L4; panel (d) "max difference (% of peak)". It prints the same 4.813/1.592/0.680/0.334% as before.
- `figures/3dtube_pointA_histories.pdf`, with `paper-data/3dtube_level3_corrected_ecore/figure/tube_histories_figure.py`: legend L1–L3. The command in its README and docstring is updated.

The Hessenthaler geometry figure contains no mesh name, so only its caption changed.

## Reference values (general request)

**New Appendix E, "Reference values"** (`appendixReferenceValues.tex`, before the provenance appendix, so that Appendices A–D keep their letters). It has compact tables, each with the definitions needed: point, window, statistic and time.

| Case | Table | Content |
|---|---|---|
| Ring | `tab:refvals_ring` | exact ω and wet/dry ratios (the ratios follow from the exact ω); our L3 values `[pending]` (IQN-ILS rerun) |
| Womersley tube | `tab:refvals_womersley` | exact k and wave speed; L3 errors from `tab:womersley` |
| Manufactured solution | `tab:refvals_mms` | L2 norms at t = 0.4, n = 128; quasi-monolithic from `tab:mms_fsi_corr`; partitioned `[pending]` |
| Collapsible channel | `tab:refvals_collapsible` | first trough, second peak, second trough and third peak (values and times) of extrapolated oomph-lib and L4, their differences in % of P = 0.216 m, and the 0.33% maximum history difference |
| Turek–Hron FSI1 | `tab:refvals_fsi1` | 8x values, conditional Richardson values and GCIs, Featflow |
| Turek–Hron FSI3 | `tab:refvals_fsi3` | L3 means and amplitudes, provisional GCIs, Failer & Richter |
| Turek–Hron FSI2 | — | stated as `[pending]` (Gmsh rerun) |
| 3-D tube | `tab:refvals_3dtube` | point-A QoIs on L3 and the L2→L3 change (C-weak, no extrapolation) |
| Mok cavity | `tab:refvals_mok` | Kratos reference (peak band; trough and mean with GCI) and our 96² values |
| Cavity, beam, Hessenthaler, blob | `tab:refvals_provisional` | finest-level values, extrapolated values with GCI or conditional label where allowed, category; the SU2–CalculiX cavity limits, the Richter values and the measured flap deflection in the text below |

The collapsible-channel extrema are the only values not printed in the paper before. They are read from the stored histories (`paper-data/collapsible_mesh_study/F4.dat`, `oomph_extrapolated_dt160.csv`) with the tutorial driver's own method: the vertex of the parabola through the three samples nearest each extremum. This reproduces the driver's L4 first trough (−0.2161029 m) to all digits and the paper's −0.21610/−0.21615 m.

Each case's Results paragraph now has one sentence pointing to its table. For the beam, the sentence is in the Section 5.4 summary. The Section 5 introduction also points to the appendix.

**Machine-readable copy.** `paper-data/reference_values/reference_values.csv` (68 rows), with columns for:
- the definition and units;
- the finest level and value;
- the extrapolated value and the uncertainty with its type;
- the reference, its source and its category;
- the status, and the paper table or paper-data source.

`README.md` there explains the columns, the collapsible extraction and the pending items.

## Pushback / for your decision

- I kept "difference" rather than "error" for the collapsible channel throughout, as argued above. If you prefer to treat extrapolated oomph-lib as the reference and say "error", the caption, the axis label and four sentences would need to change back.
- The ring table edit (mesh column and "Added mass" header) touches rows that the IQN-ILS rerun will probably also change. A merge conflict there is trivial: keep the new values with the L1–L3 names.
- The beam renumbering (0–3 → 1–4) means "the fifth mesh" is now L5. The paper-data CSVs of the beam still use 0–3.
