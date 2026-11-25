# Response to Philip's notes on the abstract and Sections 1–4

Branch `text/review-sections-1-4`. Each item below quotes or paraphrases one comment in `philip_notes.md`, then says what was done, where, and any pushback. Section 5 onwards was changed only by the global Turek–Hron rename, the "Because" rephrasings and pointers to the new appendix.

## Abstract

**Too long; check the journal limit.** The Journal of Fluids and Structures Guide for Authors limits the abstract to 250 words ("concise and factual … must not exceed 250 words"; <https://www.sciencedirect.com/journal/journal-of-fluids-and-structures/publish/guide-for-authors>, seen through a web search because direct fetches return HTTP 403; a third-party listing, scraiber.com, gives the same limit). The abstract (`article.tex`, `\abstract{…}`) has been rewritten from about 590 words to 248, not counting the `\hl` marker.

**Emphasises our solver's problems; headline the benchmark solutions.** The new abstract headlines:

- the problem the paper addresses;
- the reference-quality classification and the separate assessment of spatial, temporal and iterative errors;
- the cases, including the new manufactured solution;
- findings for people verifying FSI solvers:
  - a nominally second-order code does not reach second order for every quantity, as the manufactured solution and FSI3 show;
  - why the published references disagree: time integration in the 3-D tube, the printed inflow of the beam benchmark, and the FSI3 summary values against the detailed tables;
  - the coupling convergence test.

The solver-history sentences ("traced to path-dependent mesh motion…", "verification exposed two solver defects", the IQN-ILS and quasi-monolithic details) are gone. The FSI3 sentence keeps its `\hl{[FSI3: ring+box, to be replaced]}` marker.

**Bold O(h) and 7.5 × 10^5.** The cause is in `sn-jnl.cls`. `\Authorfont` ends in `\boldmath`, and `\@maketitle` expands `\artauthors` outside a group, so the bold maths version leaked into the affiliations and the abstract. It is fixed at the cause with an `etoolbox` `\patchcmd` in the `article.tex` preamble, which groups the author list; the class file is unchanged. I checked this on a minimal document: the abstract maths font changes from CMMIB9 (bold) to CMMI9. The new abstract happens to contain no maths.

**"Hron–Turek FSI3".** It is now "Turek–Hron" throughout the built text and tables, including Sections 5–8, the appendices and the FEniCSx blob supplement. Labels, file names and the case identifier `TurekHron2006` are unchanged. `sectionTestCases.tex` is not in the build and was not touched.

## Introduction

**Nested brackets.** The paragraph is rewritten (`sectionIntro.tex`, second paragraph) as two sentences without nested brackets. ALE is still defined there.

**Wider applicability (immersed boundary, meshless and particle methods, weakly compressible explicit solvers, beam and plate solids).** Agreed. The paragraph now separates two things: the formulation used for the calculations (strongly coupled, partitioned, body-fitted ALE fluid with a Lagrangian solid) and the cases and references, which are statements about continuum problems. It says the cases are intended for:

- immersed-boundary and fully Eulerian methods;
- meshless and particle methods;
- weakly compressible explicit solvers;
- codes with beam or plate solid models.

It also says that the model-form differences stated with each case show where a reduced model can be compared directly and where it cannot. The text makes no claim that we have tested those methods.

**Keep "verification exposes real defects…" in the contributions.** Kept unchanged, apart from the rename.

**Acknowledge the AI role in the manufactured solution.** Done in two places:

- The "Declaration of generative AI" in `article.tex` now says that the derivation, code generation and checks were developed with substantial assistance from Claude (Anthropic), under the authors' direction, and verified by the authors' symbolic and numerical checks. It is followed by `\hl{[author: insert model name and version]}`. No model name is written; please fill it in.
- A one-sentence pointer has been added at the end of the "Construction" paragraph of Section 5.1.3 (`sec:mms_fsi`).

## Section 2

**"Following OpenFOAM".** Now "Following OpenFOAM notation". The repeated definition of the ALE abbreviation in Section 2 has also been removed, because it is defined in the Introduction.

## Section 3

**"explicit non-orthogonal correction".** Now "deferred non-orthogonal correction" (Section 3.1). It is the only use in the built paper. One more occurs in the unbuilt `sectionNumericalModelQuasiMonolithic.tex`, which was left alone.

**Check the schemes in 3.1 (limited convection?).** I checked them against the fluid `fvSchemes` of the solids4foam tutorials at the commits in the provenance table, and against the per-case text in Section 5. You were right that the old wording ("central or limited linear-upwind") was wrong for most cases. Section 3.1 now states:

| Convection scheme | Cases |
|---|---|
| Central (`linear`) | ring, Womersley tube, manufactured solution, collapsible channel, corrected FSI3 |
| Linear-upwind with an unlimited least-squares gradient | FSI1 (four levels, `a421f063d`), 3-D tube, flexible-bottom cavity, blob |
| `LUST` | Hessenthaler flap |
| Gradient-limited linear-upwind (`cellLimited leastSquares 1`) | only FSI2 (`18a052637`), baseline FSI3 (`d49552f4c`), Mok cavity (`6e69ddb6c`) and beam (`aee0c8e35`) |

Section 3.1 also now mentions:

- the skew-corrected continuity flux of the cavity;
- pseudo-time stepping for the steady problems;
- a pointer to the provenance appendix for the per-case set-ups.

Two cases need a check:

- The Mok coupling-test study at `a421f063d` uses the unlimited scheme, unlike the Mok results.
- I could not see the Hessenthaler PR #676 commit `64896b0a5` locally, so its scheme is taken from Section 5.4.3.

**#748 non-orthogonal correction to-do.** The `\todo` stays. It is reworded to say that an A/B comparison with the corrected boundary conditions is running and that the paragraph will report its outcome.

**"the highest available".** Removed.

**Sentence ending in "…are."** Rephrased: "…the fluid discretisation, the BDF2 time integration and the interface load are second order, so the coupled scheme remains formally second order."

**Sentences starting with "Because".** All seven in the built paper are rephrased:

- one in Section 3.5 (Robin–Neumann);
- six in Section 5: ring (two), manufactured solution, FSI3 diagnosis, 3-D tube and blob.

**Appendix of learning points ("this difference exposed a solver defect…").** Agreed and done (see Section 4 below). The sentence in Section 3.6 now says why the compiler and optimisation level are recorded and points to the appendix. The defect detail in Section 3.4 (interface transfer) is replaced by a pointer.

## Section 4

**Repeated abbreviation definitions; QoI in Table 1.**

- FSI is no longer redefined in Section 4. That definition also said "fluid–structure" rather than "fluid–solid".
- QoI is now defined at its first use, in the opening of Section 4, so Table 1's "QoIs" is defined before the table.
- Section 4.3 no longer redefines QoI.

**"It is used here only to organise…".** Now "Here it serves only to organise the evidence that each source has been controlled."

**Conditional and provisional indicators.** Rewritten in Section 4.2 ("Richardson extrapolation") so that the test comes first and the two labels follow:

- The test: whether four or more levels are monotone with a stable but sub-nominal observed order.
- If they are, the Richardson value and the GCI are *conditional* indicators, conditional on that order persisting.
- If they are not, a GCI from the observed order is a *provisional* indicator.
- Neither is a bound.

**"The campaigns in progress…".** Rewritten as a statement of which studies use the strict test and which use the legacy test. It is followed by `\hl{[update this list if further studies are rerun with the strict test before submission]}`.

**"An audit of the solver logs…".** Cut to two sentences: what the audit found, and what that means for the legacy studies. The per-case percentages and the steady-case details are now in Appendix D, "Coupling convergence tests".

**Provenance note (re-run the cheap cases and restart the expensive ones before submission).** Noted; nothing changed in the text. The provenance table and its to-do now sit in Appendix E, so a refresh only touches that appendix.

**Table 2 out of the main body.** The provenance table is now Appendix E, "Provenance of the results" (`appendixProvenance.tex`), together with its paragraph and to-do. Section 4.3 keeps a one-sentence reference to it, and Section 3.6 points there too.

## New appendix D, "Lessons for FSI solver developers" (`appendixLessons.tex`)

This appendix is written as run-in paragraphs rather than a table, because a 12-row, 4-column table would not fit on one page. Each paragraph gives the symptom, the evidence that revealed it, the cause and its status. The lessons are tool-agnostic, with solids4foam as the example, and the text has no PR or issue numbers. The paragraphs cover:

- aliasing in field algebra under `-O3`;
- scale-dependent barycentric weights;
- coupling convergence tests (the moved audit detail, and the round-off floor of a purely relative test);
- boundary non-orthogonal corrections that vanish silently;
- interpolation of the continuity flux without skew correction;
- a steady state that depends on the pseudo-time step;
- path-dependent mesh motion;
- wall pressure and wall velocity conditions;
- errors that only an exact finite-deformation solution exposes (manufactured solution, Womersley end flux);
- robustness of the coupling and the fluid iteration (IQN-ILS on the manufactured solution, the unrelaxed PIMPLE loop in Mok, the predictor);
- mesh families and geometry (size jumps, wall arcs);
- unstated benchmark definitions.

The cross-benchmark defects section (`sec:defects`) now points to the appendix. Not included are the `#codeStream` cache, the log-line units and the cfMesh octree items from `solids4foam_learnings.md`, because they are tool-specific.

**Pushback or open points:**

- The non-orthogonal-correction paragraph of the appendix contains `\hl{[insert the outcome of the A/B comparison]}`.
- Section 5 still tells the FSI3 and 3-D tube stories partly in chronological form, for example "baseline" and "corrected" set-ups, and an earlier wall condition. I did not change that, because you have not reviewed Section 5 yet. When you do, much of the FSI3 diagnosis could move to Appendix D.

## Build

`latexmk -pdf article.tex` finishes with no errors and no undefined references. The remaining warnings are the duplicate hyperref destinations that were already there before these changes.
