# Response to Philip's notes on Section 5, up to Section 5.2

Branch `text/section5-review`. Each entry paraphrases one comment in `philip_notes_sec5onwards.md`, then says what was done and where, or why not. Items marked **in progress** are being handled separately and are flagged in the text with `\hl{[pending: …]}`. The PDF grows from 92 to 95 pages: the main text (up to the Conclusions) is two pages shorter, and the appendices are longer because material was moved there and a Supplementary Material appendix was added.

## Section 5, general

**Remove or complete the provisional cases (Section 5.4); drop "For each mature case".** In progress: a separate assessment is under way, and you will decide. The provisional cases and the "mature case" wording are unchanged, apart from their headings (below).

**Fig. 1 looks great; add some 3-D perspective images somewhere.** Not done here, because it needs new figure work. I suggest a perspective view of the 3-D tube or the Hessenthaler flap, either in Fig. 1 or with its case.

**Fig. 1 caption too long.** Shortened to one line per panel (`sectionBenchmarkSuite.tex`, `fig:solution_gallery`). The meshes, couplings, times, clipping and drawing scales are now in a new "Panel details" section of `paper-data/solution_gallery/README.md`, which the caption cites.

**Table 2: "(fluid)" undefined.** The caption now says "``(fluid)'': only the fluid mesh is refined". This is the case for the collapsible channel, where the solid is fixed at 160×8, and for the Hessenthaler flap. The "case not yet in the released suite" footnote is removed.

**Standardised headings (from the FEBio and "Role" comments).** Every case now uses the same headings, in this order: What it tests, Problem, Reference, Numerical model, Results, Independent numerical comparisons, Lessons. The Section 5 introduction lists them, and says that each case's numerical model gives only its differences from Section 3.
- In 5.1, all seven headings are used. The tube and the manufactured solution say "None; the reference is exact" under Independent numerical comparisons.
- In 5.2 onwards (collapsible channel, Turek–Hron, 3-D tube, Mok, and the provisional cases), only headings were changed and content moved, as instructed. "Evidence" became Results (or Numerical model plus Results for Turek–Hron). "Lesson" became Lessons, and "Role" became What it tests. "Reference and model-form difference" and "Reference: an independent solution…" became Reference.
- Other-code items were moved verbatim into Independent numerical comparisons: FEniCSx for the collapsible channel, preCICE/SU2–CalculiX for FSI1, Kratos for the 3-D tube, the cavity's references, and FEniCSx for the blob.
- The Hessenthaler sub-headings (interface pressure condition, convection scheme, fluid meshes, coupling iterations, validation) were folded into Numerical model and Results. The validation paragraph now opens with "Separately from the verification, the solutions are compared with the measurements."
- Where a later case has no paragraph for a heading, no content was invented.

## Fig. 2 caption

**Stretched or compressed?** `paper-data/figure_scripts/geometry_panel.py` multiplies the transverse coordinate by 3 (`s = 3.0`, `sy = 3.0`). The transverse direction is therefore stretched, so the original was correct. The wording is now "drawn at three times the scale of the axial direction".

**"The other benchmarks' geometries are shown with their sections."** Now reads "Geometries of the other benchmarks are shown in their own sections."

## 5.1.1 Ring added-mass problem

**Stub sentence in Reference.** Rewritten as full sentences. I checked the Reference paragraphs of the other 5.1 cases: the tube's is also rewritten, and the manufactured solution has none of these stubs.

**"(dry to strong)".** Now "for the dry ring and the weak, moderate and strong levels". The Problem paragraph introduces "the dry ring … compared with three wet levels, weak, moderate and strong".

**Amplitude decay: "the mechanism has not been identified".** In progress (another agent). It is replaced by `\hl{[pending: mechanism]}` in the Coupling item.

**"Independent numerical comparisons", merging FEBio and the partitioned comparison.** Done. One paragraph covers FEBio and SU2–CalculiX–preCICE. It states what they establish:
- Confidence: three different discretisations reproduce the exact ratio to 0.1% or better, and SU2–CalculiX differs from solids4foam by at most its own discretisation error (0.05%).
- Insight: FEBio verifies the viscous shift and the decay, which the solids4foam meshes do not resolve. It also shows an FE/FV difference: nodal no-slip with an unresolved layer gives −6.2%, whereas the cell-centred FV solution approaches the inviscid response on the meshes used.

**"by a decade" → "by an order of magnitude"; attribution of the K/regulariser/amplitude changes.** Now split into separate sentences, one per change, with each value attributed. Moved to Appendix B.

**"the incompressible limit is bracketed, not reached".** Made explicit in Appendix B. Over the tested range of K, compressibility changes the ratio by less than 0.01%. The K runs do not converge towards K → ∞, because the change from 10^8 to 10^9 Pa (0.0063%) is larger than the change from 10^7 to 10^8 Pa (0.0011%). So the limit can be neither reached nor extrapolated. The main text states the K, regulariser and amplitude changes as moderate-level sensitivities, not as error bounds.

**FEBio detail: keep only what supports the paper.** The main text keeps the FEBio set-up in one sentence, its headline results, and what they establish. Everything else is in a new Appendix B, "Independent computations of the ring added-mass problem":
- the model and the model-form sensitivities;
- the layer-resolution requirement;
- the time-step and resolution checks;
- the impulsive start;
- the comparison with the solids4foam damping;
- the free-slip runs;
- the SU2–CalculiX detail (meshes, start, dry-ring error, time order, amplitude loss and artificial compressibility).

The existing optional FEBio supplement (graded-mesh ladder and interface mechanics, behind the `\ifsupplementseparate` toggle) is now a sub-part of that appendix. Note: if the toggle is switched to a separate supplement, `supplement_febio_ring.tex` should be checked against the new appendix.

**Fig. 3(a) hard to understand.** Moved with the FEBio table to Appendix B, because it no longer carries a main-text point. The caption is rewritten:
- solid bars show the error against the inviscid ratio, and hatched bars the error against the finite-viscosity ratio;
- the hatched bars are near zero, so FEBio reproduces the viscous problem;
- the solid bars equal the viscous shift;
- panel (b) falls towards the inviscid ratio.

**"would need a mesh that resolves the layers".** Toned down to a caveat with its consequence, in a new "Viscous layers" item in Results and in Lessons.

Pushback on the premise that "the viscous shift is far below the discretisation error": this holds for the frequency at practical time steps, but not for the ratio on the finest mesh. The finest-mesh high-order ratio errors (+0.001 to +0.006%) are smaller than the viscous shift (−0.026% moderate, −0.08% strong). The data show something cleaner: with the layer unresolved, the solids4foam ratios approach the inviscid ratio on the meshes used. The text therefore says:
- the frequency and ratio verification is against the inviscid reference and is unaffected;
- the consequence is that the viscous shift and the viscous decay rate are verified only by FEBio.

**"Supplementary material" subsections.** The main-text paragraph is removed. A new appendix, "Contents of the Supplementary Material" (`app:supplementary`, last appendix), describes the documents: the FEniCSx collapsible channel, the flexible-bottom cavity, the two blob documents, and the optional separate FEBio document.

**"Status" paragraph.** Removed. Its essentials are elsewhere: the quantified model-form differences are in Independent numerical comparisons and Appendix B, and the free-slip runs are in Appendix B.

**"Lesson" → "Lessons"; consequence of unresolved layers; partitioned couplings easy here but hard in the collapsible channel.** Done, from the paper's own evidence and hedged:
- In the ring, both couplings converge at all levels, including m* ≈ 10 (IQN-ILS 4.9 and RN 12.3 iterations per step).
- In the collapsible channel, IQN-ILS fails with the strict test and for part of the wall-density sweep, and RN completes only for ρS/ρF ≥ 1.
- Explanation: in the ring the added mass acts on a single smooth, low-frequency mode of a linear problem, so the interface operator is fixed and effectively low-rank, which a few IQN-ILS secant updates or the Robin impedance estimate capture. The collapsible channel has a thin, very compliant wall, large flow-induced and history-dependent deformation, and many coupled modes.
- The text says that a nominal mass ratio is not a sufficient measure of difficulty, and that the explanation is consistent with the evidence but untested.

## 5.1.2 Pulsatile compliant tube

**"What it tests" unclear about thin and thick wall; the key lesson.** "What it tests" now states that the reference is the exact thick-wall solution with no thin-wall or long-wave approximation. It says that the thin-wall wave speed is 4.7% higher for this wall, so a thin-wall reference would have hidden every discretisation error studied. Lessons repeats this as the main lesson.

**"untethered".** Now "The outer surface of the wall is traction free; the tube is not tethered, i.e. its wall is not restrained axially as in Womersley's tethered-tube model."

**Reference: full sentences.** Done.

**"First Fourier coefficients over periods".** Defined in Numerical model:
- the formula: a = (2/N) Σ f(t_j) exp(−iωt_j) over the N time steps of periods 5–6, with f ≈ Re[a exp(iωt)], amplitude |a| and phase arg a;
- which histories are transformed;
- that the wave speed and Im(k) come from a least-squares fit of log p̂(x) along the tube at r = R/2.

Source: `womersleyTube/verification/README.md` and `scripts/womersley_tube_verification.py` (`fourier()`).

**Im(k) 14% of Re(k).** Spelt out: an error that is a fraction ε of |k| is a fraction of about ε/0.14 ≈ 7ε of Im(k), the small imaginary part.

**"Time" detail.** The Results bullet is shortened to the result and the reason. The O(Δt h) derivation, the 3200/6400-step estimate and the orders with the standard end condition are moved into the "Flux consistency at truncated boundaries" paragraph of Appendix D.

**Limitations: untethered and high-order solid on wedges.** Both removed: one is a software limitation, and the other is a different problem rather than a limitation. The remaining scientific limit is reworded: the attenuation can be verified only to ≈0.2%, the level at which linear theory holds at these amplitudes.

## 5.1.3 Manufactured solution

**solids4foam PR for the case.** In progress (another agent). The AUTHOR NOTE is replaced by "The case is distributed with solids4foam \hl{[pending: release]}."

**"Role" heading.** Now "What it tests". The two variants are introduced at its end, before first use.

**Hessenthaler sentence.** Partial pushback. `Hessenthaler2020` (CMAME 362, 112841) is "A class of analytic solutions for verification and convergence analysis of linear and nonlinear FSI algorithms". It is a different work from the flap experiment (`Hessenthaler2017a,b`) of Section 5.4.3, and its solutions are analytic. I could not find its PDF locally, so the shear-only characterisation is kept as previously stated, not re-checked. A reviewer would expect the closest exact class to be cited, so I kept one sentence that names it as distinct from the flap benchmark and gives the reason it does not suffice: shear motions only, no normal interface motion, no domain deformation. It is no longer presented as "exact solutions of physical FSI problems", and the paragraph then turns to the MMS literature. If you still prefer to drop it, delete that one sentence; it is also cited in the tube section.

**AUTHOR NOTE on novelty (refs 50–54; the numbering is 50 Étienne, 51 Tremblay, 52 La Spina, 53 Yu, 54 Fisch).** Resolved and removed. I read the PDFs in "Literature - sorted".
- Étienne et al. 2012 (confirmed): a St Venant–Kirchhoff large-displacement total Lagrangian solid; monolithic FE; steady 2-D and axisymmetric, and transient 2-D with ALE and implicit Euler. Traction continuity is satisfied without an interface source by solving for the fluid pressure and viscosity along the interface (μf(x,t), pf(x,t)). The solid's horizontal displacement vanishes on the interface (a(x,1,t) = 0), so the interface points move only vertically. The text now says this, and adds that our interface also slides tangentially.
- La Spina et al. 2020 uses Étienne's procedure ([52] there) to verify a hybridisable DG (HDG) fluid coupled to a CG solid. The old text said "finite element codes", which is now corrected.
- Yu et al. 2015: spring-mounted rigid bodies, 2-D and 3-D (confirmed).
- Fisch et al.: a partitioned OpenFOAM (cell-centred FV) fluid with a Carat++ pre-stressed geometrically nonlinear membrane, and MMS sources added to the interface traction balance (their Eq. 9b) (confirmed). The local file is named "2011", but the paper header is WCCM XI / ECCM V / ECFD VI, i.e. 2014, so the bib entry is correct.
- Tremblay et al. 2006 is not in the folder, so it is cited only as "see also"; the mechanism is attributed to Étienne et al.

No bib entries needed changes.

**"neo-Hookean law of the code".** Section 2 stated no law, so the compressible neo-Hookean law is added there as Eq. `eqn:neo_hookean`. It is the law used by all three neo-Hookean cases (the `neoHookeanElastic` in the tutorials). The case text and Appendix C refer to it, and Appendix C's duplicate equation is removed.

**"clamped on its top and on its ends".** Now "its top face Y = H + T is clamped, and its left and right faces, X = 0 and X = L, are clamped in the full variant and carry the exact displacement in the flat variant". The Fig. caption is fixed likewise.

**Appendix C source terms.** The Sources paragraph now gives every term explicitly:
- no mass source (∇·v = 0 identically) and no interface source;
- b_F = ∂v/∂t + v·∇v − ∇·σ*, with σ* = ν[∇v + ∇vᵀ] − pI and ∇·σ* = ν∇²v − ∇p;
- ρ_S0 b_S = ρ_S0 ∂²D/∂t² − ∇₀·Pᵀ, with P = JσF⁻ᵀ;
- the component forms, which fields and equations to use, and which variables are held fixed in each time derivative.

The closed forms run to thousands of terms after differentiation, so they are not printed. The text points to `derive/mms_fsi.py` and to the generated `derive/generated/mmsFsiExact.H` (and `generated_flat/`), whose `fluid()` and `solid()` functions return the fields and both sources at any point.

**"its trace fixes the interface motion".** Now "its value on the interface prescribes the interface motion".

**Construction too long.** Summarised in five sentences in the main text. The Möbius map, the stretched coordinate, the shear identity and the 11 checks are in Appendix C, which already held them.

**Numerical model repeats common settings.** Only the case-specific settings remain: gradients, correctors, solver tolerance, solid variant, mesh motion, couplings and tolerance, meshes and steps, and error norm. Section 3 already gives BDF2, Newton–Krylov, the direct map and the convection per case. The Section 5 introduction now says that each numerical model lists only its differences from Section 3. The diagnostic variants moved to Appendix D, and the boundary data are referred to Table C4.

**Fluid first order; the AUTHOR NOTEs on the open cause and on the pending partitioned A/B.** In progress (another agent; new final tables to come). As instructed, the current numbers are kept:
- Results now states only the results, behind one `\hl{[pending: final partitioned results after the fluid/solid fixes]}` marker.
- The diagnosis (start-up error, traction area, outlet, time-derivative flux correction, ruled-out causes, and the pressure-error figure) moved to Appendix D, "Errors that only an exact finite-deformation solution exposes", which has its own pending marker.
- Both AUTHOR NOTEs are removed.

**Partitioned solid at first order; is it already fixed?** For the quasi-monolithic coupling, yes. The Results now state its final behaviour with the three treatments (Table `mms_fsi_corr`). The earlier quasi-monolithic rows (full and flat) moved to a new Appendix D table (`tab:mms_fsi_qm`). For the partitioned couplings: in progress, as above.

**IQN-ILS paragraph.** One line in Numerical model: "IQN-ILS did not complete the case, because its least-squares update inverted the solid cells next to the interface (Appendix D)". The full paragraph is in Appendix D, "Robustness of the coupling…".

**"Partitioned solid"/"Partitioned fluid" headings.** Gone. The Results describe the solid, interface and fluid orders in prose, and the diagnosis is in Appendix D under "Solid" and "Fluid".

**"on by default in the segregated solver".** Removed, including from the moved text. "the default test" is now "the legacy test", the paper's own name for it in Section 3.

**QM non-orthogonality, QM with three corrections, moving-wall non-orthogonal correction.** All moved to Appendix D:
- the first two as "Quasi-monolithic coupling and boundary non-orthogonality", with the new table;
- the third merged into "Boundary non-orthogonal corrections that vanish silently", with the measured < 1% effect and a sentence that the limitation affects any OpenFOAM calculation with a varying wall velocity on non-orthogonal cells, including pure CFD.

The main-text Lessons keep one-line pointers.

**Limitations: don't lead with the incompressible solid.** "Limitations" is merged into Lessons. It leads with what the case exposed. The 2-D, compressible-solid point comes last, as possible extensions: a 3-D variant and an incompressible solid with a mixed formulation.

## Not changed

Sections 5.2 onwards: headings only, as instructed; no content changes. The Conclusions bullet on the manufactured solution still describes the current partitioned results and should be revised with the new tables. It also says the corrected quasi-monolithic coupling converges at second order "in space and time"; Section 5.1.3 now qualifies this, because the temporal orders of the fluid pressure and velocity fall to 1.6 and 0.9 on the smallest steps.
