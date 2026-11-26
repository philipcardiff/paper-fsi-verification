Feel free to push back on any of the comments below.

section 5 Verification benchmark suite

it feels like the provisional cases whose evidence is incomplete (Section 5.4) should be removed from the paper, given they are incomplete and the paper is already quite long; or those cases should be made complete. We should consider each case and see if it can be upgraded easily or if not then removed. For example, I thought the flexible-bottom-cavity already has good refereneces?
also, looking at the other cases, they all look complete or can be made complete before submission - please check and come back to me to discuss as needed.
When we remove this section, we can remove text like "For each mature case" since we will no longer need to distinguish between mature and incomplete cases.

Fig 1 looks great.
one minor point - not directly aimed at this figure - is that the 3d images are all shown as 2d projection - this is ok but it might be good to have some 3d perspective-style images somewhere in the paper or otherwise a reading who quickly skims might think the cases are only 2d.
The Fig 1 caption is quite comprehensive but it feels very long and tough to read - I wonder if we can find a simpler solution while still being clear - if not, it is fine.

Table 2: in the mesh column, some rows say "(fluid)" - I dont see where this is defined or what it means (maybe i an missing it)

Fig 2 caption:
- "In (c) and (d) the transverse direction is stretched three times" is it stretched or compressed? please check
- "The other benchmarks’ geometries are shown with their sections." I dont understand this English

5.1.1 Ring added-mass problem

"The exact linear eigenvalue ..." this sentence seems incomplete (it looks like "is" is missing from the end of it) - it seems like a stub of a sentence - please fix this and check for similar stub sentences in other case "reference" section.

"(dry to strong)" is a bit unclear -> do you mean, dry and wet (weak, moderate, strong)?

" the loss is unchanged when the strict-test coupling tolerance is reduced from 10−4 to
10−8, so it is not an iterative error (Section 6.3); the mechanism has not been identified, and
the frequency is unaffected at the reported precision."
interesting - ask a sub-agent to investigate this; at a minimum, we should be able to present a plausible explanation for this.

"Independent FEBio comparison."
I wonder if this section should be called "Independent numerical comparison" or similar and then we can use the same subsection title across multiple test cases. This also emphasises the fact we have an independent numerical comparison rather than emphasising the code name (which we mention anyway in the text).
In fact, the subsequent "Independent partitioned comparison." feels like it should be combined with this section into one combined "Independent numerical comparisons" or something like that. In that way, it becomes clear that these "other" code solutions are just to provide confidence (and maybe sme insight) with (into) the solids4foam solution.

"changing K by a decade" I suggest saying "by an order of magnitude", which would be more natural to me

"In the moderate regime, changing K by a decade in
either direction (107 and 109 Pa) changed the wet/dry ratio by−0.0011% and +0.0063% of the
reference ratio, increasing the regulariser stiffness tenfold by +0.0015%, and halving the amplitude
by−0.0068%;"
it is not clear to me here if "−0.0011% and +0.0063%" refers only too the K change or also to regulariser stiffness and amplitude. this needs to be makde clearer.

"and the incompressible limit is bracketed, not reached"
i dont fully understand what you mean by bracketed here -> be more explicit.

"Frequency and damping": to me, this section seems more like a description of the FEBio results than a comparison with our (solids4foam) results (as the section heading suggests).
With the purpose of the paper in mind, let us consider the purpose of including independent numerical references: to me, it is to provide confidence that our proposed solution (from solids4foam) are in fact correct. The purpose is not to provide detailed analysis of every consideration with the independent reference code to the same level as our solidsfoam solution. Possibly, on a case by case basis, the independent solution may help us provide insight or check our hypothesis or maybe just highlight differences between our FVM and a FEM approach in terms of its characteristics or solution characteristics.
Having said that, it is certainly useful to have these details for FEBio but i wonder would they be better moved to an appendix - for now - and possibly even removed from the paper later, if we decide that (then they could live in supplementary material) later.
Dont blidnly move everything to the appendix but consider what information actually supports the points we want to make and the purpose of the paper.
Related to this, I am struggling to understand Fig 3(a) - but maybe it is not needed anyway.

"and a solids4foam verification of these quantities would need a mesh that resolves the layers."
is this an important point to make? or just a note. At the moment it seems to be made quite strongly whereas it feels like it is just  caveat at most or a bound on the validity of the approach for others who will reproduce the case.

"Supplementary material" - I feel like these sections should be removed fromt the main paper; instead maybe an appendix can describe the contents of the supplementary material

"Status" section - i dont understand the purpose of this section or why it is here - what does it add? can it be removed?

"Independent partitioned comparison." see my comment above about merging all other numerical references into one section and tightening the purpose of these sections. In fact, this same approach should be adopted for all test cases.

"Lesson": I guess this should be "Lessons" - this is generally good, although for "and the solids4foam runs do not resolve the viscous layers;" i feel like we need to briefly note the consequence (or lack thereof) of this point as an uniformed reader may take this as underminng the validity of the results or something like that.
Another learning of this case for me is that the partitioned FSI solvers work just fine on this case, which was initially a surprise for me since strong added mass cases are often noted as being difficult for partitioned approaches. However, here, it seems that the combined added mass and solid deformation do not pose any challenges. In contrast, the "collapsible channel" clearly poses challenges for the partitioned solvers (whereas the monolithic solver has no issue) so it is as if the ring case still has an effective low added mass (maybe the wrong phrase) while "collapsible channel" has an effective high added mass. Or maybe our "Strong" added mass case is not actually that strong for this type of flow/deformatio - I am not sure. But i feel like a comment on this is warranted here.


5.1.2 Pulsatile compliant tube

"usually in a thin-wall, long-wave form whose modelling error is not negligible for
a thick wall. "
I am confused here by what you are saying - are you saying that you/we are using a thin wall or not? If not, make it clear wht we use (it is fine to note that most people use thin theory but we are using thick theory - if that is what you are saying)
This "What it tests" section needs to be more clear - do the same for all cases.

"an untethered axisymmetric tube" what does untethered mean here? this should be clear - it sounds like traction free on the outside.

"Reference. The exact linear solution for ..." Once again the English seems incomplete here. Use full sentences and correct all other similar reference sections.

"used as a reference, it would have hidden every discretisation error studied here."
ah ok, now I understand above that you are saying that many people use thin theory but that this is not good if you are interested in examining discretisation (and other) errors -> this is an important point from this case.

"First Fourier coefficients over periods" we should be clear what these coefficients are and how we calculate them.

"amplified about sevenfold, since Im(k) is only
14% of Re(k)."
the implication of Im(k) being 14% of Re(k) is not clear to me here - spell it out more - make it easy for the reader.

"Time." This is an important learning but given the details involved i wonder whether it would be better to briefly note it and refer to the appendix for details?

"Limitations: the wall is untethered, unlike Womersley’s constrained tube, the high-order
solid cannot be used on wedge meshes, "
How are these limitations from the paper/science perspective?
Womersley’s constrained tube is a different case which may or may not show something else - i dont see how this is a limitation.
"the high-order solid cannot be used on wedge meshes" - this is a solids4foam limitation - not a limitation of the analysis shown here - do not give such software limitations in sections like this unless there is a clear reason from the paper's perspective.

5.1.3 Manufactured solution for transient ALE FSI with a finite-strain solid

"The case is not yet part of the released solids4foam suite; until it is, the case libraries
and the case generator are in paper-data/mms fsi/case/."
Create a PR against solids4foam for this case if one does not exist

"Role." This may be an appropriate sub-section head but I think we should - for the most part - standardise the subsection headings across cases, e.g. this is the third case but it is the first one to use "Role". Standardised headings will help readers compare cases and understand how they differ and how to contribute to a suite of cases.

"Exact solutions of physical FSI problems at finite deformation are rare; those of Hessenthaler et al. [49]
include neo-Hookean solids and transient flows in channels and tubes, but the motion ..."
As it is currently phrased, it implies that Hessenthaler is an exact solution, which it is clearly not. Unless I am missing something, it should not be mentioned here. In contrast, it makes sense to refer to other similar MMS type solutions here.

"AUTHOR NOTE: Novelty statement revised from Fisch et al"
please address/resolve this author note so it can be closed - I have added the full PDFs for references 51-54 to "/Users/philipc/Documents/Papers/Literature - sorted" on my mac where the files are named as "first author surname - year - paper title"

"The solid is the compressible neo-Hookean law of the code"
"of the code" does not make sense here -> i guess it is referring to solids4foam's implementation here but instead it should just say neo-Hookean as per our definition of neo-Hookean earlier in the paper.

"clamped on its top and on its ends"
"on its ends" could be ambigous - surely it is better to say left and right or give coordinates.

"in the flat variant "
this is the first time that variants have been mentioned - up to now, the reader would presumably think there is just one solution/case/variant - we need to explain the two variants before this, maybe at the end of the role section.

The manufactured solutions really are stunning! One minor comment: Appendix C should clearly give the source terms (body forces) that are added to the various equations so that readers can easily reproduce it in their own codes.

"and its trace fixes the interface motion"
I know this may be common mathematical language but I am unsure wha you mean by trace here - I guess a projection. Revisit this description and ensure a typical computational mechanics reader can easily follow it (as well as a mathematician)

The Construction section is good but given the paper's length i wonder if some of this should be moved to the appendix and summarised here.

"Numerical model. The partitioned calculations use the procedures of Section 3. The fluid
(pimpleFluid) uses BDF2, central convection, Gauss linear gradients, one outer corrector with four
pressure correctors and one "
As noted in my comment above about having consistent subsections, it seems strange that these settings were not mentioned previously - presumably most cases use the same or similar settings and we could state it once and then only note the changes. I am not sure - reveiw this.

"The fluid does not: its velocity converges at orders 1.2–1.4 and its pressure at 0.9–1.2 on every pair of levels, " and "The
cause of the outlet layer and of the pressure-level offset is open"
We should fix this since it will overshadow the readers' confidence in our other results. Start a subagent on this to find and fix the issue. Since the QM solver works, we can surely find and fix the issue.
This will resolve the "AUTHOR NOTE: Open: cause of the first-order pressure layer ..."

"IQN-ILS. IQN-ILS did not complete this case"
For the main text, we can just briefly note that IQN-ILS struggled and hence is not shown - maybe with a very brief one line explanation. I think the text in this paragraph probably is snot needed in the paper - you can put them in the learnings appendix in sme frm.

"Partitioned solid." and "Partitioned fluid."
I know we also briefly noted the quasi-monolithic solver (which may not be needed if we resolve the order issue), but it seems strange to have sections with these names.

"At t = 0.2
the partitioned solid and interface displacements converge at first order (0.85–1.08) on every pair
of levels, with almost the same errors for the two couplings, whereas the quasi-monolithic solid
displacement converges at 2.02–2.16 "
THis seems like it should also be resolved like the fluid issue above - maybe the same aagent.
Ah wait, "one of the corrections of Table 8", maybe it was already fixed (maybe the fluid was too?); if so, we should just state the final results with the corrections. The lessons/learnings along the way can still be referred to in the dedicated appendix and mayabe briefly noted in the main text.

"AUTHOR NOTE: Pending: partitioned A/B on n= 64 and 128 "
Check if this has been done - if not, do it.

"on by default in the segregated solver"
A reader does not care about defaults in solids4foam - do not mention these things unless there is a strong reason. remember who the readers are and the focus on the paper.

"Quasi-monolithic coupling and boundary non-orthogonality" and "Quasi-monolithic coupling with three corrections."
These details are not relevant to this case/paper as far as I can see.
They can go in the learnings appendix.

"Non-orthogonal correction of the moving-wall velocity."
I think this can go in the learnings appendix. It is actually very useful for openfoam readers to know about this boundary non-orthogonal correction limitation (even for pure CFD).

"; an incompressible solid would need a mixed formulation. " I am not sure this is a limitation per se; sure I understand thaat 2D is, but incompressible solid seems more like a different feature to explore than a limitation - and especially i would not lead with this.

Comments up until before "5.2 Independently converged numerical reference"