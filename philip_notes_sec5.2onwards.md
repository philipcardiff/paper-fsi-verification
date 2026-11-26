Feel free to push back on any of the comments below.

5.2 Independently converged numerical reference

"based on the channel width" i wonder if channel "height" is a better description given its orientation

"and ReSt = 50;" is this right? is it meant to be St or the product of Re and St?

"The problem modifies that of Jensen and Heil (no pre-stress, h/a = 1/20,
clamped ends, pressure ramp, regularising density ρS = ρF, which shifts the response by 0.28% of
the peak) so that it is well posed for a continuum solid and partitioned coupling"
Reading this, it is unclear to me what Jensen and Heil did and how we modified it - i guess brackets are what they did but it is not clear how we modified it. this should be made clear.

"the published values are therefore not reproduced."
maybe it would be clearer as "the published values are therefore purposely not reproduced."

"whereas the code under test"
this is a bit unclear - i guess you mean the solids4foam FVM code we are using in the current work

"For the first trough, replaying the fluid load of
one coupled run open-loop on the continuum wall and on the beam gives at least about 0.06%,"
0.06% of what?
and what does open-loop mean?
make this clearer.

" (Supplementary Material, FEniCSx cross-code check of the collapsible channel)"
this is strange: FEniCSx is mentioned with no details -> if we are using fenicsx as further evidence of something then this should be clearer while it is still fine to leave the complete details of fenicsx elsewhere and refer to them.

"; the fluid and solid meshes are generated with blockMesh."
this clause it out of place in the "Evidence" subsection

Fig 6 caption: "maximum error against oomph-lib," is "error" the right word here or should it just be "difference"?

Fig 6: here the grids are called F1-F4 - we should be consistent in the grid names across test cases unless there is a strong reason to deviate. Please check all and update as needed.

"Second implementation. A monolithic finite element implementation by the same group "
"by the same group" could be confusing here to readers - i suggest we rephrase or clarify what we mean here - i.e. run by the authors of the current paper.

One purpose of the current work is to provide benchmark solutions so others can try reproduce them, but - at least for this case - I dont see easy to copy and compare with values, e.g. a table of values of peak/trough amplitude or similar. his may or may not be apropriate in the main body but at least it should be in an appendix. The same comment applies to all cases i.e. can a reader easily test their own code based on what we give. I know we provide all the code/cases too but surely the main results should also be in the paper.

"error falls to ≈0.5% the reference uncertainty,"
should this be "error falls to ≈0.5% of the reference uncertainty,"?

"and the high-order solid is serial only; a 4-rank
run with the linear-reconstruction solid gave an incorrect displacement [pending: re-check of the
parallel run]." do not comment on solver limitations (anyway i think this was solved or there is an open PR)
it is not relevant to the purpose of the paper and readers do not care!







