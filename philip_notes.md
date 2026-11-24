Feel free to push back on any of the comments below.

abstract

the abstract seems way too long. check what the abstract max length is for the intended journal.
it also feels like the current abstract emphasises the various challenges/issues of our chosen solids4foam solver approaches, which seem secondary to the goal of the paper.
I guess the goal is establishing these benchmark solutions, and as part of this, key findings might be that our nominally 2nd order FSI solver may not achieve monotonic convergence for some quantities in some cases, for example.
whereas comments like "This was traced to a path-dependent mesh motion and a first-order wall-pressure condition" are surely just issues in our code (which we subsequently fixed!) and are hardly of interest to the intended audience of this paper (people who want to verify FSI solvers).
Of course, we may like to comment on things like this in the main body of the work, e.g. highlighting how a code might appear to work but have issues, but it should hardly be the headline.

Also, in general, we should not give a chronological story of our own solver development for each case, e.g. we tried then X broke so we fixed it then Y broke, etc; instead, we just state the final aproach and final results and any limitations with that. Once again, think of the intended audience.

Similar comments apply to statements like "verification exposed two solver defects"

Having said all of that, if, for example, a key narrative of the paper is explaining differences seen in literature (E.g. where there is a large range in the predictions) then this is a key contribution and can be woven into the main narrative, including the abstract. But, it should not be focussed on issues in our own code which we subsequently fixed (but of course in the main body we can brielfy highlight learning points we had with our own solver as it might help others and motivate the need for true benchmark solutions). 

I dont see what some terms/symbols like O(h) and 7.5 × 10^5 are bold - they should not be bold as far as i can see.

"Hron–Turek FSI3" - I thought we already fixed this to be Turek-Hron?



Introduction

"(Eulerian, Lagrangian, or an arbitrary Lagrangian–Eulerian (ALE) fluid with a Lagrangian solid)"
It is a bit clunky to have nested brackets here

"The present work concerns strongly coupled FSI with a body-fitted ALE fluid and a Lagrangian solid,"
I agree that this is th4e approach we have adopted (and the approach of many of the references we refer to) but shouldn't we also state that the current references suitable more generally? e.g. if thinking of solids4foam, we have a new immersed boundary FSI approach and i intend to use the cases in this paper to benchmark it and assess how well it works. I imagine other FSI formulations would also find it useful, e.g. mesh-free or particle methods; also, I am a co-author on a weakly compressible fully-explicit FSI solver paper which uses cases like this. Or, in some of the cases here the solid could be represented in reduced form, like a beam or plate. Do you get my point?


"Evidence that verification exposes real defects and misleading benchmark agreement: two
solver defects found through a cross-platform discrepancy ..."
In contrast to my comments above, here i think it is appropriate to refer to issues we found in our own solver since it reports the main motivation for the paper; that is, without strong verification such errors can be missed.


"A closed-form manufactured solution ..."
If not already done, we need to acknowledge the role of Claude Opus 5.5 in the development of this solution, e.g. int he AI declaration section (although I wonder if it should also be stated in the case section itself).



2 Mathematical Models

"Following OpenFOAM, "
Should this be "Following OpenFOAM convention, " or "Following OpenFOAM notation, "?

Apart from that, section 2 reads very well.



3 Numerical methods

"with explicit non-orthogonal correction"
I think "explicit" can be confusing here as the correction terms are still implicit in time - i think "deferred correction" is less ambigous and is less likely to cause confusion for people unfamiliar with the openfoam approach.

For 3.1, we should check that this is correct since we have made several changes to the schemes used, e.g. I am not sure we use a limited convection scheme anymore.


"Pending result (A/B comparison; owner and date not yet set): on the OpenCFD
OpenFOAM releases used here, the non-orthogonal correction of the wall-normal ..."
We should fix this for this paper.


"(polynomial order 3, the highest available) "
I think we can remove "the highest available" here since it is not that relevant (I think) and we plan to remove this restriction in the code since most procedures are already fully general (the stabilisation may need to be generalised).

"the coupled scheme remains formally second order, because the fluid discretisation, the BDF2 time integration and the interface load (Section 3.4) are."
The English reads a bit strange here => it sounds strange ending the sentence with "are"

"Because it changes the iteration rather than only accelerating it, "
I am not a fan of starting sentences with because - rephrase all of them in the paper.

"and this difference exposed a solver defect (Section 5.3.2). Table 2 gives them per case."
I wonder if it would be a good idea to add an appendix to the paper with "learning points for FSI solver developers" or similar, where we catalogue the various issues we found, how we unearthed them and how we resolved them" - what do you think?


4 Verification methodology

Be careful that we are not repeatedly re-defining abbreviations (e.g. FSI).

"It is used here only to organise the evidence that each source has been controlled"
The English seems wrong here

Table 1: do we need to define "QoIs" abbreviation?

"Where a sequence of four or
more levels is monotone with a stable but sub-nominal observed order (Hron–Turek FSI1), the
Richardson value with that order and the GCI of the finest level are given as conditional indicators,
labelled as such; the GCI is centred on the finest-level value. Where a sequence does not satisfy it, a
GCI computed with the observed order is labelled a provisional indicator, and values computed with
an assumed order are labelled as such"
I found these two sentence so sound almost the same when reading them and I had to re-read them multiple times to understand the distinction. I wonder if they can be better phrased.

"The campaigns in progress ..."
we need to update this since the paper cannot refer to ongoing simulations
for text like this, it should be highlighted or flagged to avoid being missed

"An audit of the solver logs ..."
I am not sure much of the text after this is suitable for the paper - yes, there is a learning point but i dont think so many details like here are needed

Note on provence - just before submission, we can re-run all cheap cases and possibly also restart all expensive cases (these will be done by the time the review comes back).

Table 2 should not be in the main body of the paper - I think. Maybe an appendix.


5 Verification benchmark suite

