Feel free to push back on any of the comments below.


5.3.2 Three-dimensional pressure-pulse tube

Since this is C-weak, should we present it later in the paper? i.e. finish with C-weak, after we have given B and C-strong?

"The solid here is Hookean (small strain), whereas most published
calculations use St Venant–Kirchhoff (peak hoop strain ≈3%);"
It seems strange that we say most people do X but we do Y - why?
I guess the strains are small so the difference is small (but can still be large relative to discretion error on finer grids).
If we dont haev a good reason (e.g. if the reason is that the tutorial was set up this way), then we should be more careful about the phrasing here and just acknowledge the variants in literature and state the variant we look at here.

", and non-conformal (AMI) interface transfer."
We could easily have used a conformal interface here so maybe we should add a comment saying that we chose a non-conformal interface on purpose to demonstrate/test it.

"The reference is category C-weak."
I suggest we phrase it more along the lines of "hence, this reference is categorised as C-weak" just to make it clear that we are performing this categorisation, and also to briefly explain why if it is not immediately clear.

"The point-A QoIs on L3 ..."
L3 has not been defined here yet. For all test cases, check for similar mistakes where meshes or terms are mentioned before they are defined - fix them.
Or maybe you already stated before this that L1, L2, L£, etc represent the grid levels - of so, ok.

Fig 8: "smoothed and resampled to 0.1 ms, so the markers are not solver samples" - why are they smoothed and resampled? Why not present the actual extracted data? is it because noise has been added because of smoothing? if so, we should ideally comment on any error introrduced by the smoothing, e.g ideally it would be the same as or smaller than the extraction error.

Most images are collected on image pages, resulting in them being somewhat far from the accompaning text; i understand this is good for reducing the page count but it is not good for the readability. See if yuo can strike a balance where ideally most images are close to where they are referred to.

"The point-A QoIs on L3 are collected in Table E12 (Appendix E)"
it is a bit strange that the very first sentence for the reference section points readers to the appendix, especially in this case since Table 12 in the main body also has key values. It seems more appropriate that the appendix tables (with more complete or comprehensive results) are referred to after the results in the main body of the paper are introduced.

Would it make sense to run L4 on the tube case to see if monotonic and aymptotic convergence can be achieved? Maybe it would be too slow.

"Measured on a variant of the meshes whose wall arcs are defined to fewer digits"
this seems like it could be confusing for a reader - they will think "what?"
ideally all our conclusions should be presented on one consistent set of meshes, and where needed we could note this distinction only if needed.

"Platforms. This case exposed two defects ..."
I feel like this is at most a lesson for the learnings appendix, with maybe a single sentence in the lesson section of the case in the main paper.

"Implicit Euler" would it be more appropriate to refer to this as BDF1 (and 1st order implicit Euler in the first instance), since implicit Euler could possibly also refer to higher order variants (maybe not...)

"also differ in wall model, pressure switch-off "
what are these and do we need to add a brief comment to make it clear what they are and what that reference did with them?


5.3.3 Mok lid-driven cavity

The Reference paragraph is great here - it exactly highlights the need for this paper.

"The reference is an independent solution ..."
I wonder if it would be clearer here to say "The reference taken here is an indep...." just to be clear what "the reference" is here, since a reader may be confused.

Is the problem definition missing here? e.g. the informaiton required for a reader to set up the case in their own code. I am wonder the same for the other cases too.
It seems strange that much of the problem definition appears first when talking about the Kratos solution - surely is should be defined first and then for the Kratos reference, we only need to state particulars to it.

"the solids4foam meshes are uniform"
in previous cases, we seem to refer to our own (solids4foam) method by the particular scheme, e.g. IQN-ILS with TL solid etc - here seems to eb the first place where we are referring to solids4foam. I actually dont mind using solids4foam throughout (even from a branding perspective) although it is probably clearly to talk about the scheme.

"; Kratos needs no such change."
presumably Kratos does not use PIMPLE so it has a different solution algorithm - i would remove this part of the comment or qualify it with "but it uses an XXX solution algorithm"

"The study used the legacy coupling test"
i dont think we defined the "legacy" approach so just state what it is either here or earlier, ideally without using the word legacy since that sounds like something old and may have had problems

"AUTHOR NOTE: Pending results: the cause of the strict-test stalls ..."
We should resolve this point. It seems like the strict tests are not needed here so maybe we dont need to fix them here - maybe we just state that some time steps were accepted on only one residual and when we compared using the strict checks (for some meshes and time steps), it did not change results much. What do you think?


5.4 Provisional and in-development benchmarks

"Pending Philip, with Z. Tukovi´c: which constitutive law Tukovi´c et al. [12] used"
Zeljko said that neo-Hookean was not yet implemented in his code at that time so it could not have neen neo-Hookean so it was either St Venant Kirchhoff or linear elasticity - he suggested he may have accidentally used linear elasticity - we could check and see if that were to make sense. Otherwise, we can just leave it unexplained in the discussion.

"so the monitoring point must be stated."
this seems like a strange thing to say in a paper -> maybe you are trying to say that the location of the point through the thickness does have a small but measureable difference.

"of the code under test"
strange way to refer to our current code - we need to be consistent in the way we are referring to our solids4foam solver(s) in this paper.

"Two further solutions agree as closely but are not independent: a monolithic finite
element solution computed by the authors of this paper with a separate implementation (FEniCSx; Supplementary Material),"
Surely the FEniCSx solution is just as independent as the SU2 one? They are both different codes from ours and both run by co-authors? maybe I am missing something.

"the independent solutions expose an interpolation error that the four-level sequence alone does not"
What are you saying here?
if our solver has a problem, we should fix it; if this was something we fixed then is this the best way to explain it here?

"COMSOL (two-dimensional finite elements, by E. Tandis), on"
we said co-author for other results like this so why name E Tandis here? we should not

Sections 6, 7 and 8: given the length of the paper, I wonder if these 3 sections should be merged?