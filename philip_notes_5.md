Feel free to push back on any of the comments below.

5.2.4 Flexible-bottom cavity

"the outlet has p = 0 "
do we also need to state that the normal gradient of velocity is zero for the problem to be fully defined?

5.3 Partly quantified references (category C)

"Turek–Hron FSI2, the other category-C case, is treated with its family in Section 5.2.2."
To me, I think it would be better to place FSI2 here, although of course we can refer back to Section 5.2.2 for much of the case description - and in Section 5.2.2, we can state that FSI2 is not discussed there since we categorise it as Partly quantified. I think it is more important to keep cases in their correct classification order then keeping families together - otherwise a reader may incorrectly assume FSI2 is a different category - feel free to push back.

"is the only steady three-dimensional external-flow case in the suite"
is it considered external flow since it is in a channel? I am not sure if channel flow is considered internal or external.

"on the mesh of the distributed
case, including it changes the tip displacement by +1.7%, and on that mesh refined twice by
+0.35% (Section 3.5)"
what is the distributed case?
I guess the point here is that we are trying to provide evidence that ignoring that viscous term in the pressure boundary condition has only a small effect on our quoted results, right? If so, I think we should be more careful in our phrasing here since a reader who may not be concentrating my incorrectly see 1.7% and think it is large relative to later quoted confidence values.

"; the comparison was made after the sequence was characterised."
I dont know what this means

"so an independent code also reaches the
reference values only with the 0.2 m/s peak."
here i think it would be useful to briefly note again what this confirms i.e. inflow velocity must be X not Y.

", and the category is unchanged."
I guess a finer mesh in this case would possibly allow the quantification of this case to increase to category B, but that mesh would be expensive. I think it would be useful to comment on this so authors understand why we consider it C and how it might become B - in fact, I expect it can become B (and probably will in the coming 5 years) with increasing computing power and moves to faster implementations (like GPUs or parallel in time).


5.4 Unquantified references (category D)

"and for the Hessenthaler flap a separate
comparison with measurements."
I wonder should we add here something like "which as noted is not evidence of verification" to be very clear that we consider validation and verification to be different.

"The wall-front speed is"
what is the wall-front? is it clear in the text?

For this 3d pressure pulse tube, it comes across as strange that we ran this case as linear elastic when all 3 references use St Venant–Kirchhoff. I suggest you ask Amir to re-run the Kratos solution as St Venant–Kirchhoff too, and we can do the same. Until we have those results, we can keep the paper referring to linear elastic.

"A quarter of the tube is computed with symmetry planes,"
I guess we could have used axisymmetry here, right? If so, we should state that we are purposely using 3D, presumably just to test the 3d part of the code, e.g. curved interface in 3d with AMI.

"the remaining differences
in the wall law, the switch-off and the meshes are not quantified, "
what is the "switch off"?
and what do you mean by meshes here? since we did a mesh study.
we should be clear about what we did not check and if we think it is important so readers can understand if they are important or not

"(3dTube/verification/reference/ of the solids4foam tutorial)."
do not refer to solids4foam tutorials in the paper.
if you want to refer to data, it should be the paper supplementary data or the appendices.
check all cases for this.

"The raw traces were not kept, so the smoothing error cannot be measured on them directly;"
remove this statement and instead just state the bounds estimates on the error - in fact, this paragraph can be shortened appreciably - the goal of it is to just briefly state our estimated uncertanty on the extracted data - which in any case is much smaller than the differences from their data to our converged in time results.


5.4.2 Hessenthaler flap (Phase I)

"and the deformed centreline measured under the Phase I flow"
do we need to be clear here that this is the MRI experiment curve? as it stands, this may be unclear, e.g. readers may think this is our result.

It would be good if we could show the geometry as 3d for this case as the figures may incorrectly suggest a 2d case, e.g. replace or complement Fig 11 (a) with a 3D perspective view of the geometry.

"the solid uses the total Lagrangian displacement formulation on 12 ×8 ×65 hexahedra "
we should comment on the solid mesh error -> i guess it is small here right?

", where the measured plane flow rates exceed the pump rate by 8–10% "
here we are essentually pointing out an error or inconsistency in the experimental measurements - but as phrased this point may be missed (eg the measured plane flow rates could incorrectly be interpreted as from our model); I suggest this is phrased more clearly.

"Lessons. For a complex three-dimensional geometry, three affordable non-nested meshes give a
mesh sensitivity"
Ee mention "affordable" here, which is fine but we have not actually mention run-times anywhere - i wonder if we should add them to an appendix? Then, e.g. for this case, it becomes clear what we mean by affordable since L3 takes many days on HPC.

5.4.3 Blob in treacle

"An elastic half-cylinder in a highly viscous channel flow ..."
FYI, the name "blob in treacle" was introduced by me when adding the case to solids4foam - this name was not used by Liu et al - I wonder if we should note that here, e.g. "refered to here as "blob in treacle".

Section 6

In general, I really like section 6 as written: it concisely gets across the key learning points of this paper and I think it will be genuinely valuable to the community.

"For one independent cost point, FEBio’s monolithic
solver took 2 033 s"
I like that we mention this but i wonder if it is a bit unfair -> it obviously favours our solver and may come across unfavourably to reviewers/readers (as boasting of sorts). If we were to mention this, then it would seem fair to give some indication of costs for the other codes and cases too - I am not sure comprehensive details of this would be valuable here. Possibly, high level ranges across cases and codes could be useful with the caveat that they we were on different machines (where time to solution could change by factors of 2x or more).

"Its backbone is the three category-
A problems, which alone measure the total error directly; the category-B problems extend it to
large deformation within"
I guess the MMS case in category A is already large deformation right? I wonder if this text needs updating.

"The suite covers incompressible, laminar, ..."
This paragraph is fine, although it almost suggest that this suite is unsuitable for other scenarios e.g. like testing an immersed boundary code. I suggest this paragraph is rephrased or extended to say that the suite here is in fact suitable for other solver algorithms like immersed methods, and is naturally extendable to other physics like compressible or turbulent flow although each of those cases would like need their own specific cases too.

"6.7 Executable verification infrastructure"
I dont see much value to this section - it seems to repeat what was already said. I suggest merging this into previous section and extending/generalising those sections as needed.


7 Conclusions

"An exact reference at finite deformation qualifies this." this seems like a strange heading - i am not sure i understand what "this" is here

"A closed-form manufa..." this is a very long sentence

"with three unreleased corrections"
I suggest removing "unreleased" here since it is not relevant here and in any case the corrections are all public even if they are not yet on a solids4foam release.

"Independent converged references where none exists, above all for the 3-D pressure-
pulse tube,"
Since the pressure-pulse tube is axi symmetric, surely we should create and run an axisymmetric version => I even wonder if we should do this before we submit the paper since we could use the axisymmetric cheap cheap case to provide greater understanding of the 3d case. What do you think? A sub agent could do it

"Wider coverage: compressible and non-New..."
And immersed approach => i mention this since this is something we are actively working on in solids4foam and we will use this suite to test it.

"Higher-order discretisation, fo"
Do we need to briefly note here that we used higher order for solid but the overall procedure is still limited ot 2nd order, so here we are proposing higher order for solid+fluid+coupling.


Acknowledgments

Grant number "16/RC/3872" can be removed but keep the others

