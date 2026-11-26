Feel free to push back on any of the comments below.

5.3 Classical FSI benchmarks revisited

"heavier and softer flag"
i prefer "more flexible" than "softer" since softeness more correctly refers to strength than stiffness

"are category B."
I wonder if we need to continually re-define Category B (or A or C etc) since readers will quickly forgot how we defined them.
I wonder even if we can use names for the categories that are self descriptive - i cant think of suitable ones but maybe you can suggest ones.

"Its statistics are taken over t= 9–
10 s rather than five periods before t= 6 s; the definitions are the same."
What definitions?

"The Featflow tables are among the most carefully
documented FSI references available, and the Turek–Hron problems owe much of their influence
to the multi-level results published with them; the points below concern how far those results can
bound the error of another code, not their value as a benchmark. "
I agree with these comments (which are complimentary); in contrast, the earlier comments "These benchmarks are widely used, but their published values carry incompletely quantified error,
are inconsistent between sources, or both. " could come across a bit insulting to some of those authors; I agree with the text as written but I wonder if we need to be careful in how we phrase it or possibly soften the earlier comments (or sandwich them with "compliments"

"and the lift amplitude rises by ≈4 N/m per level without settling"
I wonder if it would be better to give 4 N/m as a relative or percentage term since the absolute value means nothing without a reference.

"over its 7.9 periods the per-period lift amplitude drifts from 157.7
to 153.5 N/m"
if i remember correctly, the values are going in one direction possibly suggesting some sort of bias in the featflow results, e.g. an error accumulation due of some kind - I am not sure if we should comment on this (i remember we some something like this and it was due to the use of a velocity mesh motion formulation and was resolved by changing to displacement or tigthening the tolerances reduced it)

"The detailed Featflow tables [57, 58], used here,
differ from the frequently quoted summary values of Turek and Hron [1] (Table 8), most strongly
in the FSI3 drag amplitude "
why did we use their earlier results? and which ones are more correct?

"A comparison with the summary values
would misattribute this 18% to the code under test. As a check, the driver reproduces the tabulated
statistics from the published Featflow time histories to within 0.6% (FSI3; 0.8% for drag amplitude,
0.3% for the frequencies) and 0.3% (FSI2)."
I am confused by what you are saying here - I am not sure what you are saying.

"AUTHOR NOTE: Waits on Philip (needs the Springer chapter; no open ..."
YOu can find this on my mac at "/Users/philipc/Documents/Papers/Literature - sorted" titled as "Turek\ -\ 2006\ -\ Proposal\ for\ Numerical\ Benchmarking\ of\ Fluid-Structure\ Interaction\ between\ an\ Elastic\ Object\ and\ Laminar\ Incompressible\ Flow.pdf"

Table 8: i am a bit confused about what is shown here - are you comparing the older and newer Turek and Hron results?
If so, why? for what reason? what are you trying to show?

Numerical model for FSI3: we should briefly note why we switched to an unstructured mesh

"(PIMPLE 3/3/1);" what do 3/3/1 mean here? this should be clear to the reader

"A uniform block-structured family
(1x, 2x, 4x: 5 336, 21 344 and 85 376 fluid and 630, 2 520 and 10 080 solid cells; ∆t= 0.001, 0.0005
and 0.00025 s; IQN-ILS, coupling tolerance 10−5) is used for the FSI3 coupling studies below"
this is inconsistent with the description just above it which says we use an unstructured mesh for FSI3 - make this section clear and consistent and if we have two meshes for FSI3 we need to explain why - I think the new unstructured results are meant to replace the old structured ones, right?

"Table 10 Turek–Hron FSI2: relative errors "
is "error" the right word here? or should it be difference?

It seems that FSI3 results are given first then FSI2 then FSI1 - surely it would make sure sense to start with FSI1 and then FSI2 and then FSI3?

"FSI1 with the high-order solid."
since we are trying to reduce the pape length, should we just only present the higher order results? and only brielfy comment on the 2nd order ones?

"Independent numerical comparisons."
Given that we already have FeatFlow references, I wonder if this section should be shortened dramatically unless there is something we particularly want to highlight.

