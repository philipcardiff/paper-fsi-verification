# Addendum A1 to the preregistration (2026-11-11, dated before any accepted solid-only result)
The Gate B2 Newton criterion "|r|/|r0| < 1e-9" is replaced for the solid-only static problem by
"|r| <= 1e-7 |f_ext(lambda=1)| lambda  OR  |dx|_inf < 1e-12 m", because the axially stiff / bending-soft flap has a
roundoff floor of ~2e-8 relative (observed in the first debugging runs: quadratic convergence to 3.5e-14 then stagnation). The
other thresholds (prediction interval 29.20-29.60 mm, within 0.5 mm of 29.50, monotone 3 levels) are unchanged.
Also (Gate A): the mesh-quality criterion is the gmsh minSICN of the tets; the registered threshold 0.05 is unchanged.
