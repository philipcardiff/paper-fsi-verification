Thanks, Philip. I followed up on your suggestion with one strong-M2 halved-amplitude run and a reanalysis of the stored frequency histories. The results strongly support an approximately amplitude-squared contribution to the remaining FEBio ratio-error floor.

**1. Halved-amplitude test**

I reduced the initial velocity from V0 = 1.5e-3 to 0.75e-3 m/s, keeping the remaining physical and numerical settings unchanged. The FEBio run completed all 400 steps successfully.

The wet/dry ratio error changed from +0.0399% to +0.0097%, a reduction by a factor of approximately 4.1. The corresponding two-point amplitude exponent is approximately 2.04.

A two-point extrapolation to zero amplitude gives a ratio error of about -0.0004 +/- 0.0003%. This is consistent with the disappearance of the observed finite-amplitude bias, although it is not an independent verification of the exact zero-amplitude limit.

The wet frequency changed from 0.773788 to 0.773554 rad/s. The n=2 pressure harmonic decreased proportionally with the imposed amplitude.

The original dry solution was reused. Its harmonic content is small, but the dry fundamental frequency was not independently recalculated at half amplitude.

**2. Cross-regime scaling**

The existing weak, moderate and strong results show that the inferred squared-frequency bias, when normalised by the squared oscillation amplitude, is approximately consistent across the three added-mass regimes.

This provides additional support for an amplitude-related contribution to the error floor. Only strong M2 was directly rerun at half amplitude.

**3. Start-up and frequency-fitting effects**

Reprocessing the existing histories using only the last two periods gives essentially unchanged errors:

- Weak M2: approximately +0.0041%
- Moderate M2: +0.0072% to +0.0073%
- Strong M2: +0.0399% to +0.0401%

The fitting-window spread in the strong case is approximately 0.0003 percentage points, much smaller than the roughly 0.030 percentage-point amplitude-dependent change.

Start-up contamination therefore does not appear to explain the observed floor.

A small wet-only third harmonic also changes with amplitude in a manner consistent with nonlinear behaviour. However, it does not identify the responsible fluid mechanism.

**4. Interpretation and remaining limitations**

These findings support a finite-amplitude effect in the wet response relative to the infinitesimal-amplitude analytical reference.

They do not identify whether the effect originates in fluid convection, moving-domain geometry, interface no-penetration enforcement or another aspect of the wet formulation.

Only two amplitudes were directly compared. A corresponding amplitude study in solids4foam has not been performed.

The existing original-amplitude FEBio results remain unchanged. I would retain them as the primary benchmark evidence and present this experiment separately as an explanation of the residual ratio-error floor.

The test used four cores, took 741 s wall time and consumed approximately 0.82 allocated core-hours.

The PR remains independent evidence rather than a Category A verification claim.

Would you prefer the manuscript to mention the finite-amplitude sensitivity explicitly, or keep it documented only in the accompanying independent FEBio evidence?

Thanks again for suggesting this test.

Generated with Claude Code.
