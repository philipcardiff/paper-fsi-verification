# Known limitations, unresolved questions, reopening criteria

## Limitations of the evidence
1. **No coupled run.** Every FSI1 statement about the coupled problem is NOT ESTABLISHED. The flag deflections quoted for a nodal-pressure load (11.6 mm) and an element-mean load (1.56 mm) are Euler-Bernoulli analytical diagnostics from rigid-flag loads, not simulations.
2. **Stage C is one mesh level (8,330 cells), one K (1e6), one time-step/ramp (0.5 s / 2 s), 80 steps.** The pre-registered stop rule ended it after level 1; the behaviour of the parity mode under refinement, K, dt, ramp and run length is NOT TESTED (hypothesis, UNRESOLVED: the mode is undamped and accumulates from the discrete divergence of the pseudo-transient).
3. **The two pinned junction values** (-81.60 and -123.87 Pa) were derived from neighbouring wall nodes that carry the parity mode themselves (uncertain by 10-40 Pa); the pin-value sensitivity was NOT TESTED. The pins do not cause the parity mode (cylinder amplitude identical with and without pins).
4. **Force numbers are reconstructions.** F1, F-actual, F2, F3, F3b, F4 are different post-processings of the same stored nodal fields; none is a FEBio output (reaction forces are zero in fluid and solid modules for fixed DOFs, and no coupled run produced an applied FSI load). F4 relies on my reading of the solver's weak form (interior residual closes to 0.1 %, 1.3 % at the worst node).
5. **Structural stage:** equilibrium at the planned 1e-6 level not available; the converged hex8 value is 0.30 % (u_y) beyond the Featflow value and the hex20 nt 4 value +0.08 %: the offset is unexplained; the Featflow flag geometry is read from its definition page and the tutorial, not from its mesh files.
6. **References:** Featflow tables only (retrieved 2026-10-31); the Turek-Hron 2006 paper tables for CFD1/FSI1 could not be retrieved, so any statement that they equal the Featflow numbers is UNRESOLVED. Featflow's own discretisation error of the CFD1 and CSM values is about 0.01 % or below (level differences), negligible against the effects reported.
7. **Not solids4foam evidence.** This package says nothing about the solids4foam FSI1 results; it uses its tutorial only for the case definition.
8. Data not shipped (110 MB: xplt files, multi-step logs): listed with sha256 in `provenance/retained_hpc_outputs.csv`. The final-step state of C1u and C1p is shipped and suffices to recompute all stage-C numbers; the time history of the drag (F1) over the last 25 % of the run (steady to 1e-11) is in `results/stageC_raw/*/F1_history.txt`.

## What would reopen the question (none approved or implied)
- A level-2 (and level-3) fluid-only run with the same set-up showing the odd-even amplitude falling by >= 1.5 per level toward <= 0.048 Pa and the nodal-pressure force approaching the reconstructions within 1 %;
- sensitivity runs (K x 10, dt/ramp, run length, pin values) establishing whether the amplitude is history dependent;
- a coupled level-1 FSI1 run only if the first two items pass; its comparison would still need three levels for any VERIFIED AGREEMENT;
- a FEBio-side change or option that makes the interface load insensitive to the nodal row-parity mode (not proposed here; an earlier local stabilisation test failed its own mass-conservation gate and is not part of this package).
