> **ARCHIVED WORKING REPORT (verbatim except path/host sanitisation).** Written during the investigation; file paths such as `runs/`, `tools/`, `cavity_agent/` refer to the original working tree, not to this package (see `docs/INDEX.md` for the package equivalents). Later reports supersede earlier ones where they disagree; known numerical/wording inconsistencies and their resolution are listed in `provenance/INCONSISTENCIES.md`. Nothing here claims a verified FEBio result.

# CAVITY B0a: plate-only results (FEBio v4.13, solid module, StVK)
Cases (Slurm 13809, 13810; 2 cores, main; ~1 s each, 51 MB; exit 0 and "N O R M A L TERMINATION" in the `.log`): plane strain (u_z = 0, one hex20 layer), x in [2,6], y in [-0.4,0], both ends clamped, St.Venant-Kirchhoff E = 500, nu = 0.4, uniform **follower** pressure 1.315 Pa on y = 0 (pushing down), 10 load steps, 3 Newton iterations per step, no max_refs events. Generator `tools/make_cavity_plate_feb.py`; analysis `tools/cavity_plate_analyze.py`; data `runs/cavity/B0_plate_*`.

| mesh (hex20) | nodes | u_y(A) FEBio | reference (Q9, same mesh label) | difference |
|---|---|---|---|---|
| 20x2 | 393 | -0.237313 m | -0.23887 | -0.65 % |
| 40x4 | 1343 | -0.239985 m | -0.24047 | -0.20 % |

- Load-displacement is monotone and smooth (10 steps, A: -0.0302 ... -0.2373 at 20x2; -0.0307 ... -0.2400 at 40x4).
- Two-level estimate (assuming p = 2, not verified with three levels): FEBio extrapolates to -0.24088 m; the Q9 references extrapolate to -0.24100 m (same assumption): agreement ≈ 0.05 %. The remaining 0.2-0.65 % at fixed mesh is consistent with hex20 (serendipity) vs the Q9 (Lagrange) element and not with a material-law difference: the other large-strain laws in the manuscript give +4.8 to +5.7 % (neo-Hookean/Hencky) and linear -0.288, i.e. an order of magnitude larger. **StVK equivalence is supported numerically**, with the caveat that the reference is another numerical solution, not an exact one.
- Point A is the top-surface node at (4, 0): exact node by coordinate; uy history is read at that node.
- Jacobian J (volume ratio per element point) at the end of the loading: 0.985-1.014 (20x2), 0.972-1.024 (40x4): positive, mesh valid.
- **Reaction-force equilibrium was NOT verified numerically.** The `Rx/Ry/Rz` node log returned exact zeros for the clamped nodes in this solid-module run (cause not investigated; likely the reaction vector is not populated for these BC types in the log path). The applied vertical load is known exactly by construction (a follower pressure on a surface whose end points are clamped has projected x-extent 4 m: F = 1.315 x 4 x thickness), but the balance with the supports was not measured. Not blocking for the plate check, to be revisited in the coupled model through an independent force measure.
- Pressure direction: verified by the sign of u_y (downward) and its magnitude.
- Cost: negligible. Budget used: 2 of 8 authorised cavity simulations.
**Gate B0a: PASS (with the reaction-force caveat).**
