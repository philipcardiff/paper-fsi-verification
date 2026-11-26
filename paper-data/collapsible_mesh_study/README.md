# collapsibleChannel fluid-mesh study (Fig. fig:collapsible_mesh)

Wall-midpoint histories (`F1.dat` to `F4.dat`, the `solidPointDisplacement_wallMid.dat`
of each level) and the driver's `mesh_study.csv`, from the tutorial driver
`collapsible_channel_verification.py --study mesh` (Robin–Neumann, solid 160×8 high
order, Δt = 0.00625 s, tutorial tolerances), unchanged except for the work directory.

- Build: solids4foam `a421f063d`, OpenFOAM-v2512, MeluXina, job 5411606 (2 Nov 2026).
- This run reproduces the paper's numbers, which were first computed on OpenFOAM-v2412:
  maximum error 4.650, 1.425, 0.531 and 0.496% of the reference peak (paper:
  4.65, 1.43, 0.53, 0.50%), and the finest first trough is −0.21610 m.
- Fluid cells: 3 200, 12 800, 28 800 and 51 200, so the refinement ratios are
  2, 1.5 and 1.33.

Figure: `python3 ../figure_scripts/collapsible_mesh_figure.py . REF_dt160.csv REF.csv out.pdf`,
where the two reference CSVs are `collapsibleChannel_oomph_dt160.csv` and
`collapsibleChannel_oomph_reference.csv` from the tutorial's `verification/reference/`.
The paper's figure (`fig:collapsible_mesh`) uses `oomph_extrapolated_dt160.csv` (below) as
REF_dt160.csv, with `0.05 "model form"`; the full command is in `../figure_scripts/README.md`.

## oomph-lib spatial refinement (8 Nov 2026)

`oomph_refine{1,2,3,4}_dt160.csv`: the oomph-lib driver of the tutorial
(`verification/reference/oomph-lib`, oomph-lib ca46a35, macOS, Release) at Δt = 0.00625 s,
t ≤ 1.5 s, with the tutorial arguments and `--refine N`, which multiplies the tutorial
mesh by N (refine 1: 100×16 fluid and 40 beam elements; refine 4: 400×64 and 160).
`oomph_refine1_dt160.csv` matches the tutorial's `collapsibleChannel_oomph_dt160.csv`
to 1e-8 m (rounding).

Wall-midpoint history changes (% of the 0.216 m peak): 0.284 (refine 1→2), 0.083 (2→3),
0.040 (3→4); first-trough changes −0.100, −0.031, −0.015. Observed orders in h ∝ 1/N:
1.24 and 1.07 (history), 1.11 and 1.03 (trough), so oomph-lib converges at first order
in space on this problem.

`oomph_extrapolated_dt160.csv` = 4 f(refine 4) − 3 f(refine 3), the Richardson
extrapolation with p = 1. The (2,3) extrapolation agrees with it to 0.004% of the peak.
Errors of refine 1–4 against it: 0.53, 0.24, 0.16, 0.12% (history); +0.19, +0.09,
+0.06, +0.05% (trough).

solids4foam F1–F4 against it (max error, % of peak / first-trough difference):
4.813/+0.227, 1.592/+0.107, 0.680/+0.046, 0.334/+0.020; orders 1.6, 2.1, 2.5.
Against the tutorial reference (refine 1) the same levels give 4.644, 1.423, 0.531,
0.496%: the 0.5% "floor" is the reference's spatial error.

Figure: `python3 ../figure_scripts/collapsible_mesh_figure.py . oomph_extrapolated_dt160.csv
REF.csv out.pdf 0.05 "model form"` (REF.csv: the tutorial's converged reference, used only
for the peak that normalises the error).
