# mokLidDrivenCavity: legacy vs strict coupling test (Section sec:mok)

Membrane-midpoint histories (`<variant>-<member>.csv`: time, u_y in m) of the
group-B definition, run with the tutorial driver through `mok_var.py`. The
wrapper changes only the work directory, the member selection, an extra
`mesh128` member (Δt = 0.01875 s) and, for the variants, `fsiProperties.iqnils`:

| variant | coupling test | couplingReuse |
|---|---|---|
| legacy | min(r1, r2) (tutorial) | 5 (tutorial) |
| strict | requireAllResidualMeasures yes | 5 |
| legacy0 / strict0 | as above | 0 |

Tolerance 1e-5 throughout, IQN-ILS, high-order solid 8 cells through the thickness.
Build: solids4foam `a421f063d`, OpenFOAM-v2512, MeluXina jobs 5411716 and
5411733 (2 Nov 2026).

`statistics.csv` gives the peak, trough and mean over t = 20–70 s, computed as
in the driver; it is written by `../figure_scripts/mok_figure.py` run on this
directory (command in `../figure_scripts/README.md`), whose `SERIES` table
selects the files plotted. The figure
uses legacy0 for 16², because the legacy (reuse 5) 16² run stopped at t = 60.2 s
(the solid Newton solve diverged). At 32² the two reuse settings differ by at
most 0.025%.

Completed: legacy 32², 64², 96²; legacy0 16², 32²; strict 32², 64², 64²/Δt 0.025;
strict0 32². Stopped because IQN-ILS stalled at the iteration limit (strict
test): strict 16² (t = 54.8 s), strict0 16² (64.4 s), strict 96² (14.8 s),
strict 128² (14.0 s), strict and strict0 64²/Δt 0.1 (19.4 s and 19.8 s).

Result: strict vs legacy differ by at most 0.024% (trough) on 32² and 64²,
against 0.3–2.3% between meshes. The legacy runs reproduce Table tab:mok
(first computed on OpenFOAM-v2412, macOS) to the digits given.
