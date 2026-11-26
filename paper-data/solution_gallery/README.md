# Solution gallery (fig:solution_gallery)

`figures/solution_gallery.pdf`: one representative solution of each benchmark,
drawn from exact solutions or from this paper's own partitioned runs (no image or
field of the companion quasi-monolithic paper is reused).

## Scripts

- `make_gallery.py`: assembler (layout, true aspect ratio, colour bars); drawn at
  the text width `figstyle.TEXTWIDTH` (15 cm = 427 pt) and included at natural size.
- `panel_*.py`: one script per panel; each also runs on its own and writes a
  one-panel preview (`python3 panel_X.py <data> out.png`).
- `foamio.py`: reads a reconstructed solids4foam multi-region case with
  `vtkOpenFOAMReader` (VTK Python bindings; no ParaView), cuts a plane or
  extracts a patch, returns a triangulation with point fields.
- `fsi2d.py`: planar fluid + deformed-solid snapshot (speed/vorticity,
  streamlines, solid displaced by `pointD`, undeformed outline dashed).
- `common.py`: style (`figstyle.use()`, with matplotlib's `cmr10` as the
  Computer Modern fallback), colour bars, scale bars, streamlines.

Requirements: Python 3 with numpy, matplotlib (3.10) and vtk (9.x).

```bash
python3 make_gallery.py ../../figures/solution_gallery.pdf --png preview.png \
  --womersley <s4f>/tutorials/fluidSolidInteraction/womersleyTube/constant/womersleyProperties \
  --collapsible <coll_F2> --fsi3 <fsi3_g2> --tube <tube_L1> --mok <mok64> \
  --cavity <cav_L2> --blob <blob_x2> --hessenthaler <hv_rF1>
```

## Data sources (26 Nov 2026)

All reruns used the isolated macOS build of solids4foam `206437983`
(OpenFOAM-v2512, Apple Clang `-O3`), the build of the paper's cavity and 3-D tube
results. Only the fields needed were written; the case directories are not archived
in this repository.

| Panel | Source | Snapshot | Check against the paper |
|---|---|---|---|
| (a) ring | exact: plane-strain continuum ring + potential flow, Chebyshev collocation (as `ringAddedMass/verification/scripts/ring_exact_frequencies.py`), strong level | maximum displacement | omega = 0.768654 rad/s (paper: 0.768654) |
| (b) Mok | rerun: tutorial at `6e69ddb6c` (the provenance-table commit), 64^2 member of its mesh study (64x8 high-order solid, dt 0.05 s, IQN-ILS), serial, written every 0.25 s, stopped at t = 28 s | t = 24.5 s (peak) | midpoint peak 0.2871 m (paper 64^2 peak 0.2870 m) |
| (c) Womersley | exact: travelling-wave solution, coefficients from the tutorial's `constant/womersleyProperties` (written by `womersley_exact.py`) | omega t = 140 deg | wave speed 0.872796 m/s (paper: same) |
| (d) collapsible | rerun: tutorial at `a421f063d` configured by its verification driver as mesh-study level F2 (solid 160x8 high-order, Robin-Neumann, dt 0.00625 s), serial, to the first trough | t = 0.46875 s | wall-midpoint u_y -0.21591431 m (`collapsible_mesh_study/F2.dat`: -0.21591434 m) |
| (e) FSI3 | this paper's Gmsh level-2 run on UCD Sonic, `/scratch/pcardiff/fsi3gmsh/cases/g2_np32` (read-only); meshes + U, p, D, pointD of t = 5.500711216 copied to the Mac and reconstructed there | t = 5.50 s (tip u_y 0.0340 m; cycle maximum 0.0368 m at 5.512 s) | tip displacement from the case's own postProcessing |
| (f) 3-D tube | rerun: paper level 1, corrected meshes (`robin_backward_mesh1` of the 3dTube mesh study, Robin-Neumann, backward, 8 ranks) to t = 5 ms, written every 1 ms | t = 5 ms | same case set-up as the paper's level-1 run |
| (g) cavity | rerun: skew-corrected level-2 case (`cavgap/skew_L2/.../aitken/mesh2`), serial, written at t = 200 s | t = 200 s | tip u_y -0.231268 m, identical to the stripped original run |
| (h) Hessenthaler | this paper's partitioned rF1 run on Sonic, `/scratch/pcardiff/hv/P_rF1_LUST_b02_mwp` (read-only; IQN-ILS, build 206437983 + PR 676 outlet BC); t = 45 copied and reconstructed on the Mac | t = 45 s, plane x = 0.3 mm | flap tip u_y 16.65 mm |
| (i) blob | rerun: tutorial at `1b5ab0963` (PR #694, finite-strain neo-Hookean solid), mesh factor 2 of its verification study (1 920 fluid cells, dt 0.01 s, IQN-ILS), serial | t = 3 s | top-point u_x 0.139127 m (paper level 3: 0.139075 m, 0.04%) |
