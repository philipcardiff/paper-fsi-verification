# Hessenthaler geometry figure (`figures/hessenthaler-geometry.pdf`)

- `section.py` cuts the published fluid-domain surface (solids4foam
  `tutorials/fluidSolidInteraction/hessenthalerFsi/constant/triSurface/fluidDomain.stl.gz`,
  development branch) at x = 0, giving `section_x0.npz` (segments per patch, mm).
- `rF1_x0.vtk` is the rF1 fluid mesh (312 180 cells; copied read-only from
  UCD Sonic `/scratch/pcardiff/hv/P_rF1_LUST_b02/constant/fluid/polyMesh`,
  20 Nov 2026) cut at x = 0.3 mm with OpenFOAM v2512 (`postProcess -region
  fluid`, `surfaces` function object, `cuttingPlane`, no interpolation,
  bounds y = ±41 mm, z = −32 to 102 mm); the plane is offset from x = 0 so
  that it does not coincide with octree cell faces. `checkMesh` on the copy:
  312 180 cells, maximum non-orthogonality 37.5°. It replaces the earlier
  `feF1_x0.vtk` slice, since the verification now uses the rF meshes.
- `phaseI_centreline.csv` is the measured Phase I deformed centreline
  (Hessenthaler et al. 2017, CC0; copied from the tutorial's `validation/reference`).
- `fig.py` draws the figure (matplotlib, paper style of `../figstyle`):
  `python3 fig.py --out <repo>/figures/hessenthaler-geometry.pdf` (it reads the
  data files from its own directory, or from `--data DIR`). Panel (a) is the
  section at x = 0; the mesh of panel (b) is the cut at x = 0.3 mm described above.
