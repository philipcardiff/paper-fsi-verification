# Sources, versions and provenance

## Inputs
* Published fluid-domain surface: `constant/triSurface/fluidDomain.stl.gz` of the solids4foam tutorial `hessenthalerFsi` at commit `517a44921bc9b43ec532e4b30089ef09364e64f6`, sha256 `0be13c049ddf84eb97c7d48af0b6c7680455e9b03342846f0373ad5be3ed815f` (Hessenthaler et al. 2017 data set, CC0, DOI 10.6084/m9.figshare.4141836.v1; the figshare zip, md5 `759248cd617b45682a3f3a1a4592382e`, was downloaded for the pipeline but no MRI value is used here).
* solids4foam: git archive of `517a44921bc9b43ec532e4b30089ef09364e64f6` (development, 2026-11-05; "Add backflowStabilisedOutletVelocity outlet condition (#675)"), unmodified, built with OpenFOAM.com v2512 and PETSc 3.24 (`s4f_case/BUILD_INFO.txt`, build scripts). Note: the verification driver of PR #585 (head `b65d327f1`) lives on a branch diverged from this commit (43 ahead, 582 behind the merge base `2a0eb404b`); the corrected feF1-feF3 results of the main text were not available, so the case was reconstructed from the main-text protocol (deltas in `prereg/ADDENDUM_A6b.md`, applied by `s4f_case/make_case.py`).
* FEniCSx: dolfinx 0.10.0, basix 0.10.0, ufl 2025.2.1, ffcx 0.10.1, PETSc/petsc4py 3.25.6 (MUMPS, hypre, SuperLU_DIST), MPICH 5.0.1, Python 3.11.15, gmsh 4.15.2, numpy 2.4.6; environment lock files in the campaign repository (`environment/`).
* Campaign repository (local, not public): branch `hessenthaler-phaseI`, last commit at the time of packaging in `prereg/` notes and `ledger/`; the FEniCSx code is copied to `fenicsx_source/`.

## Machine and cost
Single node (192 cores, 755 GB), Slurm. CPU-hours reconciled with sacct in `ledger/HESSENTHALER_SACCT.csv` (campaign total in the file); the solids4foam reconstruction cost 14.5 CPU-h for 0-0.75 s on 16 cores (19.3 CPU-h per simulated second; an earlier 0-0.1 s pilot of 2.1 CPU-h is separate).

## Meshes (not included; hashes)
FEniCSx: `data/checks/MESH_SHA256.txt`. solids4foam reconstruction: `s4f_case/fluid_polyMesh_SHA256.txt`, `solid_polyMesh_SHA256.txt`, `log.checkMesh.*`, `log.cartesianMesh` (cfMesh of OpenFOAM v2512 `cartesianMesh` after `surfaceFeatureEdges -angle 45`).

## Preregistration chronology
`prereg/PREREG.sha256` lists the sha256 of the preregistration and each addendum. A6 was committed at 20:09:40 and A6b at 20:10:47 (the build job had started at 20:10:04; no case had been run); both precede every solids4foam case run (first run 20:16:55). A7 (temporal control) was committed 1 s before its job and was written knowing the solids4foam Fy peak, so it is not blind; its criterion (i) covers 0.1 <= t <= 0.75 s (16.8 % at 0.1 s), the supplement headlines t >= 0.2 s (7.0 %) and states both.
