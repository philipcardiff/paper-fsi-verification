# Hron–Turek FSI3, corrected set-up: case generator

These scripts build the cases of the corrected FSI3 sequence
(Table `tab:hronturek_corrected`, levels 1x′–8x′). They are copied from the
solids4foam branch `verification/hronturek-production`, commit `636d567c6`,
directory `tutorials/fluidSolidInteraction/TurekHron2006/verification/scripts/`.
They are kept here because that branch has not been merged into solids4foam.

- `ht_coupled.py` builds a partitioned FSI3 case from the `TurekHron2006`
  tutorial with the corrections described in its header. `ht_mesh_family.py`
  is the only module it imports.
- `ht_mesh_family.py` generates the fluid blockMeshDict of the ring+box
  family. `ht_coupled.py` fixes the settings to `rows_n 16`, `tip_n 4`,
  `aft_n 60`, `wake_n 70` and `ring_n 22`, all uniform. Level 1x′ has 11 764
  cells, and each level doubles the in-plane divisions.
- `ht_coupled_summary.py` collects the QoIs of the coupled runs.

Settings that the paper's caption does not list (the `ht_coupled.py`
defaults):
- **Fluid solvers:** absolute tolerance 1e-9 for p, U and cellDisplacement,
  with non-final p solves at relTol 0.01. PIMPLE uses 3 outer correctors,
  3 pressure correctors and 1 non-orthogonal corrector. p is solved with
  GAMG.
- **Solid:** PETSc SNES with rtol 1e-9 and stol 1e-12.
- **Coupling:** IQN-ILS Dirichlet–Neumann, outerCorrTolerance 1e-5 with
  `requireAllResidualMeasures yes`, no reuse and at most 30 coupling
  iterations. A step that reaches the limit at the IQN-ILS noise floor is
  accepted and counted.
- **Protocol:** the coupling starts at t = 2 s and the runs end at 7 s.

The solver build was solids4foam `30d6e22f3`, with PR #562 added for 8x′.
`ht_coupled.py` was added to the branch after that commit. The tutorial case
files it starts from are identical at `30d6e22f3` and `636d567c6`, so the
cases it generates are the ones that were run.

Example (level 2x′, solid 4x):

    python3 ht_coupled.py --tutorial <solids4foam>/tutorials/fluidSolidInteraction/TurekHron2006 \
        --level 2 --solid-level 4 --delta-t 5e-4 --cores 16 --out c_2x
