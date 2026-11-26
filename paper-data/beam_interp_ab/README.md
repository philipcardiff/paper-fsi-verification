# beamInCrossFlow: A/B check of the interface point transfer (PR #546), 2026-11-26

Closes the todo in the beam appendix of `article.tex`. The beam runs used
solids4foam `aee0c8e35`, built with `-O2` on MeluXina, which predates the
interpolation-weight fix of PR #546 (merge `978ffc8dd`).

## Result

The defect cannot affect the beam, for two independent reasons:

1. **The beam does not use the barycentric weights.** The interface is
   conformal, and every campaign case uses `interfaceTransferMethod directMap`
   (`constant/fsiProperties.iqnils`, which `richter_case.py` copies
   unchanged). The run log confirms this: `Selecting
   interfaceToInterfaceMapping directMap`. The defective code
   (`triangle::pointToBarycentric` in `amiInterfaceToInterfaceMapping.C` and
   `amiZoneInterpolation.C`) runs only for AMI transfer and solid contact.
   PR #546 does not touch `directMapInterfaceToInterfaceMapping`.
2. **Even with AMI, no triangle would fall below the threshold.** The old
   weights used OpenFOAM v2412's `triangle::pointToBarycentric`. It returns
   (1/3, 1/3, 1/3) when `d00*d11 - d01^2 = 4A^2 < SMALL = 1e-15`, that is,
   for a triangle area `A < 1.58e-8 m^2`. The triangles are (face point i,
   face point i+1, face-point average), i.e. a quarter of each interface quad.
   Computed from the `blockMesh` output of the case generator on all four
   levels, for the solid and the fluid interface (`triangle_areas.txt`,
   `scripts/tri_areas.py`):

   | level | interface faces | min face area (m^2) | min triangle area (m^2) | min 4A^2 | ratio to SMALL | below threshold |
   | --- | ---: | ---: | ---: | ---: | ---: | ---: |
   | 0 | 256 | 3.125e-4 | 7.81e-5 | 2.4e-8 | 2.4e7 | 0 |
   | 1 | 1024 | 7.81e-5 | 1.95e-5 | 1.5e-9 | 1.5e6 | 0 |
   | 2 | 4096 | 1.95e-5 | 4.88e-6 | 9.5e-11 | 9.5e4 | 0 |
   | 3 | 16384 | 4.88e-6 | 1.22e-6 | 6.0e-12 | 6.0e3 | 0 |

   Even on level 3, the smallest triangle is 77 times the threshold area.
   The interface meshes are uniform (1.56 x 3.13 mm faces on level 3), and the
   largest displacement of about 1e-4 m does not change this. The fluid and
   solid interface meshes are identical (conformal).

**A/B run.** Level 1 at the corrected inflow was run from rest with two
builds that differ only by PR #546, on the same machine, with the same
decomposition. Over t = 0 to 1.8875 s (302 time steps, covering the inflow
ramp and the start-up transient), the displacement history at A (all
components) and the FSI residual history (1976 coupling iterations) are
**bit-identical**. The total force history `F` is also bit-identical. One
entry of one viscous force component differs by 1e-17 N in its last printed
digit. This comes from the force function object, which evaluates the
wall-velocity gradient through `newMovingWallVelocity`, a condition that
PR #546 rewrote. It does not reach the solution.

| t = 1.8875 s | A: aee0c8e35 (before #546) | B: aee0c8e35 + #546 | relative difference |
| --- | ---: | ---: | ---: |
| u_x(A) (m) | 5.08170427929e-05 | 5.08170427929e-05 | 0 |
| F_x (N) | 1.210104799124 | 1.210104799124 | 0 |

The comparison stops at t = 1.89 s, before the steady state at t = 20 s.
Reaching the steady state from rest would take about 8 h per run on 8 cores
of the shared node. Because the two operators give bit-identical solutions,
the A/B difference at the steady state is also zero, so the run was stopped.
The recorded level-1 values (`../beam_richter/r02_sequence.csv`:
u_x(A) = 5.4793363e-05 m, F_x = 1.2843702 N at the steady state) were
therefore not rerun on this platform. Such a rerun would be a cross-platform
replicate (GCC 11.4 vs 13.3), not a test of the fix. No time history at
0.2 m/s is stored in paper-data to compare against.

## Run set-up (exactly what was run)

- Case: `richter_case.py` from solids4foam branch `verification/beam-richter`
  (`7e8cd7379`), the generator of the paper's runs:

  ```
  richter_case.py R02_L1_rest --factor 2 --thickness-multiplier 2 \
      --inlet-peak 0.2 --coupling-tol 1e-5 --n-outer-corr 30 \
      --allow-unconverged --pimple-outer 3 --p-tol 1e-10 --u-tol 1e-10 \
      --motion-tol 1e-12 --end-time 20 --ranks 8
  ```

  This is the production set-up of the 0.2 m/s sequence (`R02_L1_rest`). It
  uses a 64,896-cell fluid and a 4,096-cell solid, dt = 0.00625 s, a 1 s ramp,
  IQN-ILS with the strict test at 1e-5, PIMPLE 3, a Newton solid with SNES
  rtol/stol 1e-6 and directMap transfer. The full metadata is in
  `campaign_case.json`. The only difference is 8 ranks instead of the
  original count. The case was meshed, decomposed (scotch) and run with
  `campaign/run_case.sh` (blockMesh, decomposePar, `mpirun -np 8 --bind-to
  none solids4Foam -parallel`).
- Builds (xenosim, Slurm jobs 21856 and 21857; runs 21860 and 21861).
  Isolated in `~/OpenFOAM/s4f-issues/beamab/{A,B}`, with OpenFOAM v2412
  (system), GCC 11.4 and PETSc 3.24. A private copy of `wmake` sets
  `c++OPT = -O2 -fno-tree-vectorize -march=znver2`, the MeluXina flags
  without `-fno-math-errno`.
  - A: `git archive aee0c8e35`, the paper's build commit.
  - B: `aee0c8e35` plus `git diff 978ffc8dd^1 978ffc8dd` (the whole of PR #546,
    applied cleanly). This isolates the fix from the 100+ later development
    commits.
- Evidence: `A_before546/`, `B_with546/` (`solidPointDisplacement_displacement.dat`,
  `force.dat`, `fsiResiduals.dat.gz`). Compare with `python3
  scripts/ab_compare.py A_before546 B_with546`.
