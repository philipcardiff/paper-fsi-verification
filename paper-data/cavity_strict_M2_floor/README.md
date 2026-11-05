# cavityFlexibleBottom: strict partitioned level-2 floor dump (for PR #27)

- **Case:** the solids4foam tutorial at `e8230b59f` (origin/development), run with `./verification/Allverify --variants aitken --levels 1,2`.
- **Coupling:** Aitken with `requireAllResidualMeasures yes`, nOuterCorr 20.
- **Time:** Δt = 0.25 s, end time 200 s.
- **Mesh:** level 2, Δx = 0.05 m, 6,752 fluid and 640 solid cells.
- **Build:** macOS, Apple Clang, v2512, user library built 2026-09-25.
- **Reproduction:** F_y = −5.238953 N/m against the recorded strict value −5.23881 (`5f6f7a6da`), a 2.7e-5 relative difference.

## Files

- **`cavity_strict_M2_floor_t200.dat`:**
  - per-face p and U on the fluid interface (the plate top) at t = 200 s;
  - face-area vectors computed from the written deformed points.
  - Summing p·Sf over these faces reproduces the pressure force of the `forces` function object to 1e-6.
- **`force_M2_strict.dat`:** the `forces` history, with total, pressure and viscous parts.
- **`floor_dump.py`:** the extraction script.

The per-face viscous traction is not included. OpenFOAM's `wallShearStress` function object skips the interface patch because it is of type `patch`, not `wall`. Integrated, the viscous part of F_y is −0.00167 N/m, or 0.03% of F_y.
