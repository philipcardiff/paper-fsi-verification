# Documentation index

New consolidated documents (written for this package, numbers re-derived from raw data):
- `FINDINGS_REPORT.md` - full findings.
- `SOURCE_CODE_INVESTIGATION.md` - FEBio v4.13 source analysis.
- `FAILED_STABILISATION_HISTORY.md` - the two rejected research stabilisations.
- `COMPARISON_TABLES.md` - FEBio vs the reference classes (generated).
- `TEST_MATRIX.md` (+ `../results/TEST_MATRIX.csv`) - one row per test, 61 rows.
- `REPRODUCIBILITY_AUDIT.md`, `LIMITATIONS.md`, `../provenance/INCONSISTENCIES.md`, `../provenance/BUILD_PROVENANCE.md`.

Archived working reports (`reports/`, verbatim except sanitisation, with a banner; paths in them refer to the original working tree):
| group | reports |
|---|---|
| benchmark definition | CAVITY_BENCHMARK_DEFINITION, CAVITY_EXACT_BENCHMARK_CONTRACT, CAVITY_BENCHMARK_EQUIVALENCE_REVIEW |
| structural | CAVITY_PLATE_ONLY_RESULTS |
| fluid-only | CAVITY_FLUID_ONLY_RESULTS, CAVITY_INLET_BOUNDARY_CONDITION_AUDIT, CAVITY_POISEUILLE_DIAGNOSTIC |
| coupled | CAVITY_FSI_PILOT_RESULTS, CAVITY_MESH_STUDY (interim + final), CAVITY_MESH_CONVERGENCE, CAVITY_NUMERICAL_SENSITIVITY, CAVITY_CROSS_SOLVER_COMPARISON |
| audits | CAVITY_PRESSURE_FORCE_AUDIT, CAVITY_VALIDATION_GATE_AUDIT, CAVITY_PARITY_MODE_CHECK (superseded in part, see INCONSISTENCIES I-6), CAVITY_REFERENCE_FIELD_INVENTORY, CAVITY_SOLVER_FIELD_COMPARISON |
| decisions | CAVITY_FORCE_DEFICIT_DECISION, CAVITY_NEXT_DECISION, CAVITY_SOURCE_CORRECTION_OPTIONS and CAVITY_MINIMUM_VALIDATION_PLAN (proposals, not executed) |
| source | FEBIO_V413_PRESSURE_SOURCE_AUDIT, FEBIO_FSI_TRACTION_ASSEMBLY_AUDIT, FEBIO_PRESSURE_MODE_ANALYSIS, FEBIO_PRESSURE_PROJECTION_DERIVATION |
| experimental | FEBIO_EXPERIMENTAL_BUILD, FEBIO_ALPHA_ZERO_CONTROL, FEBIO_STABILIZED_FLUID_TEST, FEBIO_FLUX_PRESERVING_PROJECTION, FEBIO_STABILIZATION_FINAL_DECISION |
Not included: campaign status files (ring work and coordination notes are unrelated to this package).
Name mapping: M1 = `C2_fsi_pilot_f1`, M2 = `C3_fsi_M2`, M3 = `C3_fsi_M4` (mesh factor 4).
