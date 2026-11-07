# Corrections to this package (audit of 2026-11-01)

Wording corrections only; no data, script logic or number was changed (the one printed number changed, below, is a rounding fix from the stored six-digit values).
1. `scripts/reconstruct_forces_from_facets.py` header: the extraction script is `scripts/b3_force_reconstruction.py` (the cited `scripts/hpc/` path does not exist in the package).
2. `docs/COMPARISON_TABLES.md` row M2: the u_y error against the partitioned reference is -2.25 % from the six-digit values (other documents print -2.2 %; see `provenance/INCONSISTENCIES.md` I-17).
3. `docs/TEST_MATRIX.md` (C1b A4 line): the 'pair-averaged -0.3 %' upstream-gradient figure is flagged as not recomputable from the retained data; the A4 criterion defect remains a diagnosis, not an independent confirmation. The criterion A4 itself stays FAILED as written.
