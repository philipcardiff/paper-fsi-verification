# Corrections to this package (post-submission audit)

Wording corrections only; no data, script or number was changed.

1. `docs/STAGE_B_STRUCTURAL.md`: the shear-locking attribution for hex8 was labelled "PROVEN by comparison". It is an inference from mesh and element-type comparison, not a direct locking test; relabelled SUPPORTED.
2. `README.md`: "four hex8 levels converge monotonically with observed order 1.84" now states that the order is computed from the three finest levels (nt 4-8-16). The asymptotic range is therefore not demonstrated, and the Richardson offset beyond the reference is an extrapolation.
3. `docs/STAGE_C_CFD1.md`: removed the sentence that a second-order refinement "could not cure" the 4 Pa odd-even amplitude. Level 2 was never run; refinement behaviour is not tested (consistent with `docs/KNOWN_LIMITATIONS.md`).

Open items noted by the audit, not changed here: the "at least 8 elements for 1 %" statement for hex8 in the structural document is not supported by the document's own table (nt 8 gives about -1.5 % in u_y and about -3 % in u_x); `docs/SCIENTIFIC_SUMMARY.md` states that coupled stage D is "not recommended for the stock formulation", which rests on one mesh level, K = 1e6 and dt = 0.5 s.
