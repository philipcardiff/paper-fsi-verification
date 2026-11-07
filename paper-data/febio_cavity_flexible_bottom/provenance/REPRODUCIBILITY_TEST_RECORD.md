# Reproducibility test record (scratch copy of the package, one login-node core, python 3.10, numpy 1.21, scipy 1.8, matplotlib 3.5)
| script | result | time |
|---|---|---|
| `scripts/check_inputs.py` | PASS: 18/20 inputs identical; 2 expected differences (INCONSISTENCIES I-9); shipped files match recorded hashes | 2 s |
| `scripts/verify_against_raw.py` | PASS 68/68 | 2 s |
| `scripts/build_tables.py` | PASS (all tables written) | 1 s |
| `scripts/make_figures.py` | PASS | 2 s |
| `scripts/reconstruct_forces_from_facets.py` | PASS: equals the stored b3 table to <=4e-7 | <1 s |
| `scripts/b4_operator.py` | PASS: identical to stored `b4_operator_spectra.json` | 38 s |
| `scripts/b4_modal_scaling.py` | part (1) PASS; part (2) needs C1b `Mon_all.txt` (HPC), skipped with a message | 22 s |
| `scripts/b4_straight_channel.py`, `b4_straight_profile.py`, `b4_transient_stokes.py` | PASS: identical to stored JSON | 13 s, 92 s, 50 s |
| `scripts/b7_ppp_check.py`, `b7_stabilisation_check.py` | PASS: identical to stored JSON | 27 s, 16 s |
| `scripts/run_plate_q9_pairmean.py` | PASS: identical to stored `plate_q9_pairmean_M3.json` (input ybranch_M3, t=500) | 22 s |
| `scripts/poiseuille_check.py pre` | PASS (analytic pre-accounting) | <1 s |
| `reproduce.sh` end to end | PASS | 3 min 12 s |
| `scripts/b2_profile_compare.py` | NOT RUN: needs a reference p(x) table that does not exist | - |
| scripts needing HPC outputs (`b3_force_reconstruction.py` etc.) | PASS in the original working tree on the retained outputs (facet tables regenerated bit-identically: `cmp`); NOT runnable from the package alone | - |
| round-1 and round-2 diffs | apply cleanly (`patch -p2`) to a scratch export of the v4.13 tag | - |
