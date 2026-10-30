#!/bin/bash
# Reproduce everything that can be reproduced from THIS PACKAGE ALONE (python3 + numpy + scipy + matplotlib; no FEBio, no HPC output).  Run from the package root.
# Steps that need the retained HPC outputs or a FEBio rerun are listed in docs/REPRODUCIBILITY_AUDIT.md and are NOT run here.
set -e
O=${1:-$(mktemp -d)}; mkdir -p "$O"; echo "scratch output directory for recomputed JSON: $O"
echo "== 1. regenerate every FEBio input from the generators; compare sha256 (expected: 2 'DIFFERENT' = early failed attempts, monitor node sets only)"
python3 -I scripts/check_inputs.py
echo "== 2. verify headline numbers against the raw stored outputs (68 checks)"
python3 -I scripts/verify_against_raw.py | tail -3
echo "== 3. rebuild the summary tables and figures"
python3 -I scripts/build_tables.py > /dev/null && python3 -I scripts/make_figures.py
echo "== 4. force reconstruction F1..F4 from the stored facet tables"
python3 -I scripts/reconstruct_forces_from_facets.py $O/force_reconstruction.csv | tail -5
echo "== 5. package-only analyses (operator spectra, stabilisation matrix checks, plate response to the pair-mean load); compare with results/raw copies"
python3 -I scripts/b4_operator.py $O/b4op > /dev/null
python3 -I scripts/b4_modal_scaling.py $O/b4modal.json 1 | tail -2
python3 -I scripts/b4_straight_channel.py $O/b4sc.json > /dev/null
python3 -I scripts/b4_transient_stokes.py $O/b4ts.json > /dev/null
python3 -I scripts/b7_ppp_check.py $O/b7ppp.json > /dev/null
python3 -I scripts/b7_stabilisation_check.py $O/b7stab.json > /dev/null
python3 -I scripts/run_plate_q9_pairmean.py results/raw/b2/ybranch_M3.json 20 2 $O/q9pm.json > /dev/null
python3 -I scripts/poiseuille_check.py pre | tail -3
python3 -I - <<PY
import json
def same(a, b):
    if isinstance(a, dict): return all(same(a[k], b[k]) for k in a if k in b)
    if isinstance(a, list): return all(same(x, y) for x, y in zip(a, b))
    if isinstance(a, (int, float)) and isinstance(b, (int, float)): return abs(a - b) <= 1e-8 * max(1, abs(b))
    return True
for new, old in [("b7ppp.json", "pressure_mode/b7_ppp_check.json"), ("b7stab.json", "pressure_mode/b7_stabilisation_check.json"), ("b4sc.json", "pressure_mode/b4_straight_channel.json"), ("b4ts.json", "pressure_mode/b4_transient_stokes.json"), ("q9pm.json", "b2/plate_q9_pairmean_M3.json"), ("b4op/b4_operator_spectra.json", "pressure_mode/b4_operator_spectra.json")]:
    print(("identical to stored " if same(json.load(open("$O/" + new)), json.load(open("results/raw/" + old))) else "DIFFERS from stored ") + old)
PY
echo "done"
