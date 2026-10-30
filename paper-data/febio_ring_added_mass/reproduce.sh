#!/bin/bash
# Reproduce (a) the FEBio inputs from the generator, (b) the tables/figures from the stored per-run outputs.
# Running FEBio itself (step c) needs FEBio v4.13 + MKL; see README.md. Run from the package root.
set -e
echo "== (a) regenerate the final inputs and compare with feb_sha256 recorded in each runs/*/model_meta.json"
T=$(mktemp -d)
gen() { # tag rho mesh_factor dir
  python3 -I tools/make_ring_feb.py --out $T/$1 --tag $1 --rho-f $2 --k 1e8 --e-mesh 0.001 --rhoi 1 --mesh-factor $3 --steps-per-period 100 --periods 4 \
    --ic potential --ic-field compatible --ic-mode pernode --wall slip --interface slip > /dev/null
  a=$(sha256sum $T/$1/$1.feb | cut -d' ' -f1); b=$(python3 -I -c "import json,sys;print(json.load(open('$4/model_meta.json'))['feb_sha256'])")
  [ "$a" = "$b" ] && echo "identical  $1" || echo "DIFFERENT  $1"; }
gen D1_weak_M1 14 1 runs/phaseD/D1_weak_M1;   gen D4_weak_M2 14 2 runs/phaseD/D4_weak_M2
gen D2_strong_M1 1400 1 runs/phaseD/D2_strong_M1; gen D3_strong_M2 1400 2 runs/phaseD/D3_strong_M2
gen D5_moderate_M3 140 4 runs/phaseD/D5_moderate_M3
echo "== (b) rebuild CSVs and figures D1-D5 from stored outputs"
python3 -I tools/phaseD_results.py
echo "== (c) to re-run a case: sbatch ... tools/run_febio_job.sh <run_dir> <tag>.feb  (needs FEBio v4.13 built per README.md), then tools/phaseD_postprocess.sh <wet> <dry> <NT> <NR>"
echo "== checksums"; sha256sum -c SHA256SUMS --quiet 2>&1 | head -5 || true
