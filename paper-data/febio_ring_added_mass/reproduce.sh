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
genx() { # tag rho mesh spp dir [extra generator args]
  tg=$1; r=$2; m=$3; sp=$4; dd=$5; shift 5
  python3 -I tools/make_ring_feb.py --out $T/$tg --tag $tg ${r:+--rho-f $r} --k 1e8 --e-mesh 0.001 --rhoi 1 --mesh-factor $m --steps-per-period $sp --periods 4 \
    --ic potential --ic-field compatible --ic-mode pernode --wall slip --interface slip "$@" > /dev/null
  a=$(sha256sum $T/$tg/$tg.feb | cut -d' ' -f1); b=$(python3 -I -c "import json;print(json.load(open('$dd/model_meta.json'))['feb_sha256'])")
  [ "$a" = "$b" ] && echo "identical  $tg" || echo "DIFFERENT  $tg"; }
echo "== (a2) new inputs (final campaign): M3 weak/strong, time study, K=1e9, loosened tolerances"
genx R_MESH1_weak_M3 14 4 100 runs/ringNext/R_MESH1_weak_M3; genx R_MESH2_strong_M3 1400 4 100 runs/ringNext/R_MESH2_strong_M3
for sp in 50 25; do genx R_TIME_dry_M2_s$sp "" 2 $sp runs/ringNext/R_TIME_dry_M2_s$sp; for c in weak:14 moderate:140 strong:1400; do genx R_TIME_${c%%:*}_M2_s$sp ${c##*:} 2 $sp runs/ringNext/R_TIME_${c%%:*}_M2_s$sp; done; done
python3 -I tools/make_ring_feb.py --out $T/R1 --tag R1_strong_M1_K1e9 --rho-f 1400 --k 1e9 --e-mesh 0.001 --rhoi 1 --mesh-factor 1 --steps-per-period 100 --periods 4 --ic potential --ic-field compatible --ic-mode pernode --wall slip --interface slip >/dev/null; a=$(sha256sum $T/R1/R1_strong_M1_K1e9.feb | cut -d' ' -f1); b=$(python3 -I -c "import json;print(json.load(open('runs/ringNext/R1_strong_M1_K1e9/model_meta.json'))['feb_sha256'])"); [ "$a" = "$b" ] && echo "identical  R1_strong_M1_K1e9" || echo "DIFFERENT  R1_strong_M1_K1e9"
genx R4_strong_M2_loosetol 1400 2 100 runs/ringNext/R4_strong_M2_loosetol --dtol 1e-4 --vtol 1e-4 --ftol 1e-5 --etol 1e-6 --rtol 1e-8
echo "== (b) rebuild CSVs and figures D1-D5 from stored outputs"
RING_ROOT=$PWD python3 -I tools/ring_build_csvs.py && RING_ROOT=$PWD python3 -I tools/ring_figures.py   # needs scipy (compressibility theory)
python3 -I tools/phaseD_results.py   # reads RING_TEMPORAL_STUDY.csv for the time-order check
echo "== (c) to re-run a case: sbatch ... tools/run_febio_job.sh <run_dir> <tag>.feb  (needs FEBio v4.13 built per README.md), then tools/phaseD_postprocess.sh <wet> <dry> <NT> <NR>"
echo "== checksums"; sha256sum -c SHA256SUMS --quiet 2>&1 | head -5 || true
