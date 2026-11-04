#!/usr/bin/env bash
# Re-run every analysis (no simulation) from the stored raw outputs.   usage: run_analyses.sh PACKAGE_ROOT WORK_DIR OUT_DIR
# WORK_DIR gets regenerated inputs (mesh_rtr.npz, model_meta.json) and the decompressed raw logs; OUT_DIR gets the result files.
set -eu
ROOT=$1; W=$2; OUT=$3; PY="python3 -I"
mkdir -p $W $OUT/campaign_analysis
$PY $ROOT/scripts/regen_inputs.py $ROOT $W/gen > $OUT/regen_inputs.txt || true
# run directories as the campaign analysis scripts expect them: regenerated mesh/meta + decompressed text logs
for r in RTA_law1c RTA_rotcheck RTA_twist_M1 RTA_twist_M2 RTA_twist_M3 RTB_fluid_M1 RTB_fluid_M2 RTB_fluid_M1_bcshift; do
  mkdir -p $W/runs/$r; cp $W/gen/$r/* $W/runs/$r/
  for f in $ROOT/results/raw/$r/*.txt; do [ -f "$f" ] && cp $f $W/runs/$r/; done
  for f in $ROOT/results/raw/$r/*.txt.xz; do [ -f "$f" ] && xz -dc $f > $W/runs/$r/$(basename $f .xz); done
done
$PY $ROOT/scripts/rtr_check_rta.py law1 $W/runs/RTA_law1c $OUT/campaign_analysis/RTA_law1c.json > /dev/null
$PY $ROOT/scripts/rtr_check_rta.py rotcheck $W/runs/RTA_rotcheck $OUT/campaign_analysis/RTA_rotcheck.json > /dev/null
for m in 1 2 3; do $PY $ROOT/scripts/rtr_check_rta.py twist $W/runs/RTA_twist_M$m $OUT/campaign_analysis/RTA_twist_M$m.json > /dev/null; done
for r in RTB_fluid_M1 RTB_fluid_M2; do
  $PY $ROOT/scripts/rtr_analyze.py fluid $W/runs/$r $OUT/campaign_analysis/${r}_as_run.json > /dev/null
  $PY $ROOT/scripts/rtr_analyze.py fluid $W/runs/$r $OUT/campaign_analysis/${r}_aligned_dt2.json 0.005 > /dev/null
done
$PY $ROOT/scripts/rtr_analyze.py fluid $W/runs/RTB_fluid_M1_bcshift $OUT/campaign_analysis/RTB_fluid_M1_bcshift_as_run.json > /dev/null
# independent re-derivation (own code); it parses the .feb files
mkdir -p $W/feb; for r in RTB_fluid_M1 RTB_fluid_M2 RTB_fluid_M1_bcshift; do cp $W/gen/$r/$r.feb $W/feb/; done
$PY $ROOT/scripts/independent_check.py $ROOT/results/raw $W/feb $OUT/independent_check.json > $OUT/independent_check.txt
$PY $ROOT/scripts/log_stats.py $ROOT/results/raw $OUT > $OUT/log_stats.txt
