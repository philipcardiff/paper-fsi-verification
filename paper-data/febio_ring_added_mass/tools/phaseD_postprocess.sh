#!/bin/bash
# usage: phaseD_postprocess.sh <wet_run_dir> <matching_dry_run_dir> <NT> <NR>   (NT = 48*f, NR = 12*f)
set -e; R=$1; DRY=$2; NT=$3; NR=$4; here=$(dirname "$(readlink -f "$0")"); tag=$(python3 -c "import json;print(json.load(open('$R/model_meta.json'))['tag'])")
python3 -I $here/run_report.py $R $DRY > /dev/null
W=$(python3 -c "import json;print(json.load(open('$R/report.json'))['omega_fit'])")
python3 -I $here/pressure_modes.py $R $W $NT $NR > /dev/null
python3 -I $here/interface_leakage.py $R > /dev/null
python3 -I $here/convergence_margins.py $R/$tag.log > $R/convergence_margins.txt
echo "$R done: omega_wet=$W"
