#!/usr/bin/env bash
# run-study.sh WORKDIR : the five FSI1 runs, one after another.
#   Mesh study 1x, 2x, 4x (fluid and solid refined together; 4x fluid on 8 ranks),
#   and solid-only refinement (2x, 4x solid) on the 1x fluid mesh.
H=$(dirname "$(readlink -f "$0")"); W=$(readlink -f "${1:?usage: run-study.sh WORKDIR}")
mkdir -p "$W"
python3 "$H/setup-case.py" 1 "$W/level1-rbf"
python3 "$H/setup-case.py" 1 "$W/fluid1-solid2" 2
python3 "$H/setup-case.py" 1 "$W/fluid1-solid4" 4
python3 "$H/setup-case.py" 2 "$W/level2"
python3 "$H/setup-case.py" 4 "$W/level4"
for spec in "level1-rbf 1" "fluid1-solid2 1" "fluid1-solid4 1" "level2 1" "level4 8"; do
    set -- $spec
    echo "=== $1 ($2 fluid ranks) started $(date)"
    bash "$H/run-case.sh" "$W/$1" $2 2>&1 | tee "$W/$1/log.run"
    echo "=== $1 finished $(date)"
done
python3 "$H/study-table.py" "$W"
