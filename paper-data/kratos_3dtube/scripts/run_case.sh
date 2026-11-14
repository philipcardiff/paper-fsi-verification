#!/bin/bash
# Run one Kratos case from an input set.
# Usage: [MAIN=MainKratos_tube.py] run_case.sh <input dir> <run dir> <threads> [extra MainKratos args]
# KRATOS_ENV may name a script that activates the Kratos 10.3.0 Python environment; without it the
# smartsimet environment is used if present, else the python on PATH must already import Kratos.
# Run from the case directory (on smartsimet ~/raid/kratos/cases/tube3d). The run dir must not exist.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
IN=$(cd "$1" && pwd)
RUN=$2
THREADS=$3
shift 3

if [ -e "$RUN" ]; then echo "run dir $RUN exists, refusing to overwrite" >&2; exit 1; fi
mkdir -p "$RUN"
cp "$IN"/* "$RUN"/
MAIN=${MAIN:-MainKratos.py}
cp "$HERE/$MAIN" "$RUN"/
cd "$RUN"
sha256sum ./* > inputs.sha256

KRATOS_ENV=${KRATOS_ENV:-/home/cloner0110/raid/kratos/kratos-env.sh}
if [ -f "$KRATOS_ENV" ]; then source "$KRATOS_ENV"; fi
{
  echo "input_set=$IN"
  echo "threads=$THREADS"
  echo "start=$(date -u +%FT%TZ)"
  echo "load_before=$(cut -d' ' -f1-3 /proc/loadavg)"
} > run_meta.txt

set +e
OMP_NUM_THREADS=$THREADS /usr/bin/time -v python "$MAIN" ProjectParameters.json "$@" > log.txt 2> time.txt
STATUS=$?
set -e

{
  echo "end=$(date -u +%FT%TZ)"
  echo "load_after=$(cut -d' ' -f1-3 /proc/loadavg)"
  echo "exit=$STATUS"
  echo "max_iterations_warnings=$(grep -c 'maximum iterations' log.txt || true)"
} >> run_meta.txt
exit $STATUS
