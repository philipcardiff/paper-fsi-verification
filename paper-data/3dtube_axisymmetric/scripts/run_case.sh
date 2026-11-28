#!/bin/bash
#SBATCH -N 1
#SBATCH --exclude=sonicgpu9
# Axisymmetric 3dTube: generate and run one case on Sonic.
# Submit with, e.g.:
#   sbatch -p shared -n 4 -t 2-00:00:00 -J axiL3 \
#     --export=ALL,NAME=L3,LEVEL=3,LAW=hooke,EXTRA="" run_case.sh
# The case goes to /scratch/pcardiff/tube_axi/runs/$NAME.
set -o pipefail
ROOT=/scratch/pcardiff/tube_axi
source /scratch/pcardiff/env-s4f671bso.sh
cd $ROOT/runs
CORES=${SLURM_NTASKS:-1}
python3 $ROOT/make_case.py --template $ROOT/ref/tutorials/fluidSolidInteraction/3dTube \
    --level $LEVEL --law ${LAW:-hooke} --cores $CORES --out $NAME ${EXTRA:-} \
    || exit 1
cd $NAME
{ date; echo "job ${SLURM_JOB_ID:-none}"; hostname; lscpu | grep -E "Model name|^CPU\(s\)";
  echo "solids4foam $(git -C /scratch/pcardiff/OpenFOAM/pcardiff-v2512/s4f-671bso log -1 --oneline)";
  which solids4Foam; echo "make_case args: --level $LEVEL --law ${LAW:-hooke} --cores $CORES ${EXTRA:-}";
  echo "load at start: $(uptime)"; } > log.environment 2>&1
echo "Started: $(date)" > log.run
./Allrun >> log.run 2>&1
echo "Finished: $(date) (exit $?)" >> log.run
grep -q "^End" log.solids4Foam && python3 $ROOT/post_axi.py . --csv qoi.csv >> log.run 2>&1
