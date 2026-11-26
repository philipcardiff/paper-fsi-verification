#!/bin/bash
#SBATCH --job-name=cavity-rigid
#SBATCH --partition=short
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
# Steady SU2 (rigid plate) only, in CASE/fluid-su2 (the checks in ../checks/rigid_plate_steady.*).
#   sbatch --export=ALL,CASE=cases/L1 scripts/run-steady.sh
set -euo pipefail
echo "Job $SLURM_JOB_ID ($SLURM_JOB_NAME)  host $(hostname)  tasks ${SLURM_NTASKS}  mem ${SLURM_MEM_PER_NODE:-?} MB  start $(date -Is)"
cd "$SLURM_SUBMIT_DIR"; CASE=$(cd "${CASE:?}" && pwd); echo "case $CASE"
set +u; source "${FSI_ENV:-$SLURM_SUBMIT_DIR/scripts/env-su2-precice.sh}"; set -u
export OMP_NUM_THREADS=1
cd "$CASE"; date +%s > wall_start.txt
set +e
(cd fluid-su2 && mpirun -np "$SLURM_NTASKS" --bind-to none SU2_CFD steady_rigid.cfg) > steady.log 2>&1
ST=$?
set -e
date +%s > wall_end.txt; echo $ST > exit_status.txt
WALL=$(( $(cat wall_end.txt) - $(cat wall_start.txt) )); echo $WALL > walltime.txt
echo "Exit $ST  end $(date -Is)  wall $WALL s  history lines $(grep -c '' fluid-su2/history_steady.csv 2>/dev/null || echo 0)"
exit $ST
