#!/bin/bash
#SBATCH --job-name=cavity-su2-calculix
#SBATCH --partition=short            # short <=12 h | cpu <=7 d | long <=14 d
#SBATCH --ntasks=3                   # = NP_FLUID SU2 MPI ranks + 1 CalculiX (serial)
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00              # set per level from the pilot (~3x projected)
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#
# cavityFlexibleBottom, SU2 (incompressible) + CalculiX (quasi-static NLGEOM) via preCICE (Slurm batch
# script; run from the study directory that holds scripts/ and cases/, or set FSI_ENV to the environment
# script). Same script as ../fsi1_su2_calculix/scripts/run-case.sh; only names and comments differ.
#   Stage 1 (steady): SU2_CFD steady_rigid.cfg on all ntasks -> restart_steady.dat (rigid plate).
#   Stage 2 (coupled): restart from stage 1; SU2 (NP_FLUID ranks) + ccx_preCICE, marched to steady state.
# Submit from the study directory:  sbatch --export=ALL,CASE=cases/L1[,STAGE=steady|coupled|all][,MAX_TIME=..] scripts/run-case.sh
#   CASE      case folder (from scripts/setup_case.py)          (required)
#   STAGE     steady | coupled | all   (default all; 'coupled' reuses an existing stage-1 restart)
#   MAX_TIME  coupled end time override for a pilot (preCICE max-time; CalculiX stops when coupling ends)
#   NP_FLUID  SU2 ranks in stage 2 (default ntasks-1);  GRACE  hang-guard grace in s (default 120)
#   STALL     progress watchdog: stop both if no new coupling window for STALL s (default 1800)
set -euo pipefail
echo "Job $SLURM_JOB_ID ($SLURM_JOB_NAME)  host $(hostname)  CPUs ${SLURM_CPUS_PER_TASK:-1}  tasks ${SLURM_NTASKS:-1}  mem ${SLURM_MEM_PER_NODE:-?} MB  GPU ${CUDA_VISIBLE_DEVICES:-none}  start $(date -Is)"
cd "$SLURM_SUBMIT_DIR"
CASE=$(cd "${CASE:?set CASE=cases/L<n>}" && pwd)
STAGE=${STAGE:-all}
NP_FLUID=${NP_FLUID:-$(( SLURM_NTASKS - 1 ))}
GRACE=${GRACE:-120}
STALL=${STALL:-1800}
[ "$NP_FLUID" -ge 1 ] && [ $(( NP_FLUID + 1 )) -le "$SLURM_NTASKS" ] || { echo "ERROR: NP_FLUID=$NP_FLUID needs ntasks >= NP_FLUID+1 (have $SLURM_NTASKS)"; exit 1; }
echo "case $CASE  stage $STAGE  steady ranks $SLURM_NTASKS  coupled: SU2 ranks $NP_FLUID + CalculiX 1"
cat "$CASE/CASE_INFO.txt"

set +u; source "${FSI_ENV:-$SLURM_SUBMIT_DIR/scripts/env-su2-precice.sh}"; set -u   # see scripts/env-su2-precice.sh
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$CASE"
date +%s > wall_start.txt
STATUS=0

# ---------------- stage 1: steady rigid plate (CFD1) ----------------
if [ "$STAGE" = steady ] || [ "$STAGE" = all ]; then
  rm -f fluid-su2/{restart_steady.dat,history_steady.csv,flow_steady*.vtu,surface_steady*,forces_breakdown_steady.dat}
  t0=$(date +%s)
  set +e
  (cd fluid-su2 && mpirun -np "$SLURM_NTASKS" --bind-to none SU2_CFD steady_rigid.cfg) > steady.log 2>&1
  ST_STEADY=$?
  set -e
  echo "$(( $(date +%s) - t0 ))" > walltime_steady.txt
  echo "stage 1 exit $ST_STEADY, wall $(cat walltime_steady.txt) s, $(grep -c '' fluid-su2/history_steady.csv 2>/dev/null || echo 0) history lines"
  [ "$ST_STEADY" -eq 0 ] && [ -s fluid-su2/restart_steady.dat ] || STATUS=1
fi

# ---------------- stage 2: coupled ----------------
if [ "$STATUS" -eq 0 ] && { [ "$STAGE" = coupled ] || [ "$STAGE" = all ]; }; then
  [ -s fluid-su2/restart_steady.dat ] || { echo "ERROR: no stage-1 restart"; exit 1; }
  if [ -n "${MAX_TIME:-}" ]; then
    sed "s|<max-time value=\"[^\"]*\"|<max-time value=\"$MAX_TIME\"|" precice-config.base.xml > precice-config.xml
  else
    cp precice-config.base.xml precice-config.xml
  fi
  diff precice-config.base.xml precice-config.xml > precice-config.diff || true
  rm -f watchdog.txt
  rm -rf precice-exports precice-run fluid-su2/{restart_flow_*.dat,history.csv,history_0*.csv,flow_*.vtu,surface_flow*,forces_breakdown.dat,example_inlet_profile*,inlet_profile_0*.dat,precice-*} \
         solid-calculix/{plate.frd,plate.dat,plate.cvg,plate.sta,plate.12d,spooles.out,precice-*}
  mkdir -p precice-exports
  # 2nd-order dual time stepping restarts from two levels: both = steady rigid solution
  cp fluid-su2/restart_steady.dat fluid-su2/restart_flow_00000.dat
  cp fluid-su2/restart_steady.dat fluid-su2/restart_flow_00001.dat
  # unsteady SU2 reads the inlet profile with a time-step suffix (inlet_profile_00000.dat)
  cp fluid-su2/inlet_profile.dat fluid-su2/inlet_profile_00000.dat
  t1=$(date +%s)
  set +e
  (cd solid-calculix && exec ccx_preCICE -i plate -precice-participant Solid) > solid.log 2>&1 &
  PID_SOLID=$!
  (cd fluid-su2 && exec mpirun -np "$NP_FLUID" --bind-to none python3 "$(command -v SU2_preCICE_FSI.py)" \
      -f coupled.cfg --parallel) > fluid.log 2>&1 &
  PID_FLUID=$!
  # Guards: (1) hang guard - once one participant exits, the other gets GRACE s, then is stopped;
  # (2) progress watchdog - if the preCICE iterations log (one line per window) has not changed
  #     for STALL s, both are stopped (a participant can abort and HANG without exiting: CalculiX
  #     SPOOLES heap corruption in the last 4x window, job 851, hung 8.5 h).
  ITLOG=solid-calculix/precice-Solid-iterations.log
  running() { local s; s=$(ps -o stat= -p "$1" 2>/dev/null); [ -n "$s" ] && [[ $s != Z* ]]; }
  stop() { for p in "$@"; do running "$p" && kill -TERM "$p" 2>/dev/null; done
           for _ in $(seq 30); do local any=0; for p in "$@"; do running "$p" && any=1; done; [ $any = 0 ] && return; sleep 1; done
           for p in "$@"; do running "$p" && { echo "pid $p ignored SIGTERM: SIGKILL"; kill -KILL "$p" 2>/dev/null; }; done; }
  REASON=""
  while :; do
    if ! running "$PID_SOLID" || ! running "$PID_FLUID"; then
      FIRST=$(running "$PID_SOLID" && echo fluid || echo solid)
      for _ in $(seq "$GRACE"); do { running "$PID_SOLID" || running "$PID_FLUID"; } || break; sleep 1; done
      if running "$PID_SOLID" || running "$PID_FLUID"; then REASON="$FIRST exited, partner still running after ${GRACE}s"; fi
      break
    fi
    ref=$(stat -c %Y "$ITLOG" 2>/dev/null || echo "$t1")
    if [ $(( $(date +%s) - ref )) -gt "$STALL" ]; then REASON="no preCICE progress for ${STALL}s (last window $(tail -1 "$ITLOG" 2>/dev/null | awk '{print $1}'))"; break; fi
    sleep 10
  done
  [ -n "$REASON" ] && { echo "WATCHDOG: $REASON: stopping both participants"; stop "$PID_SOLID" "$PID_FLUID"; }
  wait "$PID_SOLID"; ST_SOLID=$?
  wait "$PID_FLUID"; ST_FLUID=$?
  [ -n "$REASON" ] && echo "$REASON" > watchdog.txt
  set -e
  echo "$(( $(date +%s) - t1 ))" > walltime_coupled.txt
  echo "solid=$ST_SOLID fluid=$ST_FLUID" > exit_status_participants.txt
  STATUS=$(( ST_SOLID != 0 || ST_FLUID != 0 ))
  [ -f "$ITLOG" ] && awk 'NR>1{n++; s+=$3; if($3>m)m=$3; if($3>=50)u++} END{if(n) printf "windows %d  mean coupling its %.2f  max %d  at-cap %d\n", n, s/n, m, u+0}' "$ITLOG"
  echo "stage 2 exit $STATUS (solid $ST_SOLID, fluid $ST_FLUID), wall $(cat walltime_coupled.txt) s"
fi

date +%s > wall_end.txt
echo $STATUS > exit_status.txt
WALL=$(( $(cat wall_end.txt) - $(cat wall_start.txt) )); echo "$WALL" > walltime.txt
echo "Exit $STATUS  end $(date -Is)  wall $WALL s"
exit $STATUS
