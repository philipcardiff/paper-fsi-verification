#!/bin/bash
# usage (inside a Slurm job): run_febio_job.sh <run_dir> <input.feb> [extra febio args...]
# Runs FEBio 4.13 in <run_dir>, writes run.log, run_meta.json, time.txt. Exit code is FEBio's; it is NOT evidence of validity.
set -u
RUN_DIR=$1; FEB=$2; shift 2
source $HOME/work/febio-fsi-verification/tools/febio_env.sh
NT=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$NT MKL_NUM_THREADS=$NT MKL_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores
cd "$RUN_DIR" || exit 99
{
echo "{"
echo " \"slurm_job_id\": \"${SLURM_JOB_ID:-none}\", \"partition\": \"${SLURM_JOB_PARTITION:-none}\", \"cpus_per_task\": $NT,"
echo " \"host\": \"$(hostname)\", \"start\": \"$(date -Is)\", \"loadavg\": \"$(cut -d' ' -f1-3 /proc/loadavg)\","
echo " \"febio_bin\": \"$FEBIO_BIN\", \"febio_sha256\": \"$(sha256sum $FEBIO_BIN | cut -d' ' -f1)\","
echo " \"input\": \"$FEB\", \"input_sha256\": \"$(sha256sum $FEB | cut -d' ' -f1)\","
echo " \"omp_num_threads\": \"$OMP_NUM_THREADS\", \"mkl_num_threads\": \"$MKL_NUM_THREADS\""
echo "}"
} > run_meta_start.json
/usr/bin/time -v -o time.txt $FEBIO_BIN -i "$FEB" ${FEBIO_SILENT--silent} "$@" < /dev/null > run.log 2>&1
RC=$?
echo "{\"exit_code\": $RC, \"end\": \"$(date -Is)\"}" > run_meta_end.json
exit $RC
