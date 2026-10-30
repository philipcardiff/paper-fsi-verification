#!/bin/bash
# usage: submit_ring.sh <run_dir> <tag> [partition] [cpus]   (1 core by default: measured scaling, see PHASE_A)
d=$(realpath "$1"); tag=$2; part=${3:-dev}; cpus=${4:-1}
cd "$d" && sbatch --job-name="ring-$tag" --partition=$part --cpus-per-task=$cpus --mem=8G --time=00:28:00 \
  --output=slurm_%j.out --wrap="$HOME/work/febio-fsi-verification/tools/run_febio_job.sh $d $tag.feb"
