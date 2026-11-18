#!/bin/bash
#SBATCH -J hess-s4frun3
#SBATCH -p main
#SBATCH -n 16
#SBATCH -c 1
#SBATCH --mem=60G
#SBATCH -t 02:30:00
#SBATCH --export=NONE
#SBATCH -o run-%j.out
set -eo pipefail
export USER="${USER:-$(id -un)}"; export HOME="${HOME:-$(getent passwd "$USER" | cut -d: -f6)}"
R=${S4F_ROOT:?set S4F_ROOT}/build
set +eu; . $R/env.sh; set -e
cd ${S4F_ROOT:?set S4F_ROOT}/runs/feF1_strict
if [ ! -d processor0 ]; then
  decomposePar -region fluid > log.decomposePar.fluid 2>&1
  decomposePar -region solid > log.decomposePar.solid 2>&1
fi
echo START $(date +%s)
mpirun -np 16 --bind-to none solids4Foam -parallel >> log.solids4Foam 2>&1
echo END $(date +%s)
