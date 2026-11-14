#!/bin/bash
#SBATCH -J hess-s4fapp
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH --export=NONE
#SBATCH -o buildapp-%j.out
set -eo pipefail
export USER="${USER:-$(id -un)}"; export HOME="${HOME:-$(getent passwd "$USER" | cut -d: -f6)}"
R=${S4F_ROOT:?set S4F_ROOT}/build
set +eu; . $R/env.sh; set -e
export WM_NCOMPPROCS=4; export FOAM_MODULE_APPBIN=$FOAM_USER_APPBIN FOAM_MODULE_LIBBIN=$FOAM_USER_LIBBIN
cd $R/src-tree/applications/solvers/solids4Foam && wmake > $R/logs/app.log 2>&1
ls -la $FOAM_MODULE_APPBIN/solids4Foam
printf 'SHA=517a44921bc9b43ec532e4b30089ef09364e64f6 (git archive, unmodified)\nPETSc=/opt/petsc/3.24; OpenFOAM=v2512; %s\n' "$(gcc --version | head -1)" > $R/BUILD_INFO.txt
echo BUILD_OK
