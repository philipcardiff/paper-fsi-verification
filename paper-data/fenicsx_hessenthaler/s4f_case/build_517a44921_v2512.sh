#!/bin/bash
#SBATCH -J hess-s4fbuild
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=24G
#SBATCH -t 03:00:00
#SBATCH --export=NONE
#SBATCH -o build-%j.out
set -eo pipefail
export USER="${USER:-$(id -un)}"; export HOME="${HOME:-$(getent passwd "$USER" | cut -d: -f6)}"
SHA=517a44921bc9b43ec532e4b30089ef09364e64f6
R=${S4F_ROOT:?set S4F_ROOT}/build
mkdir -p $R/logs $R/src-tree
git -C ${S4F_ROOT:?set S4F_ROOT}/src.git archive $SHA | tar -x -C $R/src-tree
cat > $R/env.sh <<EOT
export PETSC_DIR=/opt/petsc/3.24
export PETSC_ARCH=
export LD_LIBRARY_PATH="/opt/petsc/3.24/lib\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
set +u
. /usr/lib/openfoam/openfoam2512/etc/bashrc WM_PROJECT_USER_DIR=$R
export PATH="$R/src-tree/applications/scripts:\$PATH"
EOT
set +eu; . $R/env.sh; set -eu
export WM_NCOMPPROCS=12
cd $R/src-tree/ThirdParty && ./Allwmake > $R/logs/thirdparty.log 2>&1 || true
cd $R/src-tree/src && ./Allwmake > $R/logs/src.log 2>&1
cd $R/src-tree/applications/solvers && ./Allwmake > $R/logs/apps.log 2>&1
ls -la $FOAM_USER_LIBBIN $FOAM_USER_APPBIN
cat > $R/BUILD_INFO.txt <<EOT
SHA=$SHA (git archive, unmodified)
PETSc=/opt/petsc/3.24; OpenFOAM=v2512 (/usr/lib/openfoam/openfoam2512); $(gcc --version | head -1)
EOT
echo BUILD_OK
