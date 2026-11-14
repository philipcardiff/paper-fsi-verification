#!/bin/bash
#SBATCH -J hess-s4fmesh
#SBATCH -p main
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 01:00:00
#SBATCH --export=NONE
#SBATCH -o mesh-%j.out
set -eo pipefail
export USER="${USER:-$(id -un)}"; export HOME="${HOME:-$(getent passwd "$USER" | cut -d: -f6)}"
R=${S4F_ROOT:?set S4F_ROOT}/build
set +eu; . $R/env.sh; set -e
export OMP_NUM_THREADS=8
cd ${S4F_ROOT:?set S4F_ROOT}/runs/feF1_strict
ln -sfn meshDict.coarse system/meshDict
(time surfaceFeatureEdges -angle 45 constant/triSurface/fluidDomain.stl.gz constant/triSurface/fluidDomain.fms) > log.surfaceFeatureEdges 2>&1
(time cartesianMesh) > log.cartesianMesh 2>&1
mkdir -p constant/fluid; mv constant/polyMesh constant/fluid/polyMesh
blockMesh -region solid > log.blockMesh 2>&1
checkMesh -region fluid > log.checkMesh.fluid 2>&1 || true
checkMesh -region solid > log.checkMesh.solid 2>&1 || true
echo MESH_DONE
