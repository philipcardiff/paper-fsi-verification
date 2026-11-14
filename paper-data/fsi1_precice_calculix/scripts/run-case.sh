#!/usr/bin/env bash
# run-case.sh CASEDIR [NP] : Hron-Turek FSI1, OpenFOAM + CalculiX via preCICE.
#   NP = number of fluid MPI ranks (default 1); CalculiX always runs as one process.
#   1) blockMesh;  2) uncoupled pimpleFoam 0 -> 2 s (inlet ramp, rigid plate);
#   3) coupled pimpleFoam + ccx_preCICE 2 -> 30 s.
source "$(dirname "$(readlink -f "$0")")/env-precice.sh"
set -e
cd "$1"; C=$PWD; NP=${2:-1}
{ date; hostname; echo "OpenFOAM $WM_PROJECT_VERSION"; ccx_preCICE -v 2>&1 | grep -i version;
  dpkg -s libprecice3 | grep Version; echo "nproc $(nproc), fluid ranks $NP + 1 solid";
  echo "load at start: $(uptime)"; env | grep -E '^(OMP|OPENBLAS)_NUM_THREADS'; } > log.environment 2>&1

cd $C/fluid-openfoam
blockMesh > log.blockMesh 2>&1
checkMesh > log.checkMesh 2>&1 || true
grep -E '^ *cells:' log.checkMesh
if [ "$NP" -gt 1 ]; then
    foamDictionary system/decomposeParDict -entry numberOfSubdomains -set $NP > /dev/null
    decomposePar -force > log.decomposePar 2>&1
    run() { mpirun -np $NP --bind-to core pimpleFoam -parallel "$@"; }
    setD() { foamDictionary system/controlDict -entry "$1" -set "$2" > /dev/null; }
else
    run() { pimpleFoam "$@"; }
    setD() { foamDictionary system/controlDict -entry "$1" -set "$2" > /dev/null; }
fi

# stage 1: uncoupled ramp
setD functions/preCICE_Adapter/enabled false
setD endTime 2
t0=$(date +%s)
run > log.pimpleFoam.ramp 2>&1
t1=$(date +%s); echo "ramp 0-2 s: $((t1-t0)) s"

# stage 2: coupled
setD functions/preCICE_Adapter/enabled true
setD endTime 30
rm -rf $C/precice-run
set +e
run > log.pimpleFoam.coupled 2>&1 &
pf=$!
( cd $C/solid-calculix && ccx_preCICE -i plate -precice-participant Solid > log.ccx 2>&1 ) &
pc=$!
wait $pf; rf=$?; wait $pc; rc=$?
t2=$(date +%s)
echo "coupled 2-30 s: $((t2-t1)) s  (fluid rc $rf, solid rc $rc)"
echo "total: $((t2-t0)) s on $((NP+1)) cores ($NP fluid rank(s) + 1 solid process)"
echo "load at end: $(uptime)" >> $C/log.environment
