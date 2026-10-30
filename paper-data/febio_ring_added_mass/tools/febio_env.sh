# source this file; sets the FEBio 4.13 user-space environment (no system changes)
export FEBIO_ROOT=$HOME/work/febio-fsi-verification
export MKL_VENV=$FEBIO_ROOT/tools/mkl-venv
export LD_LIBRARY_PATH=$MKL_VENV/lib:${LD_LIBRARY_PATH:-}
export FEBIO_BIN=${FEBIO_BIN_OVERRIDE:-$(ls $FEBIO_ROOT/build/febio-4.13/bin/febio4 2>/dev/null)}
