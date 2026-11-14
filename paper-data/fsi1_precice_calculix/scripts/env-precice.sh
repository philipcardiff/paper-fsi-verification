# Environment for OpenFOAM v2512 + preCICE 3.4.1 runs (Hron-Turek, issue #36)
# Usage:  source scripts/env-precice.sh
source /opt/openfoam/OpenFOAM-v2512/etc/bashrc
# Use the system FFTW (same 3.3.10) instead of ThirdParty's, so that the
# system PETSc pulled in by preCICE resolves all its FFTW symbols.
export LD_LIBRARY_PATH=$(echo "$LD_LIBRARY_PATH" | tr : '\n' | grep -v '/fftw-' | paste -sd:)
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
