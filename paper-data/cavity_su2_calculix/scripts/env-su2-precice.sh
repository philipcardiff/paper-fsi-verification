# Environment for the SU2 + CalculiX + preCICE runs (FSI1, cavityFlexibleBottom) (source it; edit the two paths).
#   SU2 7.5.1 built with the preCICE SU2 adapter's replacement files (CDriver.hpp,
#   python_wrapper_structure.cpp) and su2-7.5.1-cstdint.patch (GCC 13), configured with
#   ./meson.py build -Denable-pywrapper=true -Dwith-mpi=enabled
#   Python venv: pyprecice 3.5.0, mpi4py, numpy (+ setuptools for SU2's bundled Meson on Python 3.12).
SU2_INSTALL=${SU2_INSTALL:-$HOME/software/SU2-7.5.1-install}
VENV=${VENV:-$HOME/software/venv-fsi}
export SU2_RUN=$SU2_INSTALL/bin
source "$VENV/bin/activate"
export PATH=$SU2_RUN:$PATH
export PYTHONPATH=$SU2_RUN${PYTHONPATH:+:$PYTHONPATH}
export OMPI_MCA_hwloc_base_binding_policy=none
# Open MPI 4.1.6: the UCX one-sided component makes SU2_CFD abort in MPI_Win_create (SU2_MPI::Init)
export OMPI_MCA_osc=^ucx
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
