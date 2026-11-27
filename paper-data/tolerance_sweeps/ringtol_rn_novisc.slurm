#!/bin/bash
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --constraint=infiniband
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --job-name=rntolnov
#SBATCH --output=/scratch/pcardiff/mwp721/ring/quick/%x-%j.out
# A/B: same as ringtol_rn.slurm at tol 1e-6, but elasticWallPressure viscousTerm no (acceleration-only), same build/platform
cd /scratch/pcardiff/mwp721/ring
t=1e-6
PRE="sed -i -e 's/^\( *robinPressureTolerance\).*/\1 1e-05;/' -e 's/^\( *robinFluxTolerance\).*/\1 0.005;/' -e 's/^\( *writeResidualsToFile yes;\)/\1\n    requireAllResidualMeasures yes;/' constant/fsiProperties.robin && sed -i 's/^\( *type *elasticWallPressure;\)/\1\n        viscousTerm     no;/' 0/fluid/p.robin" \
  scripts/run.sh tol_rn_novisc_strong_m2_n100_tol$t strong 2 100 4 robin standard tol=$t > cases/tol_rn_novisc$t.out 2>&1
date
