#!/bin/bash
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --constraint=infiniband
#SBATCH --time=12:00:00
#SBATCH --no-requeue
#SBATCH --job-name=rntolfull
#SBATCH --output=/scratch/pcardiff/mwp721/ring/quick/%x-%j.out
# Ring strong level, mesh 2, 100 steps per period, standard solid, Robin-Neumann strict test,
# tolerance 1e-4 to 1e-8, Robin pressure/flux tolerances 10x/5000x, full elasticWallPressure (PR #777 head b9b5e2b9e, viscousTerm default yes)
cd /scratch/pcardiff/mwp721/ring
for t in 1e-4 1e-5 1e-6 1e-7 1e-8; do
  pt=$(python3 -c "print('%g' % (10*$t))"); ft=$(python3 -c "print('%g' % (5000*$t))")
  PRE="sed -i -e 's/^\( *robinPressureTolerance\).*/\1 $pt;/' -e 's/^\( *robinFluxTolerance\).*/\1 $ft;/' -e 's/^\( *writeResidualsToFile yes;\)/\1\n    requireAllResidualMeasures yes;/' constant/fsiProperties.robin" \
    scripts/run.sh tol_rn_strong_m2_n100_tol$t strong 2 100 4 robin standard tol=$t > cases/tol_rn$t.out 2>&1 &
done
wait; date
