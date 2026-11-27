#!/bin/bash
# Robin-Neumann, flat MMS variant, full interface pressure condition (MMS tree + #777, env-mmsmwp.sh), sequential serial runs
P=/Users/philipc/claude-runs/mms-partitioned
ENVF=/Users/philipc/claude-runs/env-mmsmwp.sh
cd /Users/philipc/claude-runs/quick_reruns/mms_flat_robin_full
dirs=""
for n in 16 32 64; do
  dt=$(python3 -c "print(0.32/$n)"); d=runs/n$n; dirs="$dirs $d"; rm -rf $d
  python3 $P/pkg/case/makeCase.py $d --coupling robin --variant flat --nx $n --dt $dt --end 0.4 --writeInterval 0.1 --outerTol 1e-6 --np 1 --env $ENVF \
     --pimpleExtra "fluxConsistentPatches (outlet);" --nOuter 2 --ddtCorr no --solidSurfGrad --exactOldTime > /dev/null
  echo "robin --variant flat --pimpleExtra fluxConsistentPatches (outlet); --nOuter 2 --ddtCorr no --solidSurfGrad --exactOldTime env=$ENVF" > $d/ARGS
  /usr/bin/time -p $d/Allrun 2> $d/time.log
done
python3 $P/pkg/case/orders.py $dirs --times 0.2 0.4 --cols U_L2 U_Linf p_L2 p_Linf D_L2 D_Linf Us_L2 intD_L2 intTrac_L2 intTrac_Linf intP_L2 > orders.md
cat orders.md
