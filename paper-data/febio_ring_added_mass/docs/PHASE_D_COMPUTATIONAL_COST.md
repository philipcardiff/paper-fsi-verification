# PHASE D COMPUTATIONAL COST (FEBio v4.13, 1 core each, shared node xenosim)
Source: `results/phaseD_febio_summary.csv` (time.txt via `/usr/bin/time -v`). All FEBio runs: Slurm `main`, 1 core, MKL Pardiso (Phase A measured poor thread scaling, hence 1 core; no new scaling test was run), 400 steps.

| case | nodes | wall | CPU-h | peak RSS |
|---|---|---|---|---|
| dry M1 / M2 / M3 | - | 12 s / 48 s / 244 s | 0.003 / 0.012 / 0.057 | 41 / 56 / 119 MB |
| weak M1 / M2 | 1666 / 6402 | 293 s / 1756 s | 0.081 / 0.466 | 118 / 342 MB |
| moderate M1 / M2 / M3 | 1666 / 6402 / 25090 | 273 s / 1597 s / 11538 s (3 h 12 min) | 0.076 / 0.444 / 3.19 | 118 / 357 / 1405 MB |
| strong M1 / M2 | 1666 / 6402 | 298 s / 2033 s | 0.083 / 0.543 | 121 / 346 MB |
| D6 (failed at step 75) | 6402 | partial | ~0.1 | - |

Wet cost grows ~5.8x (M1->M2) and ~7.2x (M2->M3); M3 was slower than the Phase C M3 run (2 h 35 min) because the node was more heavily loaded. Total new Phase D FEBio CPU time ≈ 4.9 core-hours (wet) + dry reuse (none new).
Caveat: solids4foam times in the README (mesh 2, 126-576 s) come from a different machine, serial OpenFOAM, a partitioned solver and a high-order solid; solids4foam memory was not recorded. **No claim is made that either solver is faster or cheaper**: hardware, discretisation, accuracy, time integration and core count differ. Figure: `figures/fig_D4_cost.png`.
