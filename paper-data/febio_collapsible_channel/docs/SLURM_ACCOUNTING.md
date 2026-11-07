# Slurm accounting (jobs 14214-14219; allocated CPU time)

Source: `sacct` queried on 2026-10-31 for the six jobs (job-level rows; steps `.batch` and `.extern` omitted); table `provenance/slurm_jobs.csv`. Allocated CPU-seconds = CPUTimeRAW = elapsed x allocated CPUs.

| job | case | cores | elapsed | CPU-s | CPU-h | state |
|---|---|---|---|---|---|---|
| 14214 | S1_80x8 | 2 | 2 s | 4 | 0.0011 | COMPLETED 0:0 |
| 14215 | S1_160x4 | 2 | 2 s | 4 | 0.0011 | COMPLETED 0:0 |
| 14216 | S1_160x8 | 2 | 4 s | 8 | 0.0022 | COMPLETED 0:0 |
| 14217 | S1_320x8 | 2 | 5 s | 10 | 0.0028 | COMPLETED 0:0 |
| 14218 | F1_unpinned | 4 | 31 s | 124 | 0.0344 | COMPLETED 0:0 |
| 14219 | D1_fsi_20steps | 4 | 2 m 18 s (138 s) | 552 | 0.1533 | COMPLETED 0:0 |
| total | | | | **702** | **0.1950** | cap 4.0 CPU-h |

Confirmed: 702 CPU-s = 0.1950 CPU-h, equal to the 0.195 allocated CPU-h of the pilot report. The shared job ledger (`provenance/ledger_extract.csv`, benchmark `collapsibleChannel`, six rows with job id and cores 2/2/2/2/4/4) has no CPU-time column, so the CPU-h figure rests on sacct. No job failed or was cancelled; no resubmission; no job outside this list belongs to the pilot. The four S1 jobs started within the same second (03:42:36 UTC) on a shared, heavily loaded node (loadavg about 105), so elapsed times are short-run wall times, not scaling measurements. The report's statement that the campaign tool charged 2.930 CPU-h in total across all agents at the end was not re-checked (it concerns other agents' jobs).
