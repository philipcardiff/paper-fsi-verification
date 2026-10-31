# Slurm accounting and job ledger

Source: `sacct -X` for the 27 jobs of this campaign (`../provenance/slurm_jobs.csv`), re-queried when this package was prepared; the ledger rows (`../provenance/ledger_extract.csv`) were compared job by job. CPUTimeRAW = allocated cores x elapsed wall seconds (verified for every row). Jobs 14185 and 14194 belong to another agent and are not counted. All jobs ran the stock binary (sha256 f5ae1ed6f347f875c392f93d847f2369f437345a9b785296d46b0bd6fbb72f02, recorded by every run in `run_meta_start.json`; see `../provenance/run_meta_sanitized.csv`) with 2 threads (4 for the RT-C M1 job) on the shared partition; no run used the heavy lock.

**Total: 27 jobs, 486 allocated CPU-seconds = 0.1350 CPU-h** (cap 3.0 CPU-h). The sacct total, the ledger job set (identical ids) and the campaign results document agree; the campaign document's own job table lists 25 jobs / 0.1306 CPU-h and adds jobs 14212 (4 CPU-s) and 14213 (12 CPU-s) in later sections.

| Job | Case | Role | sacct state | Cores | Elapsed s | CPU-s |
|---|---|---|---|---|---|---|
| 14186 | RTA_law1 | law test, design error (not evidence) | COMPLETED | 2 | 1 | 2 |
| 14187 | RTA_rotcheck | RT-A rotation expression | COMPLETED | 2 | 1 | 2 |
| 14188 | RTA_law1b | law test, design error (not evidence) | COMPLETED | 2 | 1 | 2 |
| 14189 | RTA_law1c | RT-A law identity | COMPLETED | 2 | 1 | 2 |
| 14190 | RTA_twist_M1 | RT-A twist M1, input rejected | FAILED | 2 | 1 | 2 |
| 14191 | RTA_twist_M1 | RT-A twist M1, input rejected | FAILED | 2 | 1 | 2 |
| 14192 | RTA_twist_M1 | RT-A twist M1, outer reaction not logged | COMPLETED | 2 | 2 | 4 |
| 14193 | RTA_twist_M1 | RT-A twist M1 (final) | COMPLETED | 2 | 1 | 2 |
| 14195 | RTA_twist_M2 | RT-A twist M2 | COMPLETED | 2 | 4 | 8 |
| 14196 | RTA_twist_M3 | RT-A twist M3 | COMPLETED | 2 | 19 | 38 |
| 14197 | SMK_fluid_M0 | smoke fluid, input rejected | FAILED | 2 | 0 | 0 |
| 14198 | SMK_fsi_M0 | smoke coupled MR, strict min_residual | FAILED | 2 | 1 | 2 |
| 14199 | SMK_fluid_M0 | smoke fluid | COMPLETED | 2 | 0 | 0 |
| 14200 | SMK_fsi_M0 | smoke coupled MR, stall | FAILED | 2 | 1 | 2 |
| 14201 | SMK_fsi_M0 | smoke coupled MR E_mesh 100, stall | FAILED | 2 | 1 | 2 |
| 14202 | SMK_fsi_M0_svk | diagnostic SVK, converges | COMPLETED | 2 | 2 | 4 |
| 14203 | SMK_fsi_M0_tf0 | diagnostic three-field tag, input rejected | FAILED | 2 | 1 | 2 |
| 14204 | SMK_fsi_M0_tf0 | diagnostic three-field off, stall | FAILED | 2 | 1 | 2 |
| 14205 | SMK_fsi_M0_cm | diagnostic outer mesh clamp, stall | FAILED | 2 | 1 | 2 |
| 14206 | SMK_ale_M0 | diagnostic mixture only, input rejected | FAILED | 2 | 1 | 2 |
| 14207 | SMK_ale_M0 | diagnostic mixture only, converges | COMPLETED | 2 | 1 | 2 |
| 14208 | SMK_fsi_M0_mu0p3 | diagnostic viscosity 0.3, stall | FAILED | 2 | 1 | 2 |
| 14209 | RTB_fluid_M1 | RT-B M1 (as run) | COMPLETED | 2 | 29 | 58 |
| 14210 | RTB_fluid_M2 | RT-B M2 (as run) | COMPLETED | 2 | 134 | 268 |
| 14211 | RTB_fluid_M1_bcshift | RT-B M1 dt/2 diagnostic | COMPLETED | 2 | 29 | 58 |
| 14212 | SMK_fsi_M0_rhoi0p9 | stall diagnostic rho_inf 0.9 | COMPLETED | 2 | 2 | 4 |
| 14213 | RTC_fsi_M1 | RT-C M1, diverged | FAILED | 4 | 3 | 12 |
| total | 27 jobs | | | | | **486** (= 0.1350 CPU-h) |

By stage: RT-A static 10 jobs 0.0178 CPU-h; smoke and coupled diagnostics 13 jobs (14197-14208 and 14212) 0.0072 CPU-h; RT-B 3 jobs 0.1067 CPU-h; RT-C M1 1 job 0.0033 CPU-h; sum 0.1350 CPU-h.

## Ledger versus sacct
The ledger status column of jobs 14190, 14191, 14197 and 14206 reads "completed-normal-termination" although sacct reports FAILED (exit code 1: FEBio rejected the input; duplicate node-set name, unrecognised tag, wrong material id). The ledger's validation text describes the failure correctly. sacct is authoritative; the extract adds the columns `sacct_state` and `ledger_status_note`. See `../provenance/INCONSISTENCIES.md`.

Cost basis (measured): fixed-geometry fluid 2.8e-4 core-s per cell-step at M1 and 3.3e-4 at M2, three to four times below the plan's lowest rate (1.9 Newton iterations per step in a nearly linear problem). A fixed-geometry traction study at M3 (16,384 cells, 200 steps) would cost about 0.3 CPU-h at that rate (an estimate, not run).

## Not in sacct
The campaign record notes that one FEBio process was started on the login node by mistake for under a second while diagnosing a failed XML parse; it printed nothing, wrote nothing and is not in sacct (its cost is not counted in the 0.1350 CPU-h). No other login-node run was made. The package preparation itself ran no simulation (analysis scripts only).
