# FEBio v4.13 cavityFlexibleBottom: independent diagnostic evidence and limitations

**Status: archival diagnostic package (draft). It is not a successful FEBio verification of cavityFlexibleBottom.** It preserves a closed investigation: what was run, what failed, what was learned, and how to reproduce or reopen it. No stock FEBio source change is proposed; no manuscript change is requested.

Evidence labels: PROVEN / STRONGLY SUPPORTED / HYPOTHESIS / FAILED / NOT TESTED. Numerical residual convergence, physical accuracy and benchmark verification are kept apart; agreement of numbers is never called verification.

## Executive summary
- FEBio v4.13 (stock tag, monolithic fluid-FSI, weakly compressible K=1e6, conformal hex8, StVK plate) was run at three mesh levels (dx 0.1/0.05/0.025; 1688/6752/27008 fluid cells). The model matches the benchmark in geometry, material law, fluid parameters, inlet flux, outlet datum, depth normalisation and output definitions, plus two documented differences (four pinned corner pressures; the equal-order hex8 discretisation).
- **All three coupled runs ended in an error termination** (dilatation-increment norm on the arithmetic floor at steady state); **M3 is not steady to the 1e-4 criterion**. Values are conditional estimates.
- u_y(A) = -0.2120 / -0.2262 / -0.2300 m (latest accepted t = 972.5 / 870 / 715 s); F_y (applied wall-row load, pressure part) = -4.718 / -4.934 / -4.983 N/m. Against solids4foam at the same dx the force is **-7.9 / -5.6 / -5.1 %** and the displacement +3.6 / -2.2 / -3.5 %. The gap flattens near 5 % (indicative extrapolation F_y about -5.00 vs -5.26 N/m). The displacement gap is the force gap (the plate responds proportionally).
- Excluded by measurement: compressibility, time step, tolerances, corner pins (<=0.4 %, wrong sign), viscous traction (0.03 %), shared pressure boundary conditions, plate stiffness, pressure transfer (pressure part).
- **Diagnosed numerical feature (STRONGLY SUPPORTED):** the nodal pressure alternates row by row (D = 0.36/0.17/0.13 Pa) - a soft, h^2-scaling parity mode of the unstabilised equal-order Q1/Q1 fluid element (no interior pressure stabilisation exists in v4.13). The load FEBio applies to the plate (F1) samples only the wall row; FEBio's element-mean surface-pressure definition (F2) does not (F1-F2 = 0.72/0.34/0.27 N/m, equal to the parity part). F2 matches solids4foam at M3 only (+6.1 % M1, +0.9 % M2); **F2 is not established as the physical force and this does not explain the gap.**
- **Two experimental stabilisations failed their pre-registered criteria.** Round 1 (alpha 3/4) removed the mode but broke mass conservation (-1.55 %); round 2 (flux-preserving, alpha 1e-3) failed three of six criteria on an exact Poiseuille channel. Neither is a correction or solution; neither is proposed for FEBio. (`experimental/` holds the diffs as REJECTED RESEARCH PROTOTYPES.)
- Not established: the cause of the 5 % force gap; whether the parity-free pressure is the true one; the solids4foam floor pressure/traction (never available).

## Where to read
`docs/FINDINGS_REPORT.md` (full findings), `docs/SOURCE_CODE_INVESTIGATION.md`, `docs/FAILED_STABILISATION_HISTORY.md`, `docs/COMPARISON_TABLES.md`, `docs/TEST_MATRIX.md` (61 tests; `results/TEST_MATRIX.csv`), `docs/LIMITATIONS.md`, `docs/REPRODUCIBILITY_AUDIT.md`, `provenance/INCONSISTENCIES.md` (17 discrepancies between source reports and raw data, all resolved or stated), `docs/INDEX.md` (30 archived working reports).

## Relation to other manuscript-repository work
Manuscript PRs #12 (ANSYS and COMSOL results) and #22 (COMSOL and Tukovic references, Richardson extrapolation) concern the **reference displacements** and the plate constitutive law (a 4-5 % displacement spread consistent with other plate laws, while all reference forces agree to 0.2-0.7 %). The present package concerns a different question: a **numerical force deficit of FEBio's own fluid solution** against the StVK-confirmed solids4foam values. The two should not be conflated. No existing FEBio cavity PR or fork branch was found when this was prepared.

## Layout
`docs/` reports; `results/` tables (CSV/JSON; `results/raw/` stored per-run outputs: summaries, point histories, interface logs and solver logs as `.gz`, derived analyses, facet tables); `figures/` (`fig_*` regenerated from retained data; `cavity_*` original figures from HPC outputs, labels and units checked); `scripts/` (analysis and generator scripts, relative paths); `inputs/` (`.feb` below 3 MB; the 6 MB M3 input is regenerated and hash-checked; `INPUTS.csv`); `experimental/` (rejected prototypes, hashes); `provenance/` (build, Slurm jobs and CPU-hours, ledger extract, HPC output hashes, inconsistencies); `SHA256SUMS`; `reproduce.sh`.
Note: the original figures `figures/cavity_*.png` show FEBio M1-M3 values as steady; read them with the caveat above (error-terminated, conditional estimates; M3 not steady).

## Reproduction (details and limits: `docs/REPRODUCIBILITY_AUDIT.md`)
`bash reproduce.sh` (python3, numpy, scipy, matplotlib; about 4 minutes): regenerates and hash-checks all inputs, re-derives the headline numbers from the raw stored outputs (68 checks), rebuilds the tables and figures and the F1..F4 force reconstruction, and reruns the discrete-operator analyses (identical to the stored results). This is bookkeeping on stored outputs, **not** an independent solver verification. Reproducing any FEBio result needs FEBio v4.13 (inputs shipped) and about 38 allocated CPU-h for the whole campaign (M3 alone 28 CPU-h); the 0.51 GB of HPC outputs are not shipped (names, sizes, sha256 in `provenance/hpc_retained_outputs.csv`). Binary provenance is strongly supported but not proven bit-exactly (`provenance/BUILD_PROVENANCE.md`).

## Minimal criteria for reopening the investigation
1. An independent reference **field**, not scalars: solids4foam (or another StVK code) floor p(x), viscous traction and inlet/outlet flux at M2 (and ideally the fixed-plate force). The pre-registered predictions P1-P5 in `docs/reports/CAVITY_FORCE_DEFICIT_DECISION.md` section 3 decide whether FEBio's pressure level, not only its row parity, differs.
2. FEBio output that does not exist yet: nodal w and e_f for all z=0 nodes, per-facet fluid traction (the built-in surface-force plot) and working reaction forces, at the last two states (one M1 run, about 1 CPU-h).
3. Any remedy (stabilisation, element change) must have its acceptance criteria declared first and pass on the exact Poiseuille channel (nodal flux <= 1e-4, dp/dx <= 0.2 %, alternation 10x lower) before a fixed-plate M2 fluid test and one coupled steady run; the 1.5 % conservation failure of round 1 and the nodal-flux sawtooth of round 2 are the known obstacles.
4. A clean steady termination (for example a justified change of the dilatation convergence tolerance) so results are not conditional on an error exit.
5. A clean rebuild with a binary comparison to close the provenance limitation.

## Size and integrity
See `SHA256SUMS`. Nothing here is a compiled binary or a copy of the FEBio source. Sanitised for privacy: no home paths, user names, host names or e-mail addresses.
