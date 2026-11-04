# CFD1 rigid-flag refinement (added 2026-11-01): does the nodal-pressure wall force converge?

Adds to the stage-C evidence (level 1 only: F1 drag +15.73 %, lift +157 %). Same model (FEBio `fluid` module, hex8 slab, rho 1000, mu 1 Pa s, K 1e6, 80 steps to 40 s, rigid cylinder plus flag, CFD1 Re 20), pristine v4.13, generator `scripts/make_hronturek_cfd_feb.py` unchanged. Eight additional meshes (factors 1.5 to 3.0; 12.0k to 48.0k cells) plus pin-sensitivity and unpinned runs. Reference: Featflow level 6+0, drag 14.2929, lift 1.11905 N/m (equality with the 2006 paper tables unconfirmed). **FSI1 is still not verified: no coupled run was made.**

## Result (preregistered gate B, `docs/O2_PREREGISTRATION_B.md`: FAILED)
Force errors (percent of the reference; per metre depth). F1 = integral of nodal pressure on the wall plus viscous traction (what `fluid-FSI traction` applies, reconstructed from FEBio's formula); F2 = element-mean reconstruction; F4 = variationally consistent nodal reaction (independent implementation):
| factor | cells | odd-count blocks | F1 drag / lift error % | F2 drag / lift error % | **F4 (consistent reaction) drag / lift error %** | cylinder odd-even a_t (Pa) |
|---|---|---|---|---|---|---|
| 1.25 (stage C, pins by rule) | 8330 | 19 | +15.73 / +157 | -2.12 / -1.12 | **+0.31 / -0.61** | 14.9 |
| 1.5 (pins 0) | 11996 | 6 | +23.70 / -71 | -1.64 / +0.15 | **+0.27 / +0.59** | 65.6 |
| 1.625 (pins 0) | 14191 | 29 | +12.50 / +661 | -1.71 / +0.01 | **+0.23 / +0.25** | 19.0 |
| 1.875 (pins by rule) | 18606 | 18 | -0.83 / +14 | -1.62 / +0.06 | **+0.16 / +0.29** | 6.5 |
| 2.0 (pins 0) | 21344 | 0 | +18.85 / -1 | -1.75 / +0.07 | **+0.21 / +0.37** | 29.9 |
| 2.375 (pins 0) | 29951 | 29 | +12.72 / +689 | -1.28 / +0.06 | **+0.13 / +0.21** | 16.5 |
| 2.5 (pins by rule) | 33368 | 6 | +22.89 / -27 | -1.19 / +0.08 | **+0.12 / +0.28** | 69.1 |
| 3.0 (pins 0) | 48024 | 26 | +3.68 / -292 | -1.27 / +0.14 | **+0.08 / +0.36** | 12.2 |

- Preregistered B1 (monotone decrease of the F1 errors), B2 (L2 error at most 0.7 of L1) and B3 (odd-even amplitude falls by 1.5) **failed**; B4 (soundness: normal termination, mass imbalance at most 7e-4, steady to 1e-11) met. F1 at factor 1.875 (-0.8 %) is an accident of the mode amplitude, not convergence.
- Pins are not the explanation: at factor 1.875 the F1 drag error is -0.83 % (rule pins), -0.77 % (0), -0.88 % (2x rule) and +0.02 % (level-1 values); the pinned and unpinned factor-2.5 runs agree to 0.1 pp (`docs/O2_PREREGISTRATION_B2.md`; `results/b_gate_table.md`).
- Parity test (exploratory, `docs/O2_PREREGISTRATION_B3.md`): the cylinder odd-even amplitude is highest on the two meshes with 6 odd block counts (66-69 Pa) and 6-30 Pa elsewhere (P1, P3 met; P2, P4 failed): parity correlates with the amplitude but does not control the force error.
- **F4 and F2 converge, F1 does not.** F4 is within +0.08 to +0.31 % in drag (+0.2 to +0.6 % in lift) on every mesh, falling with refinement; F2 sits at -1.2 to -2.1 % (O(h) bias).
- Mechanism (derivation, `docs/B_GATE_RESULTS.md`): F1 and F4 differ exactly by the strong-form momentum residual at the constrained wall nodes, which FEBio's discretisation does not enforce there; the nodal pressure has alternating modes that the strong gradient integral does not see but the boundary integral does (Q1/Q1 operator audit, `results/q1q1_parity_operator.json`). Established: the discrepancy is not a resolution, pin or parity effect. Not established: FEBio's coupled applied load on CFD1 (it was reconstructed, not measured), and that the formulation is wrong (the underlying pressure converges).

## Reproduce (no simulation)
`inputs/*.feb.gz` are the inputs; `results/<run>/Mon_all_final_step.txt.gz` holds the final-step nodal data (wx wy ef). To recompute F1, F2, F3b and F4 for a run: gunzip the input and the nodal file into a directory as `<tag>.feb` and `Mon_all.txt` and run `python3 ../scripts/hronturek_cfd1_independent.py <dir> out.json` (checked: it reproduces the stored `independent_analysis.json` exactly for the 1.5 mesh). `scripts/b_gate_table.py` builds the table from the analysis JSONs; the per-step steadiness and mass-balance fields in `stageC_analysis.json` came from the full logs (not stored; 40-220 MB per run). Pristine binary sha256 f5ae1ed6f347f875c392f93d847f2369f437345a9b785296d46b0bd6fbb72f02.
