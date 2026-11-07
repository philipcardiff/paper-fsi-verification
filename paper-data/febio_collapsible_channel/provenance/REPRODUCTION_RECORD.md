# Reproduction record

## Provenance of what was run
| item | value |
|---|---|
| FEBio | stock tag v4.13 (source commit 32ae206ff4881dfb54f62296cd1558e58ed9fcc6, working tree clean when checked); logs report version 4.13.0.32ae206ff; no plugins, no source modification |
| FEBio binary sha256 | f5ae1ed6f347f875c392f93d847f2369f437345a9b785296d46b0bd6fbb72f02 (re-checked on the build in the working tree; also recorded by every job in its run metadata) |
| solids4foam | `development` commit b200718f4f8390c90f9c7e3ded8b74f3e13941cc (tutorial inputs and verification reference files; see INC-9 in INCONSISTENCIES.md) |
| generators / analysis (original sha256 prefix) | make_cc_feb.py 0e0348346f1a (v1, unmodified); cc_s1_analyze.py 9302d42d2090; cc_f1_analyze.py 73f43ce439c2 (unmodified); cc_fsi_analyze.py 9d472cae5bd8 (one package edit: the reference-CSV directory path, marked in the file) |
| input sha256 (full) | S1_80x8 5c782a2d869c127b0e09e2cedc678642daee7ec64f4d38a169ab499f4cd525ae; S1_160x4 d4d97e824d612c7919d892a8a78a7ff56172c203b2197526eb6770b6503ec1ab; S1_160x8 b049c145c8c335e876e7793f8dec81ce2e0dd9296c44cfec89921614c2e54c39; S1_320x8 1d364a8b441d403321bdf6b7a3e22036756d3fe5eaf8fe70e6be1892a78c9bf3; F1_unpinned 7da9da05826c43a36d7dd178e23619892a0a0527925dee644b4ab5b3113f6558; D1_fsi_20steps 2d57045c0de92e63a6c5bc00e0ac4697b0a989309566c77b9a32262929e317e7 |
| shape file sha256 | ad885ec1d6bb8dea9a1bfde248a6d971f81abda0a4e3b89b2b2dfe1b7dc1ce81 |
| jobs | 14214-14219, see slurm_jobs.csv; 2 cores (S1) and 4 cores (F1, D1); OMP/MKL threads = cores |

## Commands used for the simulations (NOT re-run by this package)
```
python3 scripts/make_cc_feb.py --mode solid --out S1_320x8 --tag S1_320x8 --nsc 320 --ns 8       # also 80x8, 160x4, 160x8
python3 scripts/make_cc_feb.py --mode fluid --out F1_unpinned --tag F1_unpinned --nu 40 --nw 320 --nd 80 --nr 16 --k 1e7 --dt 2.5 --t-end 100 --rhoi 0 --shape inputs/shape_S1_320x8.json
python3 scripts/make_cc_feb.py --mode fsi --out D1_fsi_20steps --tag D1_fsi_20steps --nu 40 --nw 320 --nd 80 --nr 16 --nsc 320 --ns 8 --k 1e7 --dt 0.025 --t-end 0.5 --rhoi 0.9 --e-mesh 1 --plot-stride 5
febio4 -i <case>.feb -silent        # FEBio v4.13 + MKL, one Slurm job per case, 2 cores (S1) / 4 cores (F1, D1)
```
The full 140-step coupled run would be the same fsi command with `--t-end 3.5` (NOT run).

## Non-simulation checks (reproduce.sh, run in a scratch copy of the package on the login node, one core, FULL=1)
Environment: Python 3.10.12; numpy 1.21.5 scipy 1.8.0 matplotlib 3.5.1.
| check | result |
|---|---|
| regenerated inputs vs shipped S1 inputs (4) | 4/4 identical |
| regenerated inputs (6) vs ledger input_sha256 | 6/6 identical (F1 and D1 inputs are not shipped; their regenerated hashes equal the run inputs) |
| frozen-shape file regenerated from S1_320x8 Mon_bottom | identical sha256 |
| original analysis scripts re-run on stored raw outputs (S1 x4, F1 final and history, D1) vs stored results | all identical within relative tolerance 1e-9 (the 'run' path field excluded) |
| CSV tables rebuilt | 7/7 byte-identical |
| figures regenerated | 5/5 produced (inspected visually; shipped copies are the same code output) |
| independent re-computation (reverify.py, own code) | 141 quantities, 0 outside tolerance; corrections to the original report in INCONSISTENCIES.md (INC-3, 4, 5) |
| beam re-derivation (10, 20, 40 elements) | midpoints identical to the stored values to 1e-8 |
| SHA256SUMS | all files OK |
| privacy scan (text and gz) | clean: no /home paths, e-mail addresses, user or host names |
Final line of the run: ALL NON-SIMULATION CHECKS PASSED. Note: the shipped S1 and D1 per-step node logs are gzip of the original text files (sha256 of the originals in retained_hpc_outputs.csv).

## Not reproduced
The FEBio runs themselves (no new simulation was made for this package); the oomph-lib reference runs; the solids4foam runs; the claim about FEBio source behaviour of the fluid-FSI traction (carried from earlier campaign documents).
