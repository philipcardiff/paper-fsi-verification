# Build and binary provenance (what is known, and what is not)

## Stock binary used for every production run
| item | value |
|---|---|
| source | FEBio tag `v4.13`, `git describe` = `v4.13`, HEAD `32ae206ff4881dfb54f62296cd1558e58ed9fcc6` (shallow, grafted checkout) |
| source-tree cleanliness | `git status --short` empty (0 files) and `git diff` empty, checked when the build was audited and again after all experimental work. The bug-report patch for the v4.13 initial fluid velocity defect (a separate topic) is NOT applied to this tree |
| configuration | CMake Release, `-O3 -DNDEBUG`, `USE_MKL=ON` (Intel MKL Pardiso, linker-script stubs), OpenMP, no HYPRE / SuperLU / FFTW; GNU 11.4.0, CMake 3.22.1 |
| build | Slurm job 13579 (16 cores, 1 min 23 s, 1328 CPU-s allocated), finished 2026-10-29; `libfebiofluid.so` linked 18:23:58 and `febio4` 18:24:02, after the last source modification time (18:21:07). A separate patched build (job 13607, initial-velocity patch) exists and was not used for any cavity run |
| `febio4` sha256 | `f5ae1ed6f347f875c392f93d847f2369f437345a9b785296d46b0bd6fbb72f02` (recorded in every `run_meta_start.json` under `results/raw/`) |
| `libfebiofluid.so` sha256 | `178ba14e1230b5c4713a04c1725371c1a3ac51853dee6be6280c26e7ac48af2b` |
| version string in every log | `4.13.0.32ae206ff` (equals the tag commit; see the first lines of the `.log.gz` files) |
| runtime | `OMP_NUM_THREADS = MKL_NUM_THREADS =` allocated cores (2, 4 or 8), Slurm partition `main` on a shared cluster node that was heavily loaded during the runs (load average about 150), so wall times are not performance data |

### The limitation (not disguised)
That this binary was built from exactly this tree is **STRONGLY SUPPORTED, not PROVEN**: the evidence is the clean tag, the CMake cache source path, the timestamps and the version string. A bit-exact proof needs a clean rebuild and a binary comparison, which was not authorised and not done. Thread-order non-determinism (measured relative noise 7.6e-12 between two stock runs of a tiny model on 2 threads) also means that a rerun is not bitwise reproducible even with the identical binary.

## Experimental (rejected research prototype) build
Not installed, never replaced the stock binary. Exported from the pristine tag (`git archive v4.13`), patched, built in a separate directory. Round 1: Slurm 14093 (4 cores, 5 min 10 s, 0.344 CPU-h), rebuild 14101; binary `85c344a7308ecf9e1bd27a287277ae85f868f8db39db43ab7cc2c3b2649ca8a8`, `libfebiofluid.so` `1bff35b8c64cb05a068a9a5867603927f007bfa44f7429c020b4ef76440dca27`. Round 2: rebuild 14117; `febio4` unchanged, `libfebiofluid.so` `6627dadd1db0900166263cf68925d67141801bc872bb3ce5086c2ce4603118bb`. Diffs and hashes: `../experimental/README.md`. The stabilisation is activated by the environment variable `FEBIO_PPP_ALPHA` (default 0 = stock formulation); inputs differ from the stock inputs only by one XML comment line.

## What is not included
No compiled binary, no libraries, no copy of the FEBio source tree. The FEBio source lines quoted in the audits are from the unmodified v4.13 tag and should be read from that tag.
