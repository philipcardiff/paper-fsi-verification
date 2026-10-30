> **ARCHIVED WORKING REPORT (verbatim except path/host sanitisation).** Written during the investigation; file paths such as `runs/`, `tools/`, `cavity_agent/` refer to the original working tree, not to this package (see `docs/INDEX.md` for the package equivalents). Later reports supersede earlier ones where they disagree; known numerical/wording inconsistencies and their resolution are listed in `provenance/INCONSISTENCIES.md`. Nothing here claims a verified FEBio result.

# FEBIO EXPERIMENTAL BUILD febio-v4.13-ppp (Stage 1/2), 2026-10-30
Experimental, NOT installed, never replaces the stock `febio4`. Label scheme: PROVEN / SUPPORTED / UNRESOLVED / FAILED / NOT TESTED.

## 1. Source
* Base: pristine tag `v4.13` (`git describe` = v4.13, HEAD 32ae206ff4881dfb54f62296cd1558e58ed9fcc6) exported with `git -C upstream/febio-v4.13 archive v4.13 | tar -x -C src_experimental/febio-v4.13-ppp` (no `.git`, so `SET_DEVCOMMIT` adds nothing). The original tree `upstream/febio-v4.13` was only read: `git status --short` = 0 files, `git diff` empty, stock binary `build/febio-4.13/bin/febio4` sha256 f5ae1ed6...72f02 and `libfebiofluid.so` 178ba14e...af2b unchanged after all work (PROVEN).
* Experimental tree: `src_experimental/febio-v4.13-ppp`. Exact diff (6 files, 266 lines, sha256 d43859c1fee504c7a430...): `cavity_agent/ppp_experiment/ppp_v413_experimental.diff` (diff -ruN against `upstream/febio-v4.13`).

| file | change |
|---|---|
| `FEBioFluid/FEFluidPPP.h/.cpp` (new) | pure functions `Residual`, `TangentE`, `TangentX`, `Alpha()` (env `FEBIO_PPP_ALPHA`, read once, default 0), `FirstUse()`; each returns immediately if alpha == 0 or mu <= 0 |
| `FEBioFluid/FEFluidDomain3D.cpp` | after the stock loops in `ElementInternalForce` and `ElementStiffness`: block guarded by `if (febio_ppp::Alpha() != 0.0)` adding the PPP residual to the continuity row (4i+3) and the e-e tangent; log banner once in `PreSolveUpdate` |
| `FEBioFluid/FEFluidFSIDomain3D.cpp` | same in the FSI domain (7i+6), plus the d/du tangent (current-configuration quadrature weights) added to kJu |
| `FEBioFluid/FEFluidNormalVelocity.cpp` | **source change #2**: with alpha > 0 the inlet rim pressure is re-evaluated once per time step instead of every Update (alpha == 0: stock call unchanged). Needed because the stock rim constraint is lagged in the Jacobian and the PPP term makes the lag O(1) (see derivation s.8) |
| `FEBioLib/version.cpp` | version string `4.13.0-ppp` (no physics) |
Stock arithmetic in `fe`/`ke` is untouched; the experimental code only adds after the loops, so alpha = 0 executes the stock instructions (PROVEN at element level bitwise; run level: 8e-12 = thread-order noise, Stage 3).

## 2. Build
* Compiler GNU 11.4.0 (`/usr/bin/c++`), CMake 3.22.1, Release (`-O3 -DNDEBUG`), `USE_MKL=ON` (MKLROOT `tools/mkl-root`, the same linker-script stubs as the stock build), OpenMP, no HYPRE/SuperLU/FFTW; `CMAKE_{EXE,SHARED}_LINKER_FLAGS=-Wl,-rpath,tools/mkl-venv/lib`: the options of the stock `build/febio-4.13/CMakeCache.txt`. Script `jobs/build/build_febio_ppp.sbatch` (Slurm 14093, partition main, 4 cores, `cmake --build . -j 4`), incremental rebuild after source change #2 `jobs/build/rebuild_febio_ppp.sbatch` (14101; 14100 failed on a missing include, 8 CPU-s).
* Decision: a FULL separate build in its own directory (`build/febio-4.13-ppp`, own object files and libraries). Reusing production objects was rejected: absolute paths in the CMake cache/dependency files make it fragile, and the saving (stock full build = 15.8 CPU-min user) was small against the 3 CPU-h budget.
* Build wall time 5 min 9 s on 4 cores (user 17.6 min, sys 2.4 min) = **0.344 CPU-h allocated** (1240 CPU-s); the rebuilds 16 CPU-s.
* Binary `build/febio-4.13-ppp/bin/febio4` sha256 `85c344a7308ecf9e1bd27a287277ae85f868f8db39db43ab7cc2c3b2649ca8a8` (unchanged by the relink); `lib/libfebiofluid.so` (final, with change #2) `1bff35b8c64cb05a068a9a5867603927f007bfa44f7429c020b4ef76440dca27` (first build baf746b4...); `libfebiolib.so` 62ee850e...; log header `version 4.13.0-ppp` (stock `4.13.0.32ae206ff`) and, when the library loads the domains, the banner `EXPERIMENTAL febio-v4.13-ppp ... alpha = <value>` is written into every run log.
* `ldd` (env -i with only the MKL venv in LD_LIBRARY_PATH): libfebiolib, libfecore, libnumcore, libfebioopt, **libfebiofluid**, libfebiorve, libfeimglib, libfebiomix, libfebioplot, libfebiomech, libfeamr, libfebioxml all resolve to `build/febio-4.13-ppp/lib/` (the binary's RUNPATH is `tools/mkl-venv/lib:$ORIGIN/../lib/:build/febio-4.13-ppp/lib`, no production path). Run logs record `febio_bin`/`febio_sha256` in `run_meta_start.json` (PROVEN). `nm -D` shows 5 `febio_ppp` symbols in the experimental libfebiofluid and 0 in the stock one.
* Run mechanism: `FEBIO_BIN_OVERRIDE=<exp binary> FEBIO_PPP_ALPHA=<value> tools/submit_guarded.sh ...` (sbatch exports the environment; verified by the banner and `run_meta_start.json`). The stabilisation is not an input tag; inputs differ from the stock inputs only by one XML comment line (so the submission guard's hash check passes), PROVEN by `diff` (C1b: identical except line 3).

## 3. Stage 2 test record
| test | tool | result |
|---|---|---|
| element residual/tangent, 12 cases (random rotation, 15 % non-affine distortion, aspect 1:0.4:2, 2^3 and 3^3 Gauss, alphaf 1 / 0.5, linear and quadratic P(e)) against the compiled library | `cavity_agent/ppp_experiment/test_ppp_element.cpp`, output `results/stage2_element_tests.txt` | alpha = 0 bitwise no-op (12/12); constant-p residual 2.8e-16 of scale; row sums 4.6e-15; **dS/de vs central FD 2.5e-13**; **dS/dx vs FD 3.5e-10** (relative to the largest entry); y-parity energy -(alpha/mu) a^2 V/3 to 2e-16, max nodal S 3.1e-4 m^3/s for +-0.1 Pa at h = 0.1 |
| alpha = 0 on a tiny fluid-only straight channel (1400 hex8, 10 steps, 2 threads) vs stock; noise-floor control = stock repeated | T1/T2 runs, `results/T1_stock_vs_*.json` | stock vs repeated stock 7.6e-12, stock vs experimental alpha = 0 7.6e-12, identical iteration sequences (57 = 57) |
| alpha = 3/4 as built first | T1_ppp_a075 | **FAILED**: Max nr of reformations (15) at step 1, Newton contraction 0.55 per iteration |
| diagnosis | T1d_norim, T1d_maxrefs60 | without the rim constraint: quadratic, 4 it/step; max_refs 60 with rim: converges at 43 it/step. Cause: lagged rim constraint, not an element error |
| after source change #2 | T2_ppp_a075 | quadratic, 4 it/step (stock 5-7) |
| steady straight channel, exact Poiseuille known | T3 (t = 100 s) | dp/dx -0.12117 (stock) -> -0.12026 (alpha 3/4) vs exact -0.12; profile error 1.01 % -> 0.26 %; y-alternation 6.3e-4 Pa -> 1e-8 Pa; **flux 1.0000 -> 0.99264 at x = 2, outlet 0.99708** (conservation defect, derivation s.9) |
| FSI domain, 5 steps of the coupled M1 start-up (C2 generator flags, t_end 12.5) | T4_fsi_* | alpha = 0 identical to stock (9e-12, 34 = 34 it); alpha = 3/4 converges (27 it); transient fields differ strongly (start-up, not steady) |
Not tested: bitwise identity of stock vs experimental alpha = 0 at run level (thread-order noise 1e-11 makes it impossible on >1 thread; a 1-thread run is refused by the submission guard, 2 cores minimum); `FEFluidDomain3D::ElementMassMatrix` (no PPP term needed, no change).

## 4. Running CPU-hour budget (allocated cores x wall; limit 3.0 CPU-h; from `sacct`, jobs of 2026-10-30 evening)
| item | job | cores x wall | CPU-h | cumulative |
|---|---|---|---|---|
| full experimental build | 14093 | 4 x 310 s | 0.344 | 0.344 |
| rebuild attempts | 14100 (failed include), 14101 | 2 x 4 s each | 0.004 | 0.349 |
| T1 stock / alpha 0 / alpha 3/4 (failed) | 14094-14096 | 2 x (9,9,2) s | 0.011 | 0.360 |
| diagnostics T1d norim / maxrefs60 | 14098, 14099 | 2 x (4,30) s | 0.019 | 0.379 |
| T2 stock-rep / alpha 0 / alpha 3/4 | 14102-14104 | 2 x (8,8,5) s | 0.012 | 0.391 |
| T3 stock / alpha 3/4 (t = 100 s) | 14105, 14106 | 2 x (62,26) s | 0.049 | 0.440 |
| Stage 3 C1b alpha = 0 control | 14107 | 4 x 150 s | 0.167 | 0.607 |
| Stage 4 C1b alpha = 3/4 | 14108 | 4 x 48 s | 0.053 | 0.660 |
| T4 FSI-domain start-up tests stock / 0 / 3/4 | 14109-14111 | 2 x (15,15,9) s | 0.022 | 0.682 |
| **Total** | | | **0.68 CPU-h of 3.0** (sacct total 2452 CPU-s) | Stage 5 not run |
Login-node work (python analysis, one g++ compile of the element harness, ~30 s) is not allocated and is not counted.

## 5. Flux-preserving round (round 2), budget 0.75 CPU-h (separate from the 3.0 of round 1, of which 0.68 was used)
* Source change (new diff vs pristine v4.13, cumulative: `cavity_agent/ppp_experiment/ppp_v413_round2_flux_preserving.diff`, 297 lines, sha256 71bffd9b2290...; round-1 diff kept as `ppp_v413_round1_alpha075.diff`): `FEFluidPPP.h/.cpp` gain the optional constrained-node mask `con` (redistribution S~_i = S_i + (1/nF) sum_C S_c, S~_c = 0, applied to residual, e-tangent and d/du-tangent rows) and `RimOnce()` (alpha > 0 or env FEBIO_PPP_RIM_ONCE=1); `FEFluidDomain3D.cpp`/`FEFluidFSIDomain3D.cpp` build the mask from `node.get_bc(m_dofEF) != DOF_OPEN`; `FEFluidNormalVelocity.cpp` uses `RimOnce()`. alpha = 0 skips everything. The round-1 alpha = 3/4 behaviour (mask = nullptr) is no longer selectable in the binary (the mask is always passed): the round-1 results stay valid for the round-1 library (baf746b4.../1bff35b8...).
* Binary `build/febio-4.13-ppp/bin/febio4` sha256 85c344a7308ecf9e1bd27a287277ae85f868f8db39db43ab7cc2c3b2649ca8a8 (unchanged); `lib/libfebiofluid.so` **6627dadd1db0900166263cf68925d67141801bc872bb3ce5086c2ce4603118bb**; ldd: libfebiofluid and all other libraries resolve to `build/febio-4.13-ppp/lib` (checked with env -i). Stock tree: 0 changed files; stock binary f5ae1ed6... unchanged.
| item | job | cores x wall | CPU-h | cumulative |
|---|---|---|---|---|
| incremental rebuild (6 s) | 14117 | 2 x 6 s | 0.003 | 0.003 |
| T5 rim-once control, alpha = 0 (t = 100 s tiny channel) | 14118 | 2 x 40 s | 0.022 | 0.025 |
| T5 alpha = 1e-3 tiny Poiseuille channel | 14119 | 2 x 36 s | 0.020 | 0.046 |
| **Round-2 total (sacct)** | | | **0.046 of 0.75** | |
Tests run on the login node (seconds, not counted): `test_ppp_conservation.cpp`, `test_ppp_element.cpp`, `t5_analysis.py`.
