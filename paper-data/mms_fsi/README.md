# Manufactured solution for transient ALE FSI with a finite-strain solid

The case is a 2-D channel flow coupled to a compliant neo-Hookean wall.

- The solution is closed form.
- Both interface conditions hold exactly, so the interface needs no source term.
- There are two variants:
  - `full`: the interface tilts, travels and stretches;
  - `flat`: the interface moves only vertically, so the ALE mesh stays orthogonal at the interface.

The construction is described in `paper/mms_construction.tex`.

## Contents

```
derive/   mms_fsi.py          sympy/symengine construction, 11 exactness checks at 40 digits, C++ generation
          codegen.py          common-subexpression elimination -> mmsFsiExact.H
          check_generated.C   checks of the generated header, independent of the symbolic pipeline
          run_checks.sh       derive + generate + check: run_checks.sh <outdir> <rho_f> <rho_s> <full|flat> [mu_f]
          generated/          full variant (mmsFsiExact.H, reference_values.txt), with logs generated.*.log
          generated_flat/     flat variant
          generated_rho/      full variant with rho_F = 2, rho_S = 3 (density-scaling check only)
          generated_mu10/     full variant with mu_F = 1/10 (diagnostic: stronger viscous traction)
case/     src/                case libraries libmmsFsi, libmmsFsiFlat, libmmsFsiMu10:
                                source fvOption, exact boundary conditions, error function object;
                              evalExact.C (exact values at given points, used by timeOrders.py)
          makeCase.py         writes a case: coupling, variant, mesh level, time step, options
          orders.py           observed orders from the error files of a set of runs
          timeOrders.py       temporal orders from differences of successive time-step solutions
          foamio.py           minimal field reader
results/  order tables of the runs quoted in the papers (24 Nov 2026; see the build notes below)
paper/    mms_construction.tex
```

## Requirements

- **Derivation:** Python 3 with sympy and symengine (used: Python 3.14.3, sympy 1.14.0, symengine 0.15.0), and a C++17 compiler.
  - Each variant takes about 6 minutes.
- **Cases:** OpenFOAM v2512 and solids4foam built with PETSc (used: PETSc 3.22.2).
  - The quasi-monolithic (QM) runs need solids4foam with #751 (fvOptions in the QM solvers). This is merged in `feature-petsc-snes-quasi-monolithic` db3039b76.
  - The partitioned runs (`pimpleFluid` + displacement solid) need no solver change.

## Running

```bash
# 1. Optional: regenerate and check the exact solution (the generated headers are included)
derive/run_checks.sh generated 1 1 full
derive/run_checks.sh generated_flat 1 1 flat

# 2. Environment file: loads OpenFOAM and sets WM_PROJECT_USER_DIR to the solids4foam build.
#    The case libraries go to its FOAM_USER_LIBBIN.
export MMS_FSI_ENV=/path/to/env.sh
case/src/Allwmake $MMS_FSI_ENV
(cd case && clang++ -O1 -std=c++17 -I ../derive/generated src/evalExact.C -o evalExact \
         && clang++ -O1 -std=c++17 -I ../derive/generated_flat src/evalExact.C -o evalExactFlat)

# 3. A QM spatial study at a fixed small time step
cd case
for n in 16 32 64 128; do
  python3 makeCase.py runs/space-full/n$n --coupling qm --variant full --nx $n --dt 0.000625 --end 0.4 --writeInterval 0.2
  runs/space-full/n$n/Allrun
done
python3 orders.py runs/space-full/n{16,32,64,128} --times 0.2 0.4

# 4. A temporal study on n32: the exact solution is subtracted at each run's own cell centres
for dt in 0.04 0.02 0.01 0.005 0.0025 0.00125; do
  python3 makeCase.py runs/time-full/dt$dt --coupling qm --variant full --nx 32 --dt $dt --end 0.4 --writeInterval 0.2
  runs/time-full/dt$dt/Allrun
done
python3 timeOrders.py runs/time-full/dt{0.04,0.02,0.01,0.005,0.0025,0.00125} --time 0.4 --exact ./evalExact
```

**Studies.**

| Study | Δt | Options |
|---|---|---|
| Combined | 0.32/n | — |
| Partitioned | 0.32/n | `--coupling aitken\|robin\|iqnils --outerTol 1e-6` |

## Main options of `makeCase.py`

| Option | Effect |
|---|---|
| `--coupling` | `qm`, `aitken`, `robin`, `iqnils`. `fluid`: fluid alone with the exact interface motion. `solid`: solid alone with the exact interface traction |
| `--variant` | `full`, `flat`, `mu10` |
| `--fluidModel` | fluid-only mode: `newtonIcoFluid` (default) or `pimpleFluid` |
| `--outlet` | `neumannU`: exact velocity gradient and exact pressure (default). `dirichletU`: exact velocity and exact pressure |
| `--ddtCorr`, `--ddtPhiCoeff`, `--pimpleGrad`, `--interfaceP` | `pimpleFluid` settings |
| `--fsiExtra` | extra entries for the partitioned coupling coefficients, e.g. `"requireAllResidualMeasures yes;"` |
| `--rcHistory` | QM: `RhieChowTimeHistory yes`. Needs solids4foam PR #752 |
| `--grad leastSquaresS4f --uCorrected --pCorrected` | QM: the corrected-boundary recipe. Needs PR #754; without it, `fixedGradientCorrected` is not a momentum-balance pressure |
| `--solidSurfGrad` | QM: `tractionPatchSurfaceGradient yes`. Needs PR #756 |
| `--pOracle` | QM diagnostic: exact interface pressure value or exact k·∇p |
| `--np` | MPI ranks (simple decomposition in x) |

## Error output

`postProcessing/mmsFsiErrors/errors.dat` is written every time step. It holds the volume-weighted L1, L2 and L∞ errors of:

- the fluid U and p (on the current mesh);
- the solid D and U (on the initial mesh);
- the interface displacement and the traction applied to the solid (face-area weighted on the initial solid interface);
- the fluid interface pressure.

## Build notes for `results/`

- **`main_*`:** solids4foam `feature-petsc-snes-quasi-monolithic` db3039b76, OpenFOAM v2512, PETSc 3.22.2, macOS arm64.
- **`openpr_*`:** a local test merge, 682acedb3. It is db3039b76 plus three open PRs:
  - #752 (298ec4852, `RhieChowTimeHistory`);
  - #754 (f4c31f621, corrected-BC recipe);
  - #756 (41fa464c3, `tractionPatchSurfaceGradient`);

  plus one line that makes #752's `storeRhieChowHistory` call #754's `updatePressureGradient()`.

  The case options are `--grad leastSquaresS4f --uCorrected --pCorrected --solidSurfGrad --rcHistory`.
- **`p750_*`:** db3039b76. These are diagnostics of the partitioned (`pimpleFluid`) fluid; see solids4foam#750.

## Results of the verification paper: provenance

The files below support Section "Manufactured solution for transient ALE FSI with a finite-strain solid" (`sec:mms_fsi`, Table `tab:mms_fsi`) and Appendix `app:mms` of the paper.

### Problem and runs

- 2-D channel with a compressible neo-Hookean wall (plane strain), closed-form manufactured solution with exact kinematic and traction continuity; see Appendix `app:mms`.
- **Full variant:** tilting, travelling, tangentially stretching interface. **Flat variant:** the interface stays flat and moves only vertically.
- Meshes: fluid n × n/2, solid n × n/4, square cells, n = 16, 32, 64, 128 faces on the interface. Combined refinement Δt = 0.32/n (0.02 to 0.0025). Run to t = 0.4.
- Errors are logged every step by the case's error function object:
  - `U_L2`, `p_L2`: volume-weighted RMS of |U − U_exact| and |p − p_exact| at the fluid cell centres, on the moving mesh (pressure not mean-corrected; the outlet fixes the level).
  - `D_L2`, `Us_L2`: solid displacement and velocity on the reference mesh.
  - `intD_L2`, `intTrac_L2`: face-area-weighted RMS of the interface displacement and of the traction applied to the solid, against the exact traction, on the solid interface.
  - `intP_L2`: fluid interface pressure. `*_L1`, `*_Linf` as named.
- Codes: OpenFOAM-v2512, macOS (arm64), local solids4foam worktree (branch `mms-fsi/qm`, later `mms-fsi/item4`; the quasi-monolithic source-term support has since been merged into solids4foam).
  - The partitioned n16–n64 runs used commit 4071b431e and the n128 runs commit 677865f55. The source changes between the two touch only the quasi-monolithic fluid, solid and coupling (`newtonIcoFluid`, `nonLinGeomTotalLagVelocitySolid`, `newtonQuasiMonolithicCouplingInterface`), not `pimpleFluid`, the displacement solid or the partitioned couplings.
  - The case files of n128 differ from n64 only in the block divisions and Δt (checked with `diff -r`).
- Couplings:
  - `aitken_full`: Dirichlet–Neumann, Aitken, `newMovingWallVelocity` + `fixedFluxPressure` on the interface.
  - `robin_full`: Robin–Neumann (`fixedRelaxation`, factor 1; `elasticWallVelocity` + `elasticWallPressure`).
  - Both: tolerance 1e-6 with the legacy (min) test, at most 200 iterations. At acceptance r1 ≈ 1e-4 (Aitken) and ≈ 1e-6 (Robin–Neumann).
  - `qm_full`, `qm_flat`: quasi-monolithic coupling of the companion paper (same meshes, mesh motion and Δt).
- A/B on Aitken, full variant, n16–n64 (2026-11-21):
  - `aitken_full_ab_S`: `flux(HbyA) skewCorrected linear`.
  - `aitken_full_ab_SM`: S + `movingWallPressure` on the interface.
  - `aitken_full_ab_SW`: S + the exact pressure gradient on the bottom wall.
  - None changes the fluid U/p orders. S raises the interface traction order (n32 to n64) from 1.41 to 1.61.

### Files

| File | Content | Origin |
|---|---|---|
| `summary_errors_t0.2_t0.4.csv` | L2 (and U, p L∞) errors at t = 0.2 and 0.4 for every series and level; the numbers of Table `tab:mms_fsi` | `make_csv.py` |
| `coupling_iterations.csv` | coupling iterations per step (mean, max) for the partitioned runs | `make_csv.py`, from `postProcessing/fsiResiduals.dat` |
| `errors_csv/<series>_n<n>.csv` | the full error history (every time step) of each run | `make_csv.py`, from `postProcessing/mmsFsiErrors/errors.dat` |
| `order_tables/{aitken,robin}_full_n16-n128.md` | four-level order tables of the partitioned runs, t = 0.2 and 0.4 | `case/orders.py <n16> <n32> <n64> <n128> --times 0.2 0.4` |
| `order_tables/aitken_full_ab_*.md` | the A/B order tables, n16–n64 | `case/orders.py` |
| `source_tables/*.md` | the order tables of the original study, copied unchanged from `case/results/` | |
| `derive/generated*.derive.log` | output of the derivation: the 11 symbolic and 40-digit checks (full, flat, ρ_F = 2 and ρ_S = 3, μ_F × 10) | `derive/run_checks.sh` |
| `derive/generated*.check.log` | output of the independent checks of the generated C++ | `derive/check_generated.C` |
| `results/` | the final order tables of 24 Nov 2026 (QM main build and corrections under review, partitioned-fluid diagnostics `p750_*`); see the build notes above | `case/orders.py`, `case/timeOrders.py` |
| `make_csv.py` | the script that wrote the CSV files (absolute paths of the run directories inside) | |
| `figure/mms_fsi_case_figure.py` | Fig. `fig:mms_fsi_case` (`figures/mms_fsi_case.pdf`): the exact solution of the full variant at t = 0.1, 0.2, 0.3, 0.4 (fluid speed and streamlines, solid strain and deformed material lines, current configuration, to scale), the interface displacement η(x, t) over one period, and the full and flat interfaces at t = 0.2 and 0.4. See "Figure of the case" below | the exact evaluators only; no run data |
| `figure/mms_fsi_figure.py` | Fig. `fig:mms_fsi` (`figures/mms_fsi_pressure_error.pdf`): the pressure error p − p_exact at t = 0.4 on n128 for Aitken and QM. Run with Python 3 (numpy, matplotlib): `mms_fsi_figure.py mms-pfv/base/aitken-n128 mms-fsi/case/runs/qm-full/n128 OUT.pdf`. The runs wrote fields only at t = 0.4 | reads the `pError` field of the error function object and the deformed fluid mesh |

**Figure of the case** (`figure/mms_fsi_case_figure.py`). Python 3 with numpy and matplotlib; the fields come from the exact evaluators of `case/src/evalExact.C`:

```bash
cd case
clang++ -O2 -std=c++17 -I ../derive/generated      src/evalExact.C -o evalExact
clang++ -O2 -std=c++17 -I ../derive/generated_flat src/evalExact.C -o evalExactFlat
python3 ../figure/mms_fsi_case_figure.py ./evalExact ./evalExactFlat mms_fsi_case.pdf \
        [--times 0.1 0.2 0.3 0.4] [--compare 0.2 0.4] [--tmax 1] [--trun 0.4] \
        [--dpsi 0.06] [--width 5.17] [--png preview.png]
```

- Fluid: |v_F| on a grid (x, s) with y = s h(x, t), where h is interpolated from the exact interface curve. The streamlines are contours of the stream function ψ = ∫₀ʸ v_F,x dy′ (ψ = 0 on the wall), integrated by the trapezoidal rule (difference from a grid four times finer: 2e-4, against a contour spacing of 0.06), plus one closed contour around the recirculation centre.
- Solid: the Frobenius norm of the Green–Lagrange strain, E = (FᵀF − I)/2 with F_zz = 1, where F = I + (∇₀D_S)ᵀ comes from central differences (step 1e-6) of the exact displacement; the deformed material lines are the exact images of X = const and Y = const.
- The script prints, per snapshot, max |v_F|, max ‖E‖, the range of J and the interface height; at t = 0.4: 2.15, 0.263, 0.690 ≤ J ≤ 1.035, η ≤ 0.124.

**`source_tables/`:**

- `qm_combined_{full,flat}.md`: QM, Δt ∝ h (QM rows of Table `tab:mms_fsi`).
- `qm_space_{full,flat}.md`: QM, fixed Δt = 6.25e-4.
- `qm_time_*.md`: QM temporal differences on n32.
- `{aitken,robin}_combined_full.md`: partitioned n16–n64. Superseded by the four-level tables in `order_tables/`.

**Run directories:**

- partitioned n16–n64 and all QM runs: `claude-runs/mms-fsi/case/runs/{aitken-full,robin-full,qm-full,qm-flat}/n*`;
- partitioned n128: `claude-runs/mms-pfv/base/{aitken,robin}-n128`;
- A/B: `claude-runs/mms-pfv/{S,SM,SW}/n*`.

These directories are not archived here.

**Checks.** The earlier `checks/` folder is superseded by the logs in `derive/`. The ρ_F = 2, ρ_S = 3 set was regenerated on 24 Nov 2026 after the outlet-flux redesign (its velocity is now density-independent, as it should be); all 11 + 10 checks pass for every set.

**IQN-ILS runs** (Section `sec:mms_fsi`, paragraph IQN-ILS): `claude-runs/mms-fsi/final/iqn/` on the Mac (not archived here); full variant, Δt = 0.32/n, tolerance 1e-6, `couplingReuse 0`.
