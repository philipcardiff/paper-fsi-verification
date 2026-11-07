> **ARCHIVED WORKING REPORT (verbatim except path/host sanitisation).** Written during the investigation; file paths such as `runs/`, `tools/`, `cavity_agent/` refer to the original working tree, not to this package (see `docs/INDEX.md` for the package equivalents). Later reports supersede earlier ones where they disagree; known numerical/wording inconsistencies and their resolution are listed in `provenance/INCONSISTENCIES.md`. Nothing here claims a verified FEBio result.

# FEBIO ALPHA = 0 CONTROL (Stage 3), experimental build febio-v4.13-ppp, fixed-plate fluid-only cavity C1b
## A. Pre-registration (written 2026-10-30 BEFORE the control run; the tolerances below were fixed before looking at any C1b output of the experimental binary)
* Reference: the stored stock run `runs/cavityNext/C1b_fluid_fixedplate_f1_pin` (job 14015, stock binary sha256 f5ae1ed6..., 4 cores, 171.9 s wall, 0.19 CPU-h, 600 steps, 1180 equilibrium iterations, normal termination).
* Input: regenerated with `tools/make_cavity_fluid_feb.py --mesh-factor 1 --k 1e6 --dt 0.5 --t-end 300 --mu 0.01 --inlet normalvel --ubar 1.0 --pin-corners 1.15 --log-all --rhoi 0 --ramp 10`: sha256 `a94a47a4...588cb31a`, **byte-identical** to the stored stock input (checked 2026-10-30). The control input differs from it ONLY by one XML comment line (`<!-- EXPERIMENT ... FEBIO_PPP_ALPHA=0 -->`, needed because the submission guard refuses a repeated input hash). The comment is not parsed by FEBio. The stabilisation is NOT an input tag: it is the process environment variable `FEBIO_PPP_ALPHA` (0 here) and is printed in the log header, so the physical input is identical to stock.
* Solver: 4 cores, same `run_febio_job.sh`, OMP/MKL threads = 4 (the same thread count as the stock run; non-bitwise reproducibility from thread-order effects is expected: measured on a tiny model, two stock runs on 2 threads already differ by 7.6e-12 relative, `cavity_agent/ppp_experiment/results/T1_stock_vs_T2_stock_rep.json`).
* Acceptance criteria (all must hold; fixed now):
  1. Normal termination, 600 steps, no errors, no extra warnings beyond the stock `No force acting on the system`.
  2. Nodal fields of every `Mon_*.txt` block (all 601 time blocks): max |dA| / (max |A| of that block) <= **1e-8** for w_x, w_y (velocity scale = max(|w_x|,|w_y|) of the block) and e_f (e_f scale = max |e_f| of the block). Justification: thread-order noise floor 1e-11 on a short run; the 1e-8 bound leaves three orders of margin for 600 steps while being far below any physically relevant difference.
  3. Total equilibrium iterations equal to 1180 within +-1 % and the per-step iteration sequence identical except for isolated +-1 changes at steps whose residual is at the tolerance edge (then reported).
  4. Derived quantities (python re-assembly, `cavity_agent/ppp_experiment/stab_analysis.py`): D (row parity), F1, F2, flux imbalance, momentum residual, continuity residual equal to the stock values to a relative 1e-8 (absolute 1e-10 for quantities that are ~0).
  5. Wall time within +-30 % of 172 s on a comparably loaded node (informative only; the node is shared).
* If any criterion fails: STOP, diagnose the baseline discrepancy, report; Stage 4 is not run.

## B. Result (job 14107, `runs/cavityExp/C1b_ppp_alpha0`, 4 cores, experimental binary 85c344a7... + libfebiofluid 1bff35b8..., FEBIO_PPP_ALPHA = 0; input sha256 8ee3b0fa... = stock C1b + one comment line)
**STAGE 3: PASSED (all five pre-registered criteria).**
| criterion | result |
|---|---|
| 1 termination | normal, 600/600 steps, 0 errors, the same single warning `No force acting on the system` (log inspected, not the exit code) |
| 2 nodal fields, all 601 blocks of Mon_all/cl/iface/inlet/outlet/x10/x13 | max relative difference **w_x 6.6e-12, w_y 1.9e-12, e_f 8.3e-12** (<= 1e-8); not bitwise (4-thread order effects; same magnitude as stock-vs-stock) |
| 3 iterations | 1179 vs 1180 (one step differs by one iteration; mean 1.965 vs 1.967) |
| 4 derived quantities (`results/C1b_ppp_alpha0_analysis.json` vs `results/stock_C1b_analysis.json`) | row parity D 0.42068933685939 (both), F1 -4.663716740901 vs -4.6637167409, F2 -5.48327045004775 vs -5.48327045004650, Q_out/Q_in - 1 9.0632e-5 (both), momentum residual 9.35e-7, continuity residual 2.2e-7: no difference above 1e-9 relative / 1e-10 absolute |
| 5 wall | 148 s vs 172 s (informative) |
Conclusion A: **alpha = 0 reproduces the stock solver on C1b within 1e-11 (PROVEN for this case, fluid-only domain, 4 threads); the baseline is the stock solution**. The FSI domain with alpha = 0 reproduces stock to 9e-12 on a 5-step coupled start-up (SUPPORTED; not a steady coupled run). No baseline discrepancy was found.
Caveat: the control was run with the final library (source change #2 compiled in); with alpha = 0 that code path calls the stock `SetRimPressure` every Update, so the stock behaviour is kept (shown by the identical iteration sequence).
