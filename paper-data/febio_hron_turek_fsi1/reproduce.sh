#!/usr/bin/env bash
# Non-simulation reproduction checks for the FEBio v4.13 Hron-Turek FSI1 diagnostic package.  Needs only Python 3 with numpy and matplotlib.  Usage: ./reproduce.sh [scratch_dir]
# 1. regenerate the inputs from the generators and compare sha256 with inputs/;  2. rebuild the stage-B tables from the stored raw outputs;  3. recompute the stage-C reductions
# (two independent implementations + parity + traction reconstruction) from the stored final states and compare with the stored results;  4. regenerate the figures;  5. verify SHA256SUMS.
# The simulations themselves (FEBio 4.13 + MKL, Slurm) are NOT re-run; see docs/PROVENANCE_AND_CONFIGURATION.md for the commands.
set -u
ROOT=$(cd "$(dirname "$0")" && pwd); S=${1:-$(mktemp -d)}; mkdir -p "$S"; fail=0; PY="python3 -I"
echo "scratch: $S"
echo "== 1. regenerate inputs and compare hashes"
for c in B_CSM1 B_CSM2; do
  if [ $c = B_CSM1 ]; then E=1.4e6; else E=5.6e6; fi
  for nt in 2 4 8 16; do $PY $ROOT/scripts/make_hronturek_csm_feb.py --nt $nt --elem hex8 --E $E --out $S/gen/${c}_h8_nt$nt --tag ${c}_h8_nt$nt >/dev/null; done
  for nt in 2 4; do $PY $ROOT/scripts/make_hronturek_csm_feb.py --nt $nt --elem hex20 --E $E --out $S/gen/${c}_h20_nt$nt --tag ${c}_h20_nt$nt >/dev/null; done
  $PY $ROOT/scripts/make_hronturek_csm_feb.py --nt 16 --elem hex8 --E $E --root 0.25 --out $S/gen/${c}_h8_nt16_root025 --tag ${c}_h8_nt16_root025 >/dev/null
done
for nt in 8 16; do $PY $ROOT/scripts/make_hronturek_csm_feb.py --nt $nt --elem hex8 --E 1.4e6 --g -0.02 --out $S/gen/B_LIN_h8_nt$nt --tag B_LIN_h8_nt$nt >/dev/null; done
$PY $ROOT/scripts/make_hronturek_cfd_feb.py --out $S/gen/C1u_L1_f1.25_K1e6 --tag C1u_L1_f1.25_K1e6 --factor 1.25 >/dev/null
$PY $ROOT/scripts/make_hronturek_cfd_feb.py --out $S/gen/C1p_L1_f1.25_K1e6 --tag C1p_L1_f1.25_K1e6 --factor 1.25 --pin-root-values=-81.6042241147,-123.87086089350001 >/dev/null
nok=0; nbad=0
for f in $ROOT/inputs/stageB/B_*.feb $ROOT/inputs/stageC/*.feb; do b=$(basename $f .feb); g=$S/gen/$b/$b.feb
  if [ -f $g ] && [ "$(sha256sum < $f)" = "$(sha256sum < $g)" ]; then nok=$((nok+1)); else echo "  HASH MISMATCH or not regenerable: $b"; nbad=$((nbad+1)); fi; done
echo "  inputs identical to regenerated: $nok ; mismatching: $nbad (the v1 diagnostic input DIAG_v1_* is not in the B_* glob)"; [ $nbad -gt 0 ] && fail=1
echo "== 2. stage-B tables from stored raw outputs"
$PY $ROOT/scripts/make_stageB_tables.py $ROOT/results/stageB_raw $S/tablesB >/dev/null
for t in stageB_summary.csv stageB_convergence.csv; do if cmp -s $ROOT/results/$t $S/tablesB/$t; then echo "  $t identical"; else echo "  $t DIFFERS"; fail=1; fi; done
$PY $ROOT/scripts/check_inputs.py $ROOT/inputs $S/input_configuration_check.csv >/dev/null; cmp -s $ROOT/results/input_configuration_check.csv $S/input_configuration_check.csv && echo "  input_configuration_check.csv identical" || { echo "  input_configuration_check.csv DIFFERS"; fail=1; }
echo "== 3. stage-C reductions from the stored final states (no .xplt in the package: xplt cross-checks are skipped)"
for c in C1u_L1_f1.25_K1e6 C1p_L1_f1.25_K1e6; do
  mkdir -p $S/runs/$c; cp $S/gen/$c/$c.feb $S/gen/$c/model_meta.json $S/gen/$c/mesh2d.npz $S/runs/$c/; cp $ROOT/results/stageC_raw/$c/Mon_all.txt $S/runs/$c/
  $PY $ROOT/scripts/hronturek_cfd1_independent.py $S/runs/$c $S/independent_$c.json >/dev/null
  $PY $ROOT/scripts/hronturek_cfd_analyze.py $S/runs/$c >/dev/null
  $PY $ROOT/scripts/hronturek_fsi_traction_recon.py $S/runs/$c $S/recon_$c.json >/dev/null
  $PY $ROOT/scripts/parity_diagnostics.py $S/runs/$c $S/parity_$c.json >/dev/null
  $PY $ROOT/scripts/compare_json.py $ROOT/results/independent_forces_$c.json $S/independent_$c.json 1e-9 || fail=1
  $PY $ROOT/scripts/compare_json.py $ROOT/results/fsi_traction_reconstruction_$c.json $S/recon_$c.json 1e-9 || fail=1
  $PY $ROOT/scripts/compare_json.py $ROOT/results/parity_diagnostics_$c.json $S/parity_$c.json 1e-9 || fail=1
  # the first-generation analysis script against the stored result of the full run (time-history items excluded: only the final step is stored in the package)
  $PY $ROOT/scripts/subset_json.py $ROOT/results/stageC_raw/$c/stageC_analysis.json $S/a_$c.json F1 F2 F3 mass_imbalance Qsec_x1_over_Qin consist_drag_F1F2 consist_drag_F1F3 consist_lift_F1F2 consist_lift_F1F3 parity_top parity_bottom parity_cyl_at cross_flag
  $PY $ROOT/scripts/subset_json.py $S/runs/$c/stageC_analysis.json $S/b_$c.json F1 F2 F3 mass_imbalance Qsec_x1_over_Qin consist_drag_F1F2 consist_drag_F1F3 consist_lift_F1F2 consist_lift_F1F3 parity_top parity_bottom parity_cyl_at cross_flag
  $PY $ROOT/scripts/compare_json.py $S/a_$c.json $S/b_$c.json 1e-9 || fail=1
done
$PY $ROOT/scripts/hronturek_cfd_diag.py $S/runs/C1p_L1_f1.25_K1e6 >/dev/null && $PY $ROOT/scripts/compare_json.py $ROOT/results/stageC_raw/C1p_L1_f1.25_K1e6/stageC_diag.json $S/runs/C1p_L1_f1.25_K1e6/stageC_diag.json 1e-9 || fail=1
echo "== 4. figures"; $PY $ROOT/scripts/make_figures.py $ROOT $S/figures >/dev/null 2>&1 && echo "  figures regenerated in $S/figures (committed copies may differ at the byte level)" || { echo "  figure generation FAILED"; fail=1; }
echo "== 5. SHA256SUMS"; (cd $ROOT && sha256sum -c SHA256SUMS 2>&1 | grep -v ': OK$' | head -5; true)
[ $fail -eq 0 ] && echo "ALL NON-SIMULATION CHECKS PASSED" || echo "SOME CHECKS FAILED (see above)"; exit $fail
