#!/usr/bin/env bash
# Non-simulation reproduction checks for the FEBio v4.13 radialTorsionalRing (torsion member) evidence package.
# Needs Python 3 with numpy, scipy and matplotlib, plus xz and sha256sum.  usage: ./reproduce.sh [scratch_dir]   (SKIP_REF=1 skips the spectral-reference regeneration, SKIP_AE=1 the analytical coupled estimate)
# 1. verify the hashes of the vendored reference;  2. regenerate every FEBio input from scripts/make_rtr_feb.py and compare sha256 with inputs/INPUT_HASHES.csv;
# 3. re-run the campaign analysis scripts and the independent re-derivation on the stored raw outputs and compare with results/;  4. recompute the log statistics;
# 5. regenerate tables and figures;  6. regenerate the spectral fixed-geometry reference from the omega table and compare with the cached arrays;  7. re-run the analytical coupled estimate;  8. verify SHA256SUMS.
# The simulations themselves (FEBio 4.13 + MKL under Slurm) are NOT re-run; see provenance/REPRODUCTION_RECORD.md for the commands.  Nothing is written inside the package: everything is done in a scratch copy.
set -u
SRC=$(cd "$(dirname "$0")" && pwd); S=${1:-$(mktemp -d)}; mkdir -p "$S"; fail=0; PY="python3 -I"
echo "scratch: $S"; rm -rf "$S/pkg"; cp -r "$SRC" "$S/pkg"; R="$S/pkg"
echo "== 1. vendored reference hashes"
$PY - "$R" <<'PYEOF' || fail=1
import hashlib, sys
r = sys.argv[1]; ok = True
for f, h in (("rt_reference_s4f637_30ae00eb.py", "f464e35ebd5a1da87fcb52966cfe7868811e7521d2a615a4b16fa7e9e3193ba3"), ("make_reference_s4f637_30ae00eb.py", "323207dc64fcc279cdd11376a64de445e4280c87cfacae95d2b132474f77d80c"), ("omega_B_table_dt1e-3.csv", "daab5828b027d1d1cae7fca6f7c00f6be3ec2a0890a28b881a5595a932e0f889")):
    g = hashlib.sha256(open(f"{r}/reference_data/rtr/{f}", "rb").read()).hexdigest(); print("  ", f, "OK" if g == h else "MISMATCH " + g); ok &= g == h
sys.exit(0 if ok else 1)
PYEOF
echo "== 2-4. inputs, campaign analysis, independent re-derivation, log statistics"
bash "$R/scripts/run_analyses.sh" "$R" "$S/work" "$S/out" >/dev/null 2>"$S/analyses.err" || { echo "  run_analyses.sh FAILED (see $S/analyses.err)"; fail=1; }
cat "$S/out/regen_inputs.txt" 2>/dev/null; grep -q "MISMATCH" "$S/out/regen_inputs.txt" 2>/dev/null && fail=1
for f in $(cd "$R/results/campaign_analysis" && ls *.json); do printf "  campaign_analysis/%s: " $f; $PY $R/scripts/compare_json.py $R/results/campaign_analysis/$f $S/out/campaign_analysis/$f 1e-9 || fail=1; done
printf "  independent_check.json: "; $PY $R/scripts/compare_json.py $R/results/independent_check.json $S/out/independent_check.json 1e-9 || fail=1
for f in run_summary.csv newton_iterations_per_step.csv rtc_m1_step2_residuals.csv; do if cmp -s $R/results/$f $S/out/$f; then echo "  $f identical"; else echo "  $f DIFFERS"; fail=1; fi; done
echo "== 5. tables and figures"
$PY $R/scripts/make_tables_figures.py $R $S/out $S/tables $S/figures >/dev/null || fail=1
for f in $(cd $R/results && ls rta_*.csv rtb_*.csv); do if cmp -s $R/results/$f $S/tables/$f; then echo "  $f identical"; else echo "  $f DIFFERS"; fail=1; fi; done
for f in $(cd $R/figures && ls *.png); do [ -s $S/figures/$f ] && echo "  figure $f regenerated ($(stat -c %s $S/figures/$f) bytes; stored $(stat -c %s $R/figures/$f); PNG bytes may differ across matplotlib versions)" || { echo "  figure $f MISSING"; fail=1; }; done
if [ "${SKIP_REF:-0}" != 1 ]; then
echo "== 6. spectral fixed-geometry reference regenerated from the omega table (a few minutes)"
mkdir -p $S/refcheck/reference_data/rtr $S/refcheck/scripts $S/refcheck/results
cp $R/reference_data/rtr/{rt_reference_s4f637_30ae00eb.py,omega_B_table_dt1e-3.csv,make_reference_s4f637_30ae00eb.py} $S/refcheck/reference_data/rtr/; cp $R/scripts/rtr_reference.py $R/scripts/rtr_analyze.py $R/scripts/rtr_check_rta.py $S/refcheck/scripts/
$PY - $S/refcheck $R <<'PYEOF' || fail=1
import sys, numpy as np
sys.path.insert(0, sys.argv[1] + "/scripts")
import rtr_analyze as A
tt = np.round(np.arange(0, 201) * 0.01, 10); ok = True
for shift, name in ((0.0, "fluid_exact_N40_dt0.01.npz"), (0.005, "fluid_exact_N40_dt0.01_shift0.005.npz")):
    new = A.exact_fields(tt, shift); old = np.load(sys.argv[2] + "/reference_data/rtr/" + name, allow_pickle=True)
    for k in ("T_B", "H_f", "L", "p"):
        d = np.max(np.abs(new[k] - old[k])) / max(np.max(np.abs(old[k])), 1e-300); print(f"   shift {shift} {k}: max relative difference {d:.2e}"); ok &= d < 1e-9
sys.exit(0 if ok else 1)
PYEOF
fi
if [ "${SKIP_AE:-0}" != 1 ]; then
echo "== 7. analytical coupled estimate (reduced model with element-mean shear stress)"
RTR_PREDICT_OUT=$S/predicted.json $PY $R/scripts/rtr_predict_coupled.py >/dev/null && $PY $R/scripts/compare_json.py $R/results/predicted_coupled_element_mean_AE.json $S/predicted.json 1e-6 || fail=1
fi
echo "== 8. SHA256SUMS"
(cd "$SRC" && sha256sum -c --quiet SHA256SUMS) && echo "  all files match SHA256SUMS" || fail=1
[ $fail = 0 ] && echo "REPRODUCTION CHECKS PASSED" || echo "REPRODUCTION CHECKS: FAILURES (see above)"
exit $fail
