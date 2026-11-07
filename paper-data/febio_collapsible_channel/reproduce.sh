#!/usr/bin/env bash
# Non-simulation reproduction checks for the FEBio v4.13 collapsibleChannel evidence package.  Needs Python 3 with numpy, scipy (only for step 6) and matplotlib.
# Usage: ./reproduce.sh [scratch_dir]      (FULL=1 also runs the 40-element beam re-derivation, about 1.5 min on one core)
# 1. regenerate all FEBio inputs from scripts/make_cc_feb.py and compare sha256 (inputs/ and the ledger);  2. rebuild the frozen-shape file;  3. re-run the three original analysis scripts on the stored raw outputs
# and compare with the stored results;  4. rebuild the CSV tables and compare;  5. regenerate the figures;  6. independent re-computation (reverify.py) and beam re-derivation;  7. SHA256SUMS and privacy scan.
# The simulations themselves (FEBio v4.13 + MKL, Slurm) are NOT re-run.
set -u
ROOT=$(cd "$(dirname "$0")" && pwd); S=${1:-$(mktemp -d)}; mkdir -p "$S"; fail=0; PY="python3 -I"; G=$ROOT/scripts/make_cc_feb.py
echo "scratch: $S"
echo "== 1. regenerate inputs and compare hashes"
for s in "80 8" "160 4" "160 8" "320 8"; do set -- $s; t=S1_$1x$2; $PY $G --mode solid --out $S/gen/$t --tag $t --nsc $1 --ns $2 >/dev/null; done
$PY $G --mode fluid --out $S/gen/F1_unpinned --tag F1_unpinned --nu 40 --nw 320 --nd 80 --nr 16 --k 1e7 --dt 2.5 --t-end 100 --rhoi 0 --shape $ROOT/inputs/shape_S1_320x8.json >/dev/null
$PY $G --mode fsi --out $S/gen/D1_fsi_20steps --tag D1_fsi_20steps --nu 40 --nw 320 --nd 80 --nr 16 --nsc 320 --ns 8 --k 1e7 --dt 0.025 --t-end 0.5 --rhoi 0.9 --e-mesh 1 --plot-stride 5 >/dev/null
$PY $G --mode fluid --out $S/gen/tut_fluid --tag tut_fluid >/dev/null
for t in S1_80x8 S1_160x4 S1_160x8 S1_320x8; do
  if [ "$(sha256sum < $ROOT/inputs/$t.feb)" = "$(sha256sum < $S/gen/$t/$t.feb)" ]; then echo "  $t identical to shipped input"; else echo "  $t MISMATCH"; fail=1; fi; done
for t in S1_80x8 S1_160x4 S1_160x8 S1_320x8 F1_unpinned D1_fsi_20steps; do
  want=$(grep "^collapsibleChannel,$t," $ROOT/provenance/ledger_extract.csv | grep -o '[0-9a-f]\{64\}' | head -1); got=$(sha256sum $S/gen/$t/$t.feb | cut -d' ' -f1)
  if [ "$want" = "$got" ]; then echo "  $t == ledger input_sha256 ${got:0:12}"; else echo "  $t DIFFERS from ledger ($want vs $got)"; fail=1; fi; done
echo "== 2. frozen-shape file"
mkdir -p $S/runs/S1_320x8; zcat $ROOT/raw/S1_320x8/Mon_bottom.txt.gz > $S/runs/S1_320x8/Mon_bottom.txt
$PY $ROOT/scripts/make_shape.py $S/runs/S1_320x8/Mon_bottom.txt $S/shape.json; [ "$(sha256sum < $S/shape.json)" = "$(sha256sum < $ROOT/inputs/shape_S1_320x8.json)" ] && echo "  shape file identical" || { echo "  shape file DIFFERS"; fail=1; }
echo "== 3. original analysis scripts on the stored raw outputs"
for c in S1_80x8 S1_160x4 S1_160x8 S1_320x8; do d=$S/runs/$c; mkdir -p $d; zcat $ROOT/raw/$c/Mon_bottom.txt.gz > $d/Mon_bottom.txt; zcat $ROOT/raw/$c/$c.log.gz > $d/$c.log; done
(cd $S && $PY $ROOT/scripts/cc_s1_analyze.py runs/S1_80x8 runs/S1_160x4 runs/S1_160x8 runs/S1_320x8 >/dev/null)
for c in S1_80x8 S1_160x4 S1_160x8 S1_320x8; do $PY $ROOT/scripts/compare_json.py $ROOT/raw/$c/s1_result.json $S/runs/$c/s1_result.json 1e-9 run || fail=1; done
for c in F1_unpinned D1_fsi_20steps; do d=$S/runs/$c; mkdir -p $d; cp $S/gen/$c/mesh2d.json $S/gen/$c/model_meta.json $d/; for f in Mon_iface Mon_rows Mon_inlet Mon_outlet; do zcat $ROOT/raw/$c/$f.txt.gz > $d/$f.txt; done; zcat $ROOT/raw/$c/$c.log.gz > $d/$c.log; done
$PY $ROOT/scripts/cc_f1_analyze.py $S/runs/F1_unpinned >/dev/null; $PY $ROOT/scripts/compare_json.py $ROOT/raw/F1_unpinned/f1_result.json $S/runs/F1_unpinned/f1_result.json 1e-9 run || fail=1
$PY $ROOT/scripts/compare_json.py $ROOT/raw/F1_unpinned/f1_result_hist.json $S/runs/F1_unpinned/f1_result_hist.json 1e-9 || fail=1
CC_REFERENCE_DIR=$ROOT/reference $PY $ROOT/scripts/cc_fsi_analyze.py $S/runs/D1_fsi_20steps >/dev/null; $PY $ROOT/scripts/compare_json.py $ROOT/raw/D1_fsi_20steps/fsi_result.json $S/runs/D1_fsi_20steps/fsi_result.json 1e-9 run || fail=1
echo "== 4. result tables"
$PY $ROOT/scripts/make_tables.py $ROOT $S/gen $S/tables >/dev/null
for f in $ROOT/results/*.csv; do b=$(basename $f); if cmp -s $f $S/tables/$b; then echo "  $b identical"; else echo "  $b DIFFERS"; fail=1; fi; done
echo "== 5. figures"; $PY $ROOT/scripts/make_figures.py $ROOT/results $ROOT/reference $S/figures >/dev/null 2>&1 && echo "  5 figures regenerated in $S/figures (byte identity with the shipped PNG not required)" || { echo "  figure generation FAILED"; fail=1; }
echo "== 6. independent re-computation and beam re-derivation"
$PY $ROOT/scripts/reverify.py $ROOT $S/gen $S/reverification.json > $S/reverification.txt && tail -1 $S/reverification.txt || { tail -3 $S/reverification.txt; fail=1; }
$PY $ROOT/scripts/compare_json.py $ROOT/results/reverification.json $S/reverification.json 1e-9 || fail=1
if [ "${FULL:-0}" = 1 ]; then NES=10,20,40; else NES=10,20; fi
python3 $ROOT/scripts/beam_ref.py $S/beam.json $NES > /dev/null && python3 - "$ROOT/results/beam_reference_rederivation.json" "$S/beam.json" <<'PY' || fail=1
import json, sys
a, b = json.load(open(sys.argv[1]))["runs"], json.load(open(sys.argv[2]))["runs"]; ok = True
for r in b:
    s = [x for x in a if x["ne"] == r["ne"]][0]; d = abs(s["mid"] - r["mid"]); ok &= d < 1e-8; print(f"  beam ne={r['ne']}: midpoint {r['mid']:.9f} (stored {s['mid']:.9f}) {'ok' if d < 1e-8 else 'DIFFERS'}")
sys.exit(0 if ok else 1)
PY
echo "== 7. SHA256SUMS and privacy scan"
(cd $ROOT && sha256sum -c SHA256SUMS 2>&1 | grep -v ': OK$' | head -5; true)
(cd $ROOT && if { grep -rIil --exclude=reproduce.sh -E '/home/|/Users/|@[A-Za-z0-9.-]+\.(com|org|edu|ac)|xenosim|sairam' . ; for f in $(find . -name '*.gz'); do zcat $f | grep -qiE '/home/|xenosim|sairam' && echo $f; done; } | grep -q .; then echo "  PRIVACY SCAN: matches found"; exit 1; else echo "  privacy scan clean (text and gz: no /home paths, e-mail addresses, user or host names)"; fi) || fail=1
[ $fail -eq 0 ] && echo "ALL NON-SIMULATION CHECKS PASSED" || echo "SOME CHECKS FAILED (see above)"; exit $fail
