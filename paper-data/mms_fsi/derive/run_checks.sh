#!/bin/bash
# Derive, generate and check the manufactured solution.
#   ./run_checks.sh [outdir] [rho_f] [rho_s] [full|flat] [mu_f]
# e.g.  ./run_checks.sh generated            (full variant, as in the paper)
#       ./run_checks.sh generated_flat 1 1 flat
#       ./run_checks.sh generated_mu10 1 1 full 1/10
# Needs python with sympy and symengine (PYTHON, default ../venv/bin/python)
# and a C++17 compiler (CXX, default clang++).
cd "${0%/*}" || exit 1
out=${1:-generated}; rf=${2:-1}; rs=${3:-1}; variant=${4:-full}; muf=${5:-1/100}
PY=${PYTHON:-../venv/bin/python}
"$PY" -u mms_fsi.py --outdir "$out" --rho-f "$rf" --rho-s "$rs" --variant "$variant" --mu-f "$muf" \
    > "$out.derive.log" 2>&1 || { tail -5 "$out.derive.log"; echo "DERIVATION CHECKS FAILED"; }
grep -E "PASS|FAIL|checks passed" "$out.derive.log"
${CXX:-clang++} -O1 -std=c++17 -I "$out" -DREFFILE=\"$out/reference_values.txt\" check_generated.C -o "check_$out" \
    && "./check_$out" | tee "$out.check.log"
