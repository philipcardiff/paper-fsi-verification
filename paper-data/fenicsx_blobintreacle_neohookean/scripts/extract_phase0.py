#!/usr/bin/env python3
"""Regenerate data/solid_only_apex.csv and data/law_checks.csv from the Phase 0 records in the evidence tarball.

usage: python3 -B -I scripts/extract_phase0.py EVIDENCE_DIR [OUT_DIR]
Reads EVIDENCE_DIR/benchmarks/blobInTreacle/nh/phase0/{law_checks.json,s4f_transcript_comparison.json,solid_only/...}."""
import csv, json, os, sys

ev = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
P0 = os.path.join(ev, "benchmarks", "blobInTreacle", "nh", "phase0")
SOLID_CELLS = {"0.5": 30, "1": 120, "2": 480, "4": 1920}
with open(os.path.join(out, "solid_only_apex.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["load", "level", "solid_cells", "code", "ux", "uy", "closure_applied_plus_reaction_max"])
    for ld in ("L1", "L2"):
        for lv in SOLID_CELLS:
            j = json.load(open(os.path.join(P0, "solid_only", "fe_f%s_%s.json" % (lv, ld))))
            clo = max(abs(x) for x in j["closure_applied_plus_reaction"])
            w.writerow([ld, lv, SOLID_CELLS[lv], "FEniCSx", "%.10g" % j["apex_ux"], "%.10g" % j["apex_uy"], "%.3e" % clo])
            dat = os.path.join(P0, "solid_only", "s4f", "s4f_solid_x%s_%s" % (lv, ld), "solidPointDisplacement_pointDisp.dat")
            row = [l.split() for l in open(dat) if not l.startswith("#")][-1]
            w.writerow([ld, lv, SOLID_CELLS[lv], "solids4foam", row[1], row[2], ""])
law = json.load(open(os.path.join(P0, "law_checks.json"))); tr = json.load(open(os.path.join(P0, "s4f_transcript_comparison.json")))
rows = [("P: UFL-differentiated vs numpy (max rel., 20 F)", law["max_relative_errors_over_F_set"]["P_ufl_vs_np"], "law_checks.json"),
        ("P: analytic vs numpy", law["max_relative_errors_over_F_set"]["P_analytic_vs_np"], "law_checks.json"),
        ("sigma: UFL vs numpy 3-D Cauchy", law["max_relative_errors_over_F_set"]["sigma_ufl_vs_np3"], "law_checks.json"),
        ("sigma = P F^T/J vs P from Cauchy", law["max_relative_errors_over_F_set"]["P_from_cauchy"], "law_checks.json"),
        ("S = F^-1 P symmetry", law["max_relative_errors_over_F_set"]["S_sym"], "law_checks.json"),
        ("tangent: UFL derivative vs central finite difference", law["max_relative_errors_over_F_set"]["tangent_ufl_vs_fd"], "law_checks.json"),
        ("objectivity (rotation) of P", law["objectivity_rotation_1.234rad"]["P"], "law_checks.json"),
        ("objectivity (rotation) of sigma", law["objectivity_rotation_1.234rad"]["sigma"], "law_checks.json"),
        ("small-strain limit, relative error at eps=1e-5", law["small_strain_vs_linear_elastic_lambda_mu"][-1]["rel_err"], "law_checks.json"),
        ("s4f formulas (C++ transcription) vs module: sigma, s4f constants", tr["max_rel"]["sigma_with_s4f_digits_mu_K"], "s4f_transcript_comparison.json"),
        ("s4f formulas (C++ transcription) vs module: P, s4f constants", tr["max_rel"]["P_with_s4f_digits"], "s4f_transcript_comparison.json"),
        ("s4f formulas vs module: sigma, exact (50, 533.33) constants", tr["max_rel"]["sigma_with_exact_mu_K"], "s4f_transcript_comparison.json"),
        ("s4f tutorial E,nu rounding: relative deviation of mu", law["s4f_tutorial_inputs_derived"]["rel_dev_mu"], "law_checks.json"),
        ("s4f tutorial E,nu rounding: relative deviation of K", abs(law["s4f_tutorial_inputs_derived"]["rel_dev_K"]), "law_checks.json")]
with open(os.path.join(out, "law_checks.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["check", "value", "source"])
    for a, b, c in rows: w.writerow([a, "%.3e" % b, c])
    # hand-transcribed from the internal audit text (an separate re-implementation, not regenerated here)
    w.writerow(["internal audit, separate 3x3 re-implementation (F_zz=1): Cauchy stress vs module", "4e-16", "REPORTED in the internal audit report (evidence archive, NH_INDEPENDENT_AUDIT.md); not regenerated"])
    w.writerow(["internal audit, separate re-implementation: tangent vs finite difference", "1.3e-09", "REPORTED in the internal audit report (evidence archive, NH_INDEPENDENT_AUDIT.md); not regenerated"])
print("wrote", out)
