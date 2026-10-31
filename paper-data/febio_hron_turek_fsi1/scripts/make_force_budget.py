#!/usr/bin/env python3
"""results/force_budget_<run>.csv from results/independent_forces_<run>.json: term-by-term drag (x) and lift (y) budget of the reconstructions F1, F2, F3 (nodal), F3b, F4 and the
F3b - F2 difference.  usage: make_force_budget.py RESULTS_DIR RUN_TAG"""
import sys, json, csv
R, tag = sys.argv[1], sys.argv[2]; d = json.load(open(f"{R}/independent_forces_{tag}.json")); rows = []
def add(item, x, y, note=""): rows.append(dict(item=item, drag_x=x, lift_y=y, note=note))
for m in ("F1", "F2"):
    for g, nm in (("cyl", "cylinder"), ("flag", "flag")):
        p, v = d[m + "_parts"][g]["p"], d[m + "_parts"][g]["v"]; add(f"{m} {nm} pressure", p[0], p[1]); add(f"{m} {nm} viscous", v[0], v[1])
    add(f"{m} total", d[m][0], d[m][1])
for k, v in d["F3b_terms"].items(): add("F3b term " + k, v[0], v[1], "control-volume term (outer boundary, element-mean pressure); momentum terms enter the sum with a minus sign")
add("F3b total", d["F3b"][0], d["F3b"][1]); add("F3 (nodal pressure on outer boundary) total", d["F3nodal"][0], d["F3nodal"][1], "parity-contaminated")
for k, v in d["F4_parts"].items(): add("F4 part " + k, v[0], v[1], "weak-form nodal reactions, obstacle nodes")
add("F4 total", d["F4_weak_reaction"][0], d["F4_weak_reaction"][1])
w = d["weak_identity"]
for k in ("sum_obstacle_f", "sum_outer_f", "sum_interior_f", "convective_volume_integral", "closure_sum_all_minus_conv"): add("weak identity " + k, w[k][0], w[k][1])
F2, F3b, F4 = d["F2"], d["F3b"], d["F4_weak_reaction"]
t = d["F3b_terms"]; sig_facet = [sum(t[k][i] for k in t if k[0] in "pv") for i in (0, 1)]; mom = [-(t["mom_inlet"][i] + t["mom_outlet"][i]) for i in (0, 1)]
add("DIFF F3b - F2", F3b[0] - F2[0], F3b[1] - F2[1], "total difference to be explained")
add("DIFF F4 - F2", F4[0] - F2[0], F4[1] - F2[1], "obstacle-side reconstruction error of F2 relative to the weak-form force")
add("DIFF F3b - F4", F3b[0] - F4[0], F3b[1] - F4[1], "outer-boundary reconstruction error of F3b relative to the weak-form force")
add("  of which facet sigma.n on outer boundary minus weak outer nodal sum", sig_facet[0] - w["sum_outer_f"][0], sig_facet[1] - w["sum_outer_f"][1], "pressure on the outer boundary (element mean vs nodal weak form) and wall shear")
add("  of which -(boundary momentum flux - convective volume integral)", mom[0] - (-w["convective_volume_integral"][0]) * 1 if False else (mom[0] + w["convective_volume_integral"][0]) * 1, (mom[1] + w["convective_volume_integral"][1]) * 1, "momentum flux through inlet/outlet vs volume integral of rho (w.grad)w: the difference is rho w div(w) (discrete compressibility/continuity term)")
add("  of which -interior residual sum", -w["sum_interior_f"][0], -w["sum_interior_f"][1], "")
with open(f"{R}/force_budget_{tag}.csv", "w", newline="") as fh:
    wr = csv.DictWriter(fh, ["item", "drag_x", "lift_y", "note"]); wr.writeheader(); wr.writerows(rows)
for r in rows: print(f"{r['item']:75s} {r['drag_x']:12.5f} {r['lift_y']:12.5f}")
