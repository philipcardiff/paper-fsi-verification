#!/usr/bin/env python3
"""Write the CSV result tables from the stored raw outputs (no simulation).  usage: make_tables.py PACKAGE_ROOT GEN_DIR OUT_DIR
Inputs: raw/*/s1_result.json, f1_result.json, fsi_result.json (outputs of the original analysis scripts), raw/*/Mon_*.txt.gz, reference/*.csv, GEN_DIR/F1_unpinned/mesh2d.json (regenerated)."""
import sys, os, gzip, json, csv, numpy as np
ROOT, GEN, OUT = sys.argv[1:4]; os.makedirs(OUT, exist_ok=True); RAW = ROOT + "/raw"
def w(name, header, rows):
    with open(f"{OUT}/{name}", "w", newline="") as f: c = csv.writer(f); c.writerow(header); c.writerows(rows)
# S1
rows = []
for tag, asp, gate in (("80x8", 20, "3 % (tutorial mesh)"), ("160x4", 5, "3 %"), ("160x8", 5, "informative"), ("320x8", 5, "1 %")):
    r = json.load(open(f"{RAW}/S1_{tag}/s1_result.json")); verdict = {"80x8": "FAIL", "160x4": "FAIL by 0.054 pp", "160x8": "n/a", "320x8": "PASS"}[tag]
    rows.append([tag, asp, r["mid"], round(100 * r["e_mid"], 4), round(100 * r["e_q1"], 4), r["newton_total"], r["newton_avg"], gate, verdict])
w("s1_static_wall.csv", ["mesh_cols_x_rows_hex8", "aspect_ratio", "midpoint_uy_m", "midpoint_error_pct_vs_-0.14362832", "quarter_point_error_pct_vs_-0.09645925", "newton_iterations_total", "newton_per_step_avg", "pre_registered_gate", "gate_result"], rows)
# F1 summary
f = json.load(open(f"{RAW}/F1_unpinned/f1_result.json")); keep = ["K", "n_steps", "Qin", "Qout", "mass_balance", "F1y", "F2y", "rel_F1_F2_y", "G1", "G2", "rel_G", "eps_abs_equiv", "steady_dF1y_last20", "steady_dF2y_last20", "a_t_wall", "a_t_row1", "a_t_row2", "a_t_over_V", "V", "D_n_mean", "D_n_rms", "newton_avg", "junction_p", "neighbour_range_lo", "neighbour_range_hi", "pin_rule_b_lo", "pin_rule_b_hi", "pin_rule_c"]
w("f1_fixed_shape_summary.csv", ["quantity", "value"], [[k, f[k]] for k in keep])
# F1 wall rows along x (final state)
M = json.load(open(f"{GEN}/F1_unpinned/mesh2d.json")); nr = M["nr"]; i0, i1 = M["i0w"], M["i1w"]; xs = np.array(M["xs"]); F = lambda i, j: i * (nr + 1) + j; nodes = np.array(M["nodes2"])
d = {}
for l in gzip.open(f"{RAW}/F1_unpinned/Mon_rows.txt.gz", "rt"):
    if l.startswith("*Time"): t = float(l.split("=")[1]); d[t] = {}
    elif not l.startswith("*") and l.strip(): q = l.split(); d[t][int(q[0])] = float(q[3])
tl = max(d); P = lambda i, j: -1e7 * d[tl][F(i, j) + 1]
w("f1_wall_rows_final.csv", ["x_m", "y_wall_m", "p_wall_Pa", "p_row1_Pa", "p_row2_Pa"], [[nodes[F(i, nr)][0], nodes[F(i, nr)][1], P(i, nr), P(i, nr - 1), P(i, nr - 2)] for i in range(i0, i1 + 1)])
# D1
r = json.load(open(f"{RAW}/D1_fsi_20steps/fsi_result.json")); T = np.array([x["t"] for x in r["rows"]])
ref = [l.split(",") for l in open(f"{ROOT}/reference/collapsibleChannel_oomph_reference.csv") if l[0] not in "#t"]; ref = np.array(ref, float)
out = []
for x in r["rows"]:
    om = float(np.interp(x["t"], ref[:, 0], ref[:, 2]))
    out.append([x["t"], x["F1y"], x["F2y"], x["spread_F"], x["G1"], x["G2"], x["spread_G"], x["a_t"], x["D_mean"], x["D_rms"], x["Qin"], x["Qout"], x["mid"], om, (x["mid"] - om) / 0.216])
w("d1_load_spread_by_step.csv", ["t_s", "F1_y_N_per_m", "F2_y_N_per_m", "rel_F1_minus_F2", "G1", "G2", "rel_G", "streamwise_parity_a_t_Pa", "y_parity_D_mean_Pa", "y_parity_D_rms_Pa", "Q_in", "Q_out", "febio_wall_mid_uy_m", "oomph_converged_wall_mid_uy_m", "difference_over_0.216"], out)
w("d1_freebody_identification.csv", ["t_s", "p_identified_stored_sign_Pa", "p_F1_equiv_Pa", "p_F2_equiv_Pa", "p_odd_row_equiv_Pa"], [[e["t"], e["p_identified"], e["p_F1_equiv"], e["p_F2_equiv"], e["p_odd_equiv"]] for e in r["freebody_ident"]])
w("d1_newton_per_step.csv", ["step", "newton_iterations"], [[i + 1, n] for i, n in enumerate(r["newton_per_step"])])
# oomph reference features (5-point parabola)
def load(n):
    q = [l.split(",") for l in open(f"{ROOT}/reference/collapsibleChannel_oomph_{n}.csv") if l[0] not in "#t"]; a = np.array(q, float); return a[:, 0], a[:, 2]
def ext(t, u, kind, lo, hi):
    m = (t >= lo) & (t <= hi); k = int(np.argmin(np.where(m, u, 9))) if kind == "min" else int(np.argmax(np.where(m, u, -9))); s = slice(max(k - 2, 0), k + 3); c = np.polyfit(t[s], u[s], 2); tt = -c[1] / (2 * c[0]); return tt, float(np.polyval(c, tt))
rows = []
for n in ("reference", "dt40", "dt80", "dt160", "dt320"):
    t, u = load(n)
    feats = [("first_trough", "min", 0.3, 0.75), ("first_rebound_max", "max", 0.6, 1.0), ("second_trough", "min", 0.9, 1.4), ("third_trough", "min", 1.5, 2.0)]
    ff = {nm: ext(t, u, k, lo, hi) for nm, k, lo, hi in feats if t[-1] > hi}
    for nm, (tt, uu) in ff.items(): rows.append([n, nm, round(tt, 4), round(uu, 5)])
    if "second_trough" in ff: rows.append([n, "period_trough1_to_trough2_s", round(ff["second_trough"][0] - ff["first_trough"][0], 4), ""])
    if "third_trough" in ff: rows.append([n, "period_trough2_to_trough3_s", round(ff["third_trough"][0] - ff["second_trough"][0], 4), ""])
    if "first_rebound_max" in ff: rows.append([n, "spacing_trough1_to_rebound_max_s", round(ff["first_rebound_max"][0] - ff["first_trough"][0], 4), ""])
w("oomph_reference_features.csv", ["csv_file_collapsibleChannel_oomph_<name>", "feature", "time_s", "wall_mid_uy_m"], rows)
