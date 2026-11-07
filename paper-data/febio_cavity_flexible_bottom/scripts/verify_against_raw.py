#!/usr/bin/env python3
"""Independent re-derivation of the headline numbers from the RAW stored outputs in results/raw (pointA.txt, Mon_iface.txt, Mon_inlet/outlet.txt, logs, facet tables) and comparison with
(a) the stored summaries (fsi_check.json, b3_forces_summary.csv) and (b) the values quoted in docs/*.md.  Package-only (numpy).  Usage from the package root:
    python3 -I scripts/verify_against_raw.py
Writes results/VERIFICATION_LOG.csv.  Exit status 0 if every check is within its stated tolerance.
This checks that numbers quoted in the text equal numbers in the stored outputs; it does NOT verify the FEBio solution (no independent solver is involved)."""
import os, re, io, sys, csv, gzip, json, math, tempfile, subprocess
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cavity_traction as ct
RAW = "results/raw"; checks = []
def chk(name, got, want, tol, kind="abs"):
    err = abs(got - want) if kind == "abs" else abs(got / want - 1)
    ok = bool(err <= tol); checks.append((name, got, want, err, tol, "PASS" if ok else "FAIL")); return ok
def txt(p): return gzip.open(p + ".gz", "rt").read() if not os.path.exists(p) else open(p).read()
def read_cols(path):
    out, t = {}, None
    for ln in txt(path).splitlines():
        if ln.startswith("*Time"): t = float(ln.split("=")[1]); out[t] = {}
        elif not ln.startswith("*") and ln.strip(): f = ln.split(); out[t][int(f[0])] = tuple(map(float, f[1:]))
    return out
def feb_for(tag, meta):
    f = f"inputs/{tag}.feb"
    if os.path.exists(f): return f
    tmp = tempfile.mkdtemp(); args = json.loads([r for r in csv.DictReader(open("inputs/INPUTS.csv")) if r["tag"] == tag][0]["args"])
    subprocess.run([sys.executable, "-I", "scripts/make_cavity_fsi_feb.py", "--out", tmp, "--tag", tag] + args, check=True, capture_output=True); return f"{tmp}/{tag}.feb"
for tag in ("C2_fsi_pilot_f1", "C3_fsi_M2", "C3_fsi_M4", "C4_dt1p25_M1", "C4_K1e7_M1", "C4_tol_M1"):
    rd = f"{RAW}/coupled/{tag}"; fc = json.load(open(f"{rd}/fsi_check.json")); meta = json.load(open(f"{rd}/model_meta.json")); K = meta["k"]
    nodes = ct.read_feb_nodes(feb_for(tag, meta))
    A = read_cols(f"{rd}/pointA.txt"); ts = sorted(A); aid = meta["point_A_node"]
    chk(f"{tag}: u_y(A) final vs pointA.txt", A[ts[-1]][aid][3], fc["uy_A_final"], 1e-12); chk(f"{tag}: t_last", ts[-1], fc["t_last"], 1e-9)
    I = read_cols(f"{rd}/Mon_iface.txt"); tI = max(I); ids = sorted(I[tI], key=lambda i: nodes[i][0])
    X = np.array([nodes[i][0] + I[tI][i][0] for i in ids]); Y = np.array([nodes[i][1] + I[tI][i][1] for i in ids]); P = np.array([-K * I[tI][i][4] for i in ids])
    Fy = float(-np.sum(0.5 * (P[1:] + P[:-1]) * np.diff(X))); chk(f"{tag}: F_y pressure (trapezoid of nodal p, deformed polyline) vs fsi_check", Fy, fc["interface"]["Fy_pressure_raw"], 1e-9)
    chk(f"{tag}: interface mean p", float(P.mean()), fc["interface"]["p_mean"], 1e-9)
    Q = {}
    for nm in ("inlet", "outlet"):
        d = read_cols(f"{rd}/Mon_{nm}.txt"); t = max(d); ii = sorted(d[t], key=lambda i: nodes[i][1]); Q[nm] = float(np.trapz([d[t][i][0] for i in ii], [nodes[i][1] for i in ii]))
    chk(f"{tag}: mass imbalance Qout/Qin-1", Q["outlet"] / Q["inlet"] - 1, fc["mass_imbalance_rel"], 1e-12)
# facet-table reconstruction vs the stored b3 table (same states; never mixing times)
import reconstruct_forces_from_facets as rf
b3 = {r["run_state"]: r for r in csv.DictReader(open(f"{RAW}/pressure_mode/b3_forces_summary.csv"))}
for k in ("M1", "M2", "M3", "M3_t717p5_failure_exit_state"):
    d = np.load(f"results/raw/facets/facets_{k}.npz"); r = rf.recompute(d)
    chk(f"{k}: F1 from facet table", r["F1_Fy"], float(b3[k]["F1_applied_nodal_Fy"]), 5e-7); chk(f"{k}: F2 from facet table", r["F2_Fy"], float(b3[k]["F2_element_mean_Fy"]), 5e-7)
    chk(f"{k}: F4 from facet table (F2=F4 identity for planar facets)", r["F4_Fy"], r["F2_Fy"], 2e-7); chk(f"{k}: D", r["D_row1_minus_row0_Pa"], float(b3[k]["D_row1_minus_row0_Pa"]), 5e-7)
# values quoted in the reports
Q_ = lambda n, v, w, tol, kind="abs": chk("QUOTED " + n, v, w, tol, kind)
fcs = {t: json.load(open(f"{RAW}/coupled/{t}/fsi_check.json")) for t in ("C2_fsi_pilot_f1", "C3_fsi_M2", "C3_fsi_M4")}
Q_("M1 u_y(A) -0.2120 (t=972.5)", fcs["C2_fsi_pilot_f1"]["uy_A_final"], -0.2120, 5e-5); Q_("M2 u_y(A) -0.2262", fcs["C3_fsi_M2"]["uy_A_final"], -0.2262, 5e-5); Q_("M3 u_y(A) -0.2300", fcs["C3_fsi_M4"]["uy_A_final"], -0.2300, 5e-5)
Q_("M1 F_y -4.718", fcs["C2_fsi_pilot_f1"]["interface"]["Fy_pressure_raw"], -4.718, 5e-4); Q_("M2 F_y -4.934", fcs["C3_fsi_M2"]["interface"]["Fy_pressure_raw"], -4.934, 5e-4); Q_("M3 F_y -4.983", fcs["C3_fsi_M4"]["interface"]["Fy_pressure_raw"], -4.983, 5e-4)
Q_("M1 F1 t=900 -4.7176", float(b3["M1"]["F1_applied_nodal_Fy"]), -4.7176, 5e-5); Q_("M1 F2 -5.4375", float(b3["M1"]["F2_element_mean_Fy"]), -5.4375, 5e-5)
Q_("M2 F1 t=750 -4.9339", float(b3["M2"]["F1_applied_nodal_Fy"]), -4.9339, 5e-5); Q_("M2 F2 -5.2748", float(b3["M2"]["F2_element_mean_Fy"]), -5.2748, 5e-5)
Q_("M3 F1 t=500 -4.9809", float(b3["M3"]["F1_applied_nodal_Fy"]), -4.9809, 5e-5); Q_("M3 F2 t=500 -5.2469", float(b3["M3"]["F2_element_mean_Fy"]), -5.2469, 5e-5)
Q_("M3 t=717.5 F2 -5.2494 (separate state)", float(b3["M3_t717p5_failure_exit_state"]["F2_element_mean_Fy"]), -5.2494, 5e-5)
yb = json.load(open(f"{RAW}/b2/ybranch_M3_t717p5.json")); Q_("M3 t=717.5 row parity D = 0.1335 Pa", yb["D_row1_minus_row0_mean_interior"], 0.1335, 1e-4)
Q_("M3 deficit vs QM same dx -5.1 %", fcs["C3_fsi_M4"]["interface"]["Fy_pressure_raw"] / -5.250356 - 1, -0.051, 5e-4); Q_("M1 deficit -7.9 %", fcs["C2_fsi_pilot_f1"]["interface"]["Fy_pressure_raw"] / -5.122626 - 1, -0.079, 5e-4)
Q_("M2 deficit -5.6 %", fcs["C3_fsi_M2"]["interface"]["Fy_pressure_raw"] / -5.226559 - 1, -0.056, 5e-4)
S1 = {k: float(open(f"{RAW}/plate_only/S1_hex8_plates/A{i}.txt").read().split("\n")[-2].split()[2]) for i, k in enumerate(("40x4", "40x8", "80x8", "160x16"))}
Q_("S1 hex8 40x4 -0.23026", S1["40x4"], -0.23026, 5e-6); Q_("S1 hex8 160x16 -0.24023", S1["160x16"], -0.24023, 5e-6)
w = csv.DictReader(open("provenance/slurm_jobs.csv")); cav = [r for r in w if r["group"] in ("cavity-FEBio-run", "cavity-experimental")]
prod = [r for r in cav if r["group"] == "cavity-FEBio-run"]; early_failed = {"13856", "13864", "13875", "13881", "13888"}
Q_("cavity FEBio jobs 13809..14025 excluding the five failed early inlet attempts: 37.9 CPU-h quoted by the campaign status", sum(float(r["cpu_hours_allocated"]) for r in prod if r["job_id"] not in early_failed), 37.9, 0.1)
chk("sacct total of ALL cavity FEBio jobs 13809..14025 (incl. early failed attempts) [stated in INCONSISTENCIES I-1]", sum(float(r["cpu_hours_allocated"]) for r in prod), 38.09, 0.01)
Q_("experimental PPP rounds 1+2 CPU-h (0.68 + 0.046)", sum(float(r["cpu_hours_allocated"]) for r in cav if r["group"] == "cavity-experimental"), 0.727, 0.005)
with open("results/VERIFICATION_LOG.csv", "w", newline="") as fh:
    wr = csv.writer(fh); wr.writerow(["check", "recomputed_or_stored", "reference_or_quoted", "abs_error", "tolerance", "result"])
    for c in checks: wr.writerow([c[0], f"{c[1]:.10g}", f"{c[2]:.10g}", f"{c[3]:.3g}", f"{c[4]:.3g}", c[5]])
bad = [c for c in checks if c[5] == "FAIL"]
for c in checks: print(f"{c[5]}  {c[0]}  got={c[1]:.8g} want={c[2]:.8g} err={c[3]:.2e}")
print(f"{len(checks) - len(bad)}/{len(checks)} checks passed"); sys.exit(1 if bad else 0)
