#!/usr/bin/env python3
"""Steady values of the Kratos beamInCrossFlow runs and their comparison with solids4foam and Richter (2012).

Reads each run's history.csv, run_info.json, run_meta.txt and time.txt, and the solids4foam sequences
in paper-data/beam_richter/ (r02_sequence.csv for the 0.2 m/s peak, production_sequence.csv for the
printed 0.3 m/s peak). Writes a Markdown summary and a CSV of the steady values.

The steady value of a run is its last row (the run stops when u_x(A) and Fc_x change by less than 1e-6,
relative, over 1 s). The repeat check compares every history column except the step wall time.

Usage (from this folder):
  python scripts/beam_qoi.py --runs runs --s4f ../beam_richter --static checks/solid_static \
      --out results/summary.md --csv results/steady_values.csv
"""
import argparse
import csv
import json
import math
import os

RICHTER = {"ux": 5.95e-5, "Fx": 1.33}          # Richter (2012), peak 0.2 m/s
KLEV = (0, 1, 2)                                # Kratos levels
RUNS = {(0, 0.2): "f1_U0.2_dt0.05", (1, 0.2): "f2_U0.2_dt0.05", (2, 0.2): "f4_U0.2_dt0.05_amg",
        (0, 0.3): "f1_U0.3_dt0.05", (1, 0.3): "f2_U0.3_dt0.05", (2, 0.3): "f4_U0.3_dt0.05_amg"}
MAIN = [RUNS[(lev, U)] for U in (0.2, 0.3) for lev in KLEV]
DT_CHECK = ("f1_U0.2_dt0.1", "f1_U0.2_dt0.05")
SOLVER_CHECK = ("f1_U0.2_dt0.05_amg", "f1_U0.2_dt0.05")   # AMG plate solve against the direct one
REPEATS = ["f1_U0.2_dt0.05", "f1_U0.2_dt0.1", "f1_U0.3_dt0.05"]
LEVEL = {1: 0, 2: 1, 4: 2, 8: 3}                # factor -> solids4foam level


def read_rows(run):
    return list(csv.DictReader(open(os.path.join(run, "history.csv"))))


def meta(run):
    out = {}
    for line in open(os.path.join(run, "run_meta.txt")):
        k, _, v = line.strip().partition("=")
        out[k] = v
    return out


def peak_memory_gb(run):
    for line in open(os.path.join(run, "time.txt")):
        if "Maximum resident set size" in line:
            return int(line.split(":")[-1]) / 1024 ** 2
    return float("nan")


def amgcl_warnings(run):
    """AMGCL 'Non converged linear solution' warnings in log.txt, as (time, final residual) pairs."""
    out, t = [], float("nan")
    for line in open(os.path.join(run, "log.txt"), errors="replace"):
        if "TIME =" in line:
            t = float(line.split("=")[1])
        elif "Non converged linear solution" in line:
            out.append((t, float(line.split("[")[-1].split(">")[0])))
    return out


def steady(run):
    rows = read_rows(run)
    warn = amgcl_warnings(run)
    case = json.load(open(os.path.join(run, "case.json")))
    info = json.load(open(os.path.join(run, "run_info.json")))
    last = rows[-1]
    t_end = float(last["time"])
    old = min(rows, key=lambda r: abs(float(r["time"]) - (t_end - 1.0)))
    v = {k: float(last[k]) for k in last if k not in ("converged",)}
    rel = lambda k: abs(float(last[k]) - float(old[k])) / abs(float(last[k]))
    solves = [int(r["coupling_solves"]) for r in rows]
    return {
        "run": os.path.basename(run.rstrip("/")), "factor": case["factor"], "level": LEVEL[case["factor"]],
        "U": case["inlet_peak"], "dt": case["dt"], "fluid_tets": case["fluid_tets"],
        "fluid_nodes": case["fluid_nodes"], "solid_hexes": case["solid_hexes"],
        "t_stop": t_end, "steps": len(rows), "coupling_solves": sum(solves),
        "solves_per_step": sum(solves) / len(rows), "all_converged": all(r["converged"] == "1" for r in rows),
        "max_rms_residual": max(float(r["rms_residual_at_acceptance"]) for r in rows),
        "change_ux_1s": rel("ux_A"), "change_Fcx_1s": rel("Fcx"),
        "ux": v["ux_A"], "uy": v["uy_A"], "uz": v["uz_A"],
        "Fcx": v["Fcx"], "Fcy": v["Fcy"], "Fcz": v["Fcz"], "Ffx": v["Ffx"], "Ffy": v["Ffy"], "Ffz": v["Ffz"],
        "Frx": v["Frx"], "Fry": v["Fry"], "Frz": v["Frz"],
        "threads": info["kratos_threads"], "wall_h": info["wall_total_s"] / 3600.0,
        "peak_memory_GB": peak_memory_gb(run), "exit": meta(run).get("exit"), "stop_reason": info["stop_reason"],
        "amgcl_warnings": len(warn), "amgcl_max_residual": max((w for _, w in warn), default=0.0),
    }


def s4f_sequences(d):
    r02 = {}
    for r in csv.DictReader(open(os.path.join(d, "r02_sequence.csv"))):
        r02[int(r["level"][1:])] = {"ux": float(r["ux"]), "uy": float(r["uy"]), "uz": float(r["uz"]),
                                     "Fx": float(r["Fx"]), "Fy": float(r["Fy"]), "Fz": float(r["Fz"])}
    r03 = {}
    for r in csv.DictReader(open(os.path.join(d, "production_sequence.csv"))):
        r03[int(r["level"])] = {"ux": float(r["ux"]), "uy": float(r["uy"]), "uz": float(r["uz"]),
                                "Fx": float(r["fx"]), "Fy": float(r["fy"]), "Fz": float(r["fz"])}
    return {0.2: r02, 0.3: r03}


def richardson(f1, f2, f3):
    """Extrapolated value and observed order from three levels at refinement ratio 2."""
    p = math.log(abs((f1 - f2) / (f2 - f3))) / math.log(2.0)
    return f3 + (f3 - f2) / (2.0 ** p - 1.0), p


def conditional(f1, f2, f3):
    """richardson() for a monotone sequence with decreasing changes, else None."""
    return richardson(f1, f2, f3) if (f1 - f2) / (f2 - f3) > 1.0 else None


def pct(a, b):
    return f"{100.0 * (a / b - 1.0):+.2f}%"


def repeat_check(runs_dir, name):
    a = read_rows(os.path.join(runs_dir, "_first_pass", name))
    b = read_rows(os.path.join(runs_dir, name))
    if len(a) != len(b):
        return f"different lengths ({len(a)}, {len(b)} steps)"
    for ra, rb in zip(a, b):
        for k in ra:
            if k != "wall_step_s" and ra[k] != rb[k]:
                return f"differs at t = {ra['time']} in {k}"
    return f"identical, {len(a)} steps"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs")
    ap.add_argument("--s4f", default="../beam_richter")
    ap.add_argument("--static", default="checks/solid_static")
    ap.add_argument("--out", default="results/summary.md")
    ap.add_argument("--csv", default="results/steady_values.csv")
    a = ap.parse_args()

    names = MAIN + [DT_CHECK[0], SOLVER_CHECK[0]]
    res = {n: steady(os.path.join(a.runs, n)) for n in names}
    s4f = s4f_sequences(a.s4f)

    keys = list(res[MAIN[0]].keys())
    with open(a.csv, "w") as fo:
        w = csv.DictWriter(fo, fieldnames=keys)
        w.writeheader()
        for n in names:
            w.writerow({k: (f"{v:.10g}" if isinstance(v, float) else v) for k, v in res[n].items()})

    L = ["# beamInCrossFlow, Kratos 10.3.0, steady values\n",
         "Generated by `scripts/beam_qoi.py`. Steady value = last step of each run. Fc: from the interface "
         "reactions minus the floor share at the plate base edge (primary). Ff: direct integral of the P1 "
         "traction over the wetted plate faces (cross-check). Forces on the half plate.\n"]

    L += ["## Runs\n",
          "| Run | Level | U (m/s) | Δt (s) | Steps | Stop (s) | Change over the last 1 s, u_x / Fc_x | "
          "Coupling solves per step | All steps converged | Max RMS residual (m) | Fluid AMGCL warnings (largest) | "
          "Threads | Wall (h) | Peak memory (GB) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for n in names:
        r = res[n]
        warn = f"{r['amgcl_warnings']} ({r['amgcl_max_residual']:.1e})" if r["amgcl_warnings"] else "0"
        L.append(f"| `{n}` | {r['level']} | {r['U']} | {r['dt']} | {r['steps']} | {r['t_stop']:.6g} | "
                 f"{r['change_ux_1s']:.1e} / {r['change_Fcx_1s']:.1e} | {r['solves_per_step']:.2f} | "
                 f"{'yes' if r['all_converged'] else 'NO'} | {r['max_rms_residual']:.1e} | {warn} | {r['threads']} | "
                 f"{r['wall_h']:.2f} | {r['peak_memory_GB']:.2f} |")

    L += ["\n## Steady values\n",
          "| Run | u_x(A) (m) | u_y(A) (m) | u_z(A) (m) | Fc_x (N) | Ff_x (N) | Fc_x − Ff_x (N) | Fc_y (N) | Ff_y (N) | Fc_z (N) | Ff_z (N) |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for n in names:
        r = res[n]
        L.append(f"| `{n}` | {r['ux']:.6e} | {r['uy']:.5e} | {r['uz']:.5e} | {r['Fcx']:.6f} | {r['Ffx']:.6f} | "
                 f"{r['Fcx'] - r['Ffx']:.4f} | {r['Fcy']:.5f} | {r['Ffy']:.5f} | {r['Fcz']:.5f} | {r['Ffz']:.5f} |")

    L += ["\n## Against solids4foam and Richter (2012)\n",
          "solids4foam levels 0 to 3 from `paper-data/beam_richter/` (Table `tab:beam_seq`). Kratos levels 0 "
          "to 2 use the same fluid and plate block counts as solids4foam levels 0 to 2. A conditional Richardson "
          "value is given only for a monotone sequence with decreasing changes.\n"]
    for U in (0.2, 0.3):
        k = [res[RUNS[(lev, U)]] for lev in KLEV]
        s = s4f[U]
        ext_u, p_u = richardson(s[1]["ux"], s[2]["ux"], s[3]["ux"])
        ext_f, p_f = richardson(s[1]["Fx"], s[2]["Fx"], s[3]["Fx"])
        kx = {q: conditional(*(r[q] for r in k)) for q in ("ux", "Fcx", "Ffx")}
        fmt = lambda e, scale=1.0, d=4: "not monotone" if e is None else f"{scale * e[0]:.{d}f} (order {e[1]:.2f})"
        L += [f"### Peak {U} m/s\n",
              "| | u_x(A) (1e-5 m) | F_x (N) | u_y(A) (1e-7 m) | u_z(A) (1e-7 m) | F_y (N) | F_z (N) |",
              "|---|---|---|---|---|---|---|"]
        for lev, r in zip(KLEV, k):
            L.append(f"| Kratos {lev} | {1e5 * r['ux']:.4f} | Fc {r['Fcx']:.4f}, Ff {r['Ffx']:.4f} | "
                     f"{1e7 * r['uy']:.4f} | {1e7 * r['uz']:.4f} | Fc {r['Fcy']:.4f}, Ff {r['Ffy']:.4f} | "
                     f"Fc {r['Fcz']:.4f}, Ff {r['Ffz']:.4f} |")
        L.append(f"| Kratos, conditional Richardson from levels 0–2 | {fmt(kx['ux'], 1e5)} | "
                 f"Fc {fmt(kx['Fcx'])}, Ff {fmt(kx['Ffx'])} | | | | |")
        for lev in range(4):
            q = s[lev]
            L.append(f"| solids4foam {lev} | {1e5 * q['ux']:.4f} | {q['Fx']:.4f} | {1e7 * q['uy']:.4f} | "
                     f"{1e7 * q['uz']:.4f} | {q['Fy']:.4f} | {q['Fz']:.4f} |")
        L.append(f"| solids4foam, conditional Richardson from levels 1–3 | {1e5 * ext_u:.4f} (order {p_u:.2f}) | "
                 f"{ext_f:.4f} (order {p_f:.2f}) | | | | |")
        if U == 0.2:
            L.append(f"| Richter (2012) | {1e5 * RICHTER['ux']:.2f} | {RICHTER['Fx']:.2f} | | | | |")
        q3 = s[3]
        L.append("")
        for lev in (1, 2):
            r = k[lev]
            L.append(f"- Kratos level {lev} against solids4foam level 3: u_x(A) {pct(r['ux'], q3['ux'])}, "
                     f"Fc_x {pct(r['Fcx'], q3['Fx'])}, Ff_x {pct(r['Ffx'], q3['Fx'])}, "
                     f"u_y(A) {pct(r['uy'], q3['uy'])}, u_z(A) {pct(r['uz'], q3['uz'])}, "
                     f"Fc_y {pct(r['Fcy'], q3['Fy'])}, Ff_z {pct(r['Ffz'], q3['Fz'])}.")
        L.append(f"- Against the solids4foam conditional Richardson values: Kratos level 1 u_x(A) "
                 f"{pct(k[1]['ux'], ext_u)}, level 2 u_x(A) {pct(k[2]['ux'], ext_u)}, Fc_x {pct(k[2]['Fcx'], ext_f)}, "
                 f"Ff_x {pct(k[2]['Ffx'], ext_f)}"
                 + (f"; the Kratos conditional values u_x(A) {pct(kx['ux'][0], ext_u)}" if kx["ux"] else "")
                 + (f", Fc_x {pct(kx['Fcx'][0], ext_f)}" if kx["Fcx"] else "") + ".")
        L.append(f"- Kratos level 0 to 1 and 1 to 2: u_x(A) {pct(k[1]['ux'], k[0]['ux'])}, {pct(k[2]['ux'], k[1]['ux'])}; "
                 f"Fc_x {pct(k[1]['Fcx'], k[0]['Fcx'])}, {pct(k[2]['Fcx'], k[1]['Fcx'])}; "
                 f"Ff_x {pct(k[1]['Ffx'], k[0]['Ffx'])}, {pct(k[2]['Ffx'], k[1]['Ffx'])}; Fc_x − Ff_x "
                 + " → ".join(f"{r['Fcx'] - r['Ffx']:.4f}" for r in k) + " N.")
        L.append(f"- solids4foam level 0 to 1: u_x(A) {pct(s[1]['ux'], s[0]['ux'])}, F_x {pct(s[1]['Fx'], s[0]['Fx'])}.")
        if U == 0.2:
            for lev in (1, 2):
                r = k[lev]
                L.append(f"- Kratos level {lev} against Richter (2012): u_x(A) {pct(r['ux'], RICHTER['ux'])}, "
                         f"Fc_x {pct(r['Fcx'], RICHTER['Fx'])}, Ff_x {pct(r['Ffx'], RICHTER['Fx'])}.")
        L.append("")

    L += ["### Ratio of the two inflows (0.3 / 0.2)\n",
          "| | u_x(A) | F_x |", "|---|---|---|"]
    for lev in KLEV:
        a2, a3 = res[RUNS[(lev, 0.2)]], res[RUNS[(lev, 0.3)]]
        L.append(f"| Kratos {lev} | {a3['ux'] / a2['ux']:.4f} | Fc {a3['Fcx'] / a2['Fcx']:.4f}, "
                 f"Ff {a3['Ffx'] / a2['Ffx']:.4f} |")
    for lev in range(4):
        L.append(f"| solids4foam {lev} | {s4f[0.3][lev]['ux'] / s4f[0.2][lev]['ux']:.4f} | "
                 f"{s4f[0.3][lev]['Fx'] / s4f[0.2][lev]['Fx']:.4f} |")

    L += ["\n## Checks\n"]
    d1, d2 = res[DT_CHECK[0]], res[DT_CHECK[1]]
    L.append(f"- **Time step.** Level 0, peak 0.2: Δt 0.1 against 0.05 changes u_x(A) by "
             f"{abs(d1['ux'] / d2['ux'] - 1):.1e}, Fc_x by {abs(d1['Fcx'] / d2['Fcx'] - 1):.1e} and Ff_x by "
             f"{abs(d1['Ffx'] / d2['Ffx'] - 1):.1e} (relative).")
    ra = read_rows(os.path.join(a.runs, SOLVER_CHECK[0]))
    rb = read_rows(os.path.join(a.runs, SOLVER_CHECK[1]))
    m1, m2 = res[SOLVER_CHECK[0]], res[SOLVER_CHECK[1]]
    du = max(abs(float(x["ux_A"]) - float(y["ux_A"])) for x, y in zip(ra, rb))
    same = sum(x["coupling_solves"] == y["coupling_solves"] for x, y in zip(ra, rb))
    L.append(f"- **Plate solver.** Level 0, peak 0.2: CG with AMG (`{SOLVER_CHECK[0]}`) against the direct solve "
             f"(`{SOLVER_CHECK[1]}`), {len(ra)} and {len(rb)} steps. u_x(A) differs by at most {du:.1e} m over the "
             f"run; the steady u_x(A), Fc_x and Ff_x differ by {abs(m1['ux'] / m2['ux'] - 1):.1e}, "
             f"{abs(m1['Fcx'] / m2['Fcx'] - 1):.1e} and {abs(m1['Ffx'] / m2['Ffx'] - 1):.1e} (relative); the "
             f"coupling solves are the same in {same} of {min(len(ra), len(rb))} steps.")
    for n in names:
        w = amgcl_warnings(os.path.join(a.runs, n))
        if w:
            big = sorted({t for t, v in w if v > 1e-8})
            L.append(f"- **Fluid linear solver** `{n}`: {len(w)} AMGCL warnings (tolerance 1e-10) in "
                     f"{len({t for t, _ in w})} steps, final residuals {min(v for _, v in w):.1e} to "
                     f"{max(v for _, v in w):.1e}; above 1e-8 only at t = "
                     + ", ".join(f"{t:.4g} s" for t in big) + ". The other runs have none.")
    for n in REPEATS:
        L.append(f"- **Repeat** `_first_pass/{n}` against `{n}`: {repeat_check(a.runs, n)} "
                 "(all columns except the step wall time).")
    st = json.load(open(os.path.join(a.static, "solid_static_TotalLagrangianElement3D8N.json")))
    L += ["\n### Plate alone, 100 Pa in +x on the upstream face (`scripts/solid_static_check.py`)\n",
          "| Factor | Plate cells | Load (N) | Kratos u_x(A) (m) | solids4foam m2 u_x(A) (m) | Kratos / solids4foam − 1 |",
          "|---|---|---|---|---|---|"]
    for r in st:
        L.append(f"| {r['factor']} | {' × '.join(str(c) for c in r['cells_xyz'])} | {r['load_N']:.5f} | "
                 f"{r['ux_A']:.5e} | {r['s4f_m2_ux_A']:.5e} | {pct(r['ux_A'], r['s4f_m2_ux_A'])} |")

    text = "\n".join(L) + "\n"
    with open(a.out, "w") as fo:
        fo.write(text)
    print(text)


if __name__ == "__main__":
    main()
