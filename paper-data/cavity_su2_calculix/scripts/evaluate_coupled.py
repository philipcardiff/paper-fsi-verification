#!/usr/bin/env python3
"""Steady cavityFlexibleBottom values of SU2 + CalculiX coupled cases.

    evaluate_coupled.py CASEDIR [CASEDIR ...] [--csv OUT.csv] [--hist-dir DIR] [--tail 20]

* u_x(A), u_y(A): preCICE watch-point PointA (4, 0) on the solid mesh (time = coupled time);
  cross-checked against CalculiX *NODE PRINT of set NA (last increment, both z-faces).
* F_x, F_y on the plate: SU2 history CD/CL (force reference 1, so N per metre depth), converged
  (last) row of each time step; window k = Time_Iter k+1.
* Steady value = final window.  Steadiness: max |change| of each quantity over the last --tail
  windows, relative to the final value; and the first time after which every later window stays
  within 1e-6 (relative) of the final value.
"""
import argparse
import csv
import re
from pathlib import Path


def watchpoint(p):
    L = [l.split() for l in open(p) if l.strip()]
    h = L[0]
    it, ix, iy = h.index("Time"), h.index("Displacement0"), h.index("Displacement1")
    return [(float(r[it]), float(r[ix]), float(r[iy])) for r in L[1:]]


def su2(p):
    R = [{k.strip().strip('"'): v.strip() for k, v in r.items()} for r in csv.DictReader(open(p))]
    last = {}
    for r in R:
        last[int(r["Time_Iter"])] = r
    inner = [int(r["Inner_Iter"]) + 1 for r in R]
    return [(float(r["Cur_Time"]), float(r["CD"]), float(r["CL"])) for _, r in sorted(last.items())], inner


def ccx_na(p):
    """last NODE PRINT block of set NA in plate.dat -> (u_x, u_y) of both nodes."""
    txt = open(p).read()
    blocks = re.findall(r"displacements \(vx,vy,vz\) for set NA and time\s+(\S+)\s*\n\s*\n((?:\s+\d+\s+\S+\s+\S+\s+\S+\s*\n)+)", txt)
    if not blocks:
        return None
    rows = [l.split() for l in blocks[-1][1].strip().split("\n")]
    return float(blocks[-1][0]), [(float(r[1]), float(r[2])) for r in rows]


def settle(t, v, rel=1e-6):
    f = v[-1]
    for i in range(len(v) - 1, -1, -1):
        if abs(v[i] - f) > rel * abs(f):
            return t[i + 1] if i + 1 < len(t) else None
    return t[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--csv")
    ap.add_argument("--hist-dir")
    ap.add_argument("--tail", type=int, default=20)
    a = ap.parse_args()
    table = []
    for c in map(Path, a.cases):
        info = (c / "CASE_INFO.txt").read_text()
        level = int(re.search(r"fluid level (\d+)", info).group(1))
        dt = float(re.search(r"dt ([0-9.e-]+) s", info).group(1))
        wp = watchpoint(c / "solid-calculix/precice-Solid-watchpoint-PointA.log")
        fo, inner = su2(sorted((c / "fluid-su2").glob("history_0*.csv"))[-1])
        its = [int(l.split()[2]) for l in open(c / "solid-calculix/precice-Solid-iterations.log")
               if l.split() and l.split()[0].isdigit()]
        tw = [w[0] for w in wp]
        tf = [f[0] for f in fo]
        q = {"ux": (tw, [w[1] for w in wp]), "uy": (tw, [w[2] for w in wp]),
             "Fx": (tf, [f[1] for f in fo]), "Fy": (tf, [f[2] for f in fo])}
        print(f"\n=== {c}  (level {level}, dt {dt}; solid to t = {tw[-1]:g}, fluid to t = {tf[-1]:g})")
        print(f"{'qty':4} {'final':>15} {'max rel change, last '+str(a.tail):>26} {'within 1e-6 from t':>20}")
        row = {"case": c.name, "level": level, "dt": dt, "t_end": tw[-1]}
        for k, (t, v) in q.items():
            tail = v[-a.tail - 1:]
            ch = max(abs(x - v[-1]) for x in tail) / abs(v[-1])
            s = settle(t, v)
            row.update({k: v[-1], f"{k}_tailchange": ch, f"{k}_settled_t": s})
            print(f"{k:4} {v[-1]:15.8e} {ch:26.2e} {str(s):>20}")
        na = ccx_na(c / "solid-calculix/plate.dat") if (c / "solid-calculix/plate.dat").exists() else None
        if na:
            print(f"CalculiX NA (last print, step time {na[0]:g}): u_y = " + ", ".join(f"{u[1]:.8e}" for u in na[1])
                  + f"  (watch-point {q['uy'][1][-1]:.8e})")
        print(f"windows {len(its)}, coupling its mean {sum(its)/len(its):.2f}, max {max(its)}, at cap (50) {sum(i >= 50 for i in its)}; "
              f"last {a.tail} windows mean {sum(its[-a.tail:])/len(its[-a.tail:]):.2f}")
        print(f"SU2 inner its per solve: mean {sum(inner)/len(inner):.1f}, max {max(inner)}, solves {len(inner)}")
        wall = c / "walltime_coupled.txt"
        if wall.exists():
            print(f"coupled wall {wall.read_text().strip()} s")
            row["wall_coupled_s"] = wall.read_text().strip()
        row.update({"windows": len(its), "coupling_its_mean": sum(its) / len(its), "inner_its_mean": sum(inner) / len(inner)})
        table.append(row)
        if a.hist_dir:
            Path(a.hist_dir).mkdir(parents=True, exist_ok=True)
            F = {round(t, 9): (fx, fy) for t, fx, fy in fo}
            with open(Path(a.hist_dir) / f"{c.name}.csv", "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["time_s", "ux_m", "uy_m", "Fx_N_per_m", "Fy_N_per_m"])
                for t, ux, uy in wp:
                    if round(t, 9) in F:
                        w.writerow([f"{t:.6f}", ux, uy, *F[round(t, 9)]])
    if a.csv and table:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(table[0]))
            w.writeheader()
            w.writerows(table)


if __name__ == "__main__":
    main()
