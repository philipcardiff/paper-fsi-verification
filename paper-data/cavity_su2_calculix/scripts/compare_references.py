#!/usr/bin/env python3
"""Compare the SU2 + CalculiX steady cavity values with the other cavityFlexibleBottom solutions in this repository.

Our values were fixed and posted blind on #36 (issuecomment-6928621052) before any of these references was opened.

Sources (all relative to this package):
  cavity_su2_steady.csv                                   this solution, levels 1-16 (dx = 0.1 / L m)
  comparison/s4f_partitioned_issue36.csv                  solids4foam partitioned St Venant-Kirchhoff, strict coupling
                                                          test, dt halved per level, as posted by Philip Cardiff on #36
                                                          (issuecomment-7015729752)
  ../cavity_monolithic_618/quasi_monolithic_618.csv       solids4foam quasi-monolithic, #618 build (and #597 columns)
  ../fenicsx_cavityflexiblebottom/data/coupled_levels.csv FEniCSx monolithic (same research group)
  ../cavity_plate_material_law/reference_levels.csv       Tukovic et al. 2018, ANSYS, COMSOL (paper's reference levels)

The limit of each sequence is the Richardson value of its three finest levels (ours: L4-L8-L16), with the observed
order; it is not given where the triple is not monotone.  Differences are (reference / ours - 1) in percent.

Writes comparison/reference_comparison.csv (levels) and comparison/reference_limits.csv (limits) and prints both.
Usage (from the package folder):  python3 scripts/compare_references.py
"""
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PD = ROOT.parent


def rows(path):
    return list(csv.DictReader(open(path, encoding="utf-8-sig")))


def richardson(v):
    """Observed order and limit of the three finest values (coarse to fine), or (None, None) if not monotone."""
    a, b, c = v[-3:]
    r = (a - b) / (b - c) if b != c else -1
    if r <= 0:
        return None, None
    p = math.log2(r)
    return p, c + (c - b) / (2 ** p - 1)


def main():
    seqs = {}   # name -> {dx: (uy, Fy)}
    ours = {0.1 / int(r["level"]): (float(r["uy"]), float(r["Fy"])) for r in rows(ROOT / "cavity_su2_steady.csv")}
    seqs["SU2 + CalculiX (this work)"] = ours
    seqs["solids4foam partitioned (strict, #36)"] = {float(r["dx_m"]): (float(r["uy_A_m"]), float(r["Fy_N_per_m"]))
                                                     for r in rows(ROOT / "comparison/s4f_partitioned_issue36.csv")}
    m = rows(PD / "cavity_monolithic_618/quasi_monolithic_618.csv")
    seqs["solids4foam quasi-monolithic #618"] = {float(r["delta_x_m"]): (float(r["uy_top_m"]), float(r["fy_N_per_m"])) for r in m}
    seqs["solids4foam quasi-monolithic #597"] = {float(r["delta_x_m"]): (float(r["uy_top_597_m"]), float(r["fy_597_N_per_m"])) for r in m}
    seqs["FEniCSx monolithic"] = {float(r["dx_m"]): (float(r["uy_A_m"]), float(r["Fy_N_per_m"]))
                                  for r in rows(PD / "fenicsx_cavityflexiblebottom/data/coupled_levels.csv")}
    for r in rows(PD / "cavity_plate_material_law/reference_levels.csv"):
        if r["source"] != "solids4foam":   # the paper's rounded solids4foam levels: superseded by the #36 sequence
            seqs.setdefault(r["source"], {})[float(r["dx_m"])] = (float(r["uy_A_m"]), float(r["Fy_N_per_m"]))

    dxs = [0.1, 0.05, 0.025, 0.0125]
    lev = []
    for name, s in seqs.items():
        for dx in dxs:
            if dx in s and name != "SU2 + CalculiX (this work)":
                uy, fy = s[dx]
                lev.append(dict(source=name, dx_m=dx, uy_A_m=uy, Fy_N_per_m=fy, our_uy_A_m=ours[dx][0], our_Fy_N_per_m=ours[dx][1],
                                uy_diff_pct=100 * (uy / ours[dx][0] - 1), Fy_diff_pct=100 * (fy / ours[dx][1] - 1)))
    lim = []
    po, uo = richardson([ours[d][0] for d in sorted(ours, reverse=True)])
    qo, fo = richardson([ours[d][1] for d in sorted(ours, reverse=True)])
    for name, s in seqs.items():
        d = sorted(s, reverse=True)
        pu, lu = richardson([s[x][0] for x in d])
        pf, lf = richardson([s[x][1] for x in d])
        fin = d[-1]
        lim.append(dict(source=name, finest_dx_m=fin, finest_uy_A_m=s[fin][0], finest_Fy_N_per_m=s[fin][1],
                        uy_order=pu, uy_limit=lu, Fy_order=pf, Fy_limit=lf,
                        uy_limit_diff_pct=None if lu is None else 100 * (lu / uo - 1),
                        Fy_limit_diff_pct=None if lf is None else 100 * (lf / fo - 1),
                        finest_uy_diff_pct=100 * (s[fin][0] / uo - 1), finest_Fy_diff_pct=100 * (s[fin][1] / fo - 1)))

    def fmt(v):
        return "" if v is None else (f"{v:.7g}" if isinstance(v, float) else v)
    for path, data in (("comparison/reference_comparison.csv", lev), ("comparison/reference_limits.csv", lim)):
        with open(ROOT / path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(list(data[0]))
            w.writerows([[fmt(r[k]) for k in data[0]] for r in data])
        print(f"wrote {path}")
    print("\nlimits (differences against our limit; finest = against our limit too)")
    for r in lim:
        o = lambda v: "   -  " if v is None else f"{v:5.2f}"
        l = lambda v: "      -     " if v is None else f"{v:+.6f}"
        d = lambda v: "     -  " if v is None else f"{v:+7.3f}%"
        print(f"{r['source']:40s} u_y p {o(r['uy_order'])} lim {l(r['uy_limit'])} {d(r['uy_limit_diff_pct'])} | "
              f"F_y p {o(r['Fy_order'])} lim {l(r['Fy_limit'])} {d(r['Fy_limit_diff_pct'])} | finest {d(r['finest_uy_diff_pct'])} {d(r['finest_Fy_diff_pct'])}")
    print("\nlevels (reference / ours at the same dx - 1)")
    for r in lev:
        print(f"{r['source']:40s} dx {r['dx_m']:<7} u_y {r['uy_diff_pct']:+7.3f}%  F_y {r['Fy_diff_pct']:+7.3f}%")


if __name__ == "__main__":
    main()
