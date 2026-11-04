#!/usr/bin/env python3
"""Steady tip displacement of the coupled FSI1 runs from pointA.txt (columns x y ux uy at the tip node A = (0.6, 0.2)); displacement = x - 0.6, y - 0.2 is not used: ux, uy are logged."""
import sys, json, os
REF = dict(ux=2.270493e-5, uy=8.208773e-4)
def tip(run):
    cur = None; rows = []
    for l in open(run + "/pointA.txt"):
        if l.startswith("*Time"): cur = float(l.split("=")[1])
        elif not l.startswith("*") and l.strip(): v = l.split(); rows.append((cur, float(v[3]), float(v[4])))
    t = [r[0] for r in rows]; ux = [r[1] for r in rows]; uy = [r[2] for r in rows]
    last = [i for i, tt in enumerate(t) if tt >= t[-1] - 2.0 - 1e-9]
    spread = (max(uy[i] for i in last) - min(uy[i] for i in last)) / abs(uy[-1])
    return dict(t_end=t[-1], ux=ux[-1], uy=uy[-1], spread_last2s=spread, ux_err_pct=100 * (ux[-1] / REF["ux"] - 1), uy_err_pct=100 * (uy[-1] / REF["uy"] - 1))
if __name__ == "__main__":
    out = {}
    for r in sys.argv[1:]:
        try: out[os.path.basename(r.rstrip("/"))] = tip(r)
        except Exception as e: out[os.path.basename(r.rstrip("/"))] = str(e)
    for k, v in out.items(): print(k, v if isinstance(v, str) else {a: (round(b, 7) if isinstance(b, float) else b) for a, b in v.items()})
    json.dump(out, open("coupled_tip_results.json", "w"), indent=1)
