"""Extract per-time-step data from our solids4foam feF1 run (log.solids4Foam + flapTip function object). usage: extract_s4f.py CASEDIR OUT.csv
Columns: t, tip_dy_m, tip_dz_m, Ffluid(x,y,z) 'Total force on fluid interface', Fsolid(x,y,z) 'Total force on solid interface' (last coupling iteration of the step), n_coupling_its, r2 (last FSI residual2)."""
import re, sys, csv, glob
case, out = sys.argv[1], sys.argv[2]
rows = {}; cur = None; its = 0
pat = lambda s: re.compile(r'^Total force on %s interface 0: \((\S+) (\S+) (\S+)\)' % s)
pf, ps = pat("fluid"), pat("solid")
for fn in sorted(glob.glob(case + "/log.solids4Foam*")):
    if fn.endswith("pilot_0_0.1"): continue         # the pilot log is the first 0.1 s of the same run (restart appends to log.solids4Foam)
    for l in open(fn):
        m = re.match(r'^Time = (\S+)$', l)
        if m: cur = round(float(m.group(1)), 6); its = 0; rows[cur] = rows.get(cur, {}); continue
        if cur is None: continue
        if l.startswith("Time =") and "iteration:" in l: its += 1; rows[cur]["its"] = its
        m = pf.match(l)
        if m: rows[cur]["Ff"] = tuple(map(float, m.groups()))
        m = ps.match(l)
        if m: rows[cur]["Fs"] = tuple(map(float, m.groups()))
        m = re.search(r'FSI residual2 norm for interface 0: (\S+)', l)
        if m: rows[cur]["r2"] = float(m.group(1))
tip = {}
for l in open(case + "/postProcessing/0/solidPointDisplacement_flapTip.dat"):
    if l.startswith("#"): continue
    a = l.split(); tip[round(float(a[0]), 6)] = (float(a[2]), float(a[3]))
with open(out, "w") as f:
    w = csv.writer(f); w.writerow(["t", "tip_dy_m", "tip_dz_m", "Ff_x", "Ff_y", "Ff_z", "Fs_x", "Fs_y", "Fs_z", "its", "r2"])
    for t in sorted(rows):
        r = rows[t]
        if t not in tip or "Ff" not in r: continue
        w.writerow([t, tip[t][0], tip[t][1], *r["Ff"], *r.get("Fs", (float("nan"),) * 3), r.get("its", ""), r.get("r2", "")])
print("rows", len(rows))
