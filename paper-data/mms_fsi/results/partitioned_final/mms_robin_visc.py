import csv, os, sys, math
R = "/Users/philipc/claude-runs/mms-robin-viscous/runs"
OUT = sys.argv[1]
runs = {16: f"{R}/cp-robin-visc-src/n16", 32: f"{R}/cp-robin-visc-src/n32",
        64: f"{R}/cp-robin-visc-src/n64", 128: f"{R}/cp-robin-visc-src128/n128"}
cols = ["U_L2","U_Linf","p_L2","p_Linf","D_L2","D_Linf","Us_L2","intD_L2","intTrac_L2","intTrac_Linf","intP_L2"]
data = {}
rows = []
iters = []
for n, run in runs.items():
    lines = open(os.path.join(run, "postProcessing", "mmsFsiErrors", "errors.dat")).read().splitlines()
    hdr = lines[0].lstrip("#").split()
    for l in lines[1:]:
        if not l.strip(): continue
        d = dict(zip(hdr, l.split())); t = float(d["time"])
        for T in (0.1, 0.2, 0.3, 0.4):
            if abs(t - T) < 1e-6:
                rows.append([f"n{n}", f"{T:g}"] + [d[c] for c in cols])
                data[(n, T)] = {c: float(d[c]) for c in cols}
    m = {}
    for l in open(os.path.join(run, "postProcessing", "fsiResiduals.dat")):
        if l.startswith("#") or not l.strip() or not l.split()[0][0].isdigit(): continue
        p = l.split(); t = float(p[0]); k = int(float(p[1])); m[t] = max(m.get(t, 0), k)
    vals = [m[t] for t in sorted(m) if t > 0]
    iters.append(["robin_full_visc", n, len(vals), f"{sum(vals)/len(vals):.1f}", max(vals)])
with open(os.path.join(OUT, "robin_full_visc_errors.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["mesh","time"] + cols); w.writerows(rows)
with open(os.path.join(OUT, "robin_full_visc_orders.md"), "w") as fh:
    for T in (0.2, 0.4):
        fh.write(f"### t = {T:g}\n\n| quantity | n16 | n32 | n64 | n128 | orders |\n|---|---|---|---|---|---|\n")
        for c in cols:
            e = [data[(n, T)][c] for n in (16, 32, 64, 128)]
            o = ", ".join(f"{math.log(e[i]/e[i+1], 2):.2f}" for i in range(3))
            fh.write(f"| {c} | " + " | ".join(f"{x:.3e}" for x in e) + f" | {o} |\n")
        fh.write("\n")
for r in iters: print(r)
