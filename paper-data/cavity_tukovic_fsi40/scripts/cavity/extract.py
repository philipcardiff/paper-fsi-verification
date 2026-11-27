"""Point A (top mid-span vertex of the plate) u_y at every written time, and
the last-iteration fluid force per time step.  Usage: python3 -I extract.py CASE R"""
import os, re, sys
case, R = sys.argv[1], int(sys.argv[2])
nx, ny = 40 * R, 4 * R
iA = nx // 2 + ny * (nx + 1)
fl = os.path.join(case, "fluid")
times = sorted((d for d in os.listdir(fl) if re.fullmatch(r"[0-9.]+", d) and d != "0"), key=float)
for t in times:
    p = os.path.join(fl, t, "solid", "pointD")
    if not os.path.exists(p):
        continue
    txt = open(p).read()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*\n?(\d+)\s*\n?\(", txt)
    vals = re.findall(r"\(([^()]*)\)", txt[m.end():])
    print("t", t, "uA", vals[iA])
log = open(os.path.join(fl, "log.fsiFoam")).read()
blocks = re.split(r"\nTime = ([0-9.e+-]+) \(dt = [^)]*\)\n", log)
F = {}
for t, body in zip(blocks[1::2], blocks[2::2]):
    f = re.findall(r"Total force \(fluid\) = \(([^)]*)\)", body)
    it = re.findall(r"Time = [0-9.e+-]+, iteration: (\d+)", body)
    if f:
        F[t] = (f[-1], it[-1] if it else "?")
ks = list(F)
for t in ks[-3:]:
    print("t", t, "Ffluid", F[t][0], "iters", F[t][1])
print("max iters", max(int(v[1]) for v in F.values() if v[1] != "?"))
