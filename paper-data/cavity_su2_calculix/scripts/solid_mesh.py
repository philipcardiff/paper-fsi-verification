#!/usr/bin/env python3
"""CalculiX mesh for the cavityFlexibleBottom plate: x 2-6 m, y -0.4-0 m, one C3D8I layer of
depth DEPTH (u_z = 0 on all nodes -> plane strain), (40 S) x (4 S) elements (solids4foam
tutorial solid blockMesh is 40 x 4 at S = 1, i.e. dx = 0.1 m).
Sets: Nall; Nfix = end faces x = 2 and x = 6 (tutorial 'left'/'right': fixedDisplacement 0);
Nsurface = top nodes y = 0 (the wetted interface; the bottom is traction-free);
Etop = top element row (face S5 is y = 0); NA = the two nodes at point A = (4, 0).

Usage: solid_mesh.py S OUT.inp [--depth 1.0] [--etype C3D8I]
"""
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("s", type=int)
ap.add_argument("out")
ap.add_argument("--depth", type=float, default=1.0)
ap.add_argument("--etype", default="C3D8I")
a = ap.parse_args()
x0, x1, y0, y1 = 2.0, 6.0, -0.4, 0.0
nx, ny = 40 * a.s, 4 * a.s
nid = lambda i, j, k: 1 + i + (nx + 1) * (j + (ny + 1) * k)
L = ["*NODE, NSET=Nall"]
fix, surf = [], []
for k in range(2):
    for j in range(ny + 1):
        for i in range(nx + 1):
            L.append(f"{nid(i,j,k)}, {x0 + (x1-x0)*i/nx:.12e}, {y0 + (y1-y0)*j/ny:.12e}, {a.depth*k:.12e}")
            if i in (0, nx):
                fix.append(nid(i, j, k))
            if j == ny:
                surf.append(nid(i, j, k))
L.append(f"*ELEMENT, TYPE={a.etype}, ELSET=Eall")
e, top = 0, []
for j in range(ny):
    for i in range(nx):
        e += 1
        L.append(f"{e}, " + ", ".join(map(str, [nid(i, j, 0), nid(i+1, j, 0), nid(i+1, j+1, 0), nid(i, j+1, 0),
                                                nid(i, j, 1), nid(i+1, j, 1), nid(i+1, j+1, 1), nid(i, j+1, 1)])))
        if j == ny - 1:
            top.append(e)


for kw, name, ids in (("NSET", "Nfix", fix), ("NSET", "Nsurface", surf), ("ELSET", "Etop", top),
                      ("NSET", "NA", [nid(nx // 2, ny, 0), nid(nx // 2, ny, 1)])):
    L.append(f"*{kw}, {kw}={name}")
    for b in range(0, len(ids), 10):
        L.append(", ".join(map(str, ids[b:b + 10])) + ",")
open(a.out, "w").write("\n".join(L) + "\n")
print(f"{a.out}: {nx}x{ny} {a.etype} ({e} elements), {len(surf)//2} interface points, depth {a.depth} m")
