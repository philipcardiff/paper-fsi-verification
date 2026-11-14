#!/usr/bin/env python3
"""Create one Hron-Turek FSI1 preCICE case (OpenFOAM + CalculiX) for a mesh level.

    setup-case.py LEVEL OUTDIR [SOLID_LEVEL]
        LEVEL = fluid refinement factor 1, 2 or 4 (also sets the time step)
        SOLID_LEVEL = solid refinement factor (default LEVEL), for solid-only refinement

Fluid: the solids4foam TurekHron2006 blockMesh with every in-plane block division
multiplied by LEVEL. Solid: (105*LEVEL) x (6*LEVEL) x 1 C3D8I elements, 0.015 m deep
(the fluid mesh depth). Time step as in the solids4foam FSI1 verification.
"""
import re, shutil, sys
from pathlib import Path

level = int(sys.argv[1]); out = Path(sys.argv[2])
slevel = int(sys.argv[3]) if len(sys.argv) > 3 else level
dt = {1: 0.025, 2: 0.025, 4: 0.0125}[level]
tpl = Path(__file__).resolve().parent.parent / "case"
if out.exists():
    sys.exit(f"{out} exists")
shutil.copytree(tpl, out)

# coupling / solid time step
for f in (out / "precice-config.xml", out / "solid-calculix/plate.inp"):
    f.write_text(f.read_text().replace("@DT@", repr(dt)))
cd = out / "fluid-openfoam/system/controlDict"
cd.write_text(re.sub(r"deltaT\s+[^;]+;.*", f"deltaT          {dt};", cd.read_text()))

# fluid mesh refinement: hex (...) (nx ny 1)
bm = out / "fluid-openfoam/system/blockMeshDict"
bm.write_text(re.sub(r"(hex\s*\([^)]*\)\s*)\((\d+)\s+(\d+)\s+1\)",
                     lambda m: f"{m.group(1)}({int(m.group(2))*level} {int(m.group(3))*level} 1)",
                     bm.read_text()))

# solid mesh
x0, x1, y0, y1, depth = 0.2489898, 0.6, 0.19, 0.21, 0.015
nx, ny = 105 * slevel, 6 * slevel
nid = lambda i, j, k: 1 + i + (nx + 1) * (j + (ny + 1) * k)
lines = ["*NODE, NSET=Nall"]
fix, surf = [], []
for k in range(2):
    for j in range(ny + 1):
        for i in range(nx + 1):
            x = x0 + (x1 - x0) * i / nx
            y = y0 + (y1 - y0) * j / ny
            lines.append(f"{nid(i,j,k)}, {x:.12e}, {y:.12e}, {depth*k:.12e}")
            if i == 0:
                fix.append(nid(i, j, k))
            if j == 0 or j == ny or i == nx:
                surf.append(nid(i, j, k))
lines.append("*ELEMENT, TYPE=C3D8I, ELSET=Eall")
e = 0
for j in range(ny):
    for i in range(nx):
        e += 1
        n = [nid(i,j,0), nid(i+1,j,0), nid(i+1,j+1,0), nid(i,j+1,0),
             nid(i,j,1), nid(i+1,j,1), nid(i+1,j+1,1), nid(i,j+1,1)]
        lines.append(f"{e}, " + ", ".join(map(str, n)))
def nset(name, ids):
    lines.append(f"*NSET, NSET={name}")
    for a in range(0, len(ids), 10):
        lines.append(", ".join(map(str, ids[a:a+10])) + ",")
nset("Nfix", fix)
nset("Nsurface", surf)   # top, bottom and tip faces (the fluid 'plate' patch)
(out / "solid-calculix/mesh.inp").write_text("\n".join(lines) + "\n")
print(f"fluid level {level}, solid level {slevel}: dt {dt}, solid {nx}x{ny} C3D8I ({e} elements, "
      f"{len(surf)//2} interface points), case in {out}")
