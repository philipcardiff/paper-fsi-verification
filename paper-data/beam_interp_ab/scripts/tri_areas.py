#!/usr/bin/env python3
"""Smallest AMI point-weight triangle on the 'interface' patch of each region.

The pre-#546 weights (triangle::pointToBarycentric) used the triangle
(face point i, face point i+1, face-point average) and fell back to 1/3 when
d00*d11 - d01^2 (= 4 A^2) < SMALL = 1e-15.
"""
import re, sys
from pathlib import Path
import numpy as np

def body(path):
    t = Path(path).read_text()
    t = re.sub(r"/\*.*?\*/", "", t, flags=re.S)
    t = re.sub(r"//[^\n]*", "", t)
    return t[t.index("}", t.index("FoamFile")) + 1:]

def points(mesh):
    t = body(mesh / "points")
    return np.array([[float(v) for v in m.split()] for m in re.findall(r"\(([^()]+)\)", t)])

def faces(mesh):
    t = body(mesh / "faces")
    t = t[t.index("(") + 1:]
    return [[int(v) for v in m.split()] for m in re.findall(r"\d+\(([^()]+)\)", t)]

def patch(mesh, name):
    t = body(mesh / "boundary")
    m = re.search(name + r"\s*\{([^}]*)\}", t)
    n = int(re.search(r"nFaces\s+(\d+)", m.group(1)).group(1))
    s = int(re.search(r"startFace\s+(\d+)", m.group(1)).group(1))
    return s, n

for case in sys.argv[1:]:
    for region in ("solid", "fluid"):
        mesh = Path(case) / "constant" / region / "polyMesh"
        p, f = points(mesh), faces(mesh)
        s, n = patch(mesh, "interface")
        areas, fareas = [], []
        for face in f[s:s + n]:
            x = p[face]
            c = x.mean(axis=0)
            for i in range(len(face)):
                a, b = x[i], x[(i + 1) % len(face)]
                v0, v1 = b - a, c - a
                d00, d01, d11 = v0 @ v0, v0 @ v1, v1 @ v1
                areas.append(np.sqrt(max(d00 * d11 - d01 * d01, 0.0)) / 2)
            fareas.append(0.5 * np.linalg.norm(np.cross(x[2] - x[0], x[3] - x[1])))
        areas = np.array(areas)
        print(f"{case} {region}: faces {n}, face area min {min(fareas):.4e} max {max(fareas):.4e}; "
              f"triangle area min {areas.min():.4e}; 4A^2 min {4*areas.min()**2:.4e}; "
              f"below threshold {(4*areas**2 < 1e-15).sum()}")
