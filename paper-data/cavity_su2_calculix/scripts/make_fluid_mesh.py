#!/usr/bin/env python3
"""2-D blockMesh -> SU2 mesh (as ../fsi1_su2_calculix/scripts/make_fluid_mesh.py; used here for
cavityFlexibleBottom with the solids4foam tutorial's fluid blockMeshDict, case/fluid-blockMeshDict).

Takes a 2-D (one cell deep) blockMeshDict, multiplies every in-plane block division by LEVEL
(the same substitution as PR #75's setup-case.py), runs OpenFOAM
v2412 blockMesh, and converts the one-cell-deep polyMesh to a 2-D .su2 mesh:
  * elements = the z = zmin faces of the frontAndBackPlanes patch (quads, made CCW);
  * nodes    = the z = zmin points;
  * markers  = each non-empty patch, as the z = zmin edges of its faces
               (cavity: inlet, outlet, interface, walls).

Usage: make_fluid_mesh.py BLOCKMESHDICT LEVEL OUT.su2 [--workdir DIR]
"""
import argparse
import re
import subprocess
from pathlib import Path

import numpy as np

OF_BASHRC = "/usr/lib/openfoam/openfoam2412/etc/bashrc"
CONTROL_DICT = """FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }
application blockMesh; startFrom startTime; startTime 0; stopAt endTime; endTime 1;
deltaT 1; writeControl timeStep; writeInterval 1; writeFormat ascii; writePrecision 16;
writeCompression off; timeFormat general; timePrecision 6; runTimeModifiable false;
"""


def refine_blockmesh(text, level):
    # identical substitution to PR #75 scripts/setup-case.py
    return re.sub(r"(hex\s*\([^)]*\)\s*)\((\d+)\s+(\d+)\s+1\)",
                  lambda m: f"{m.group(1)}({int(m.group(2))*level} {int(m.group(3))*level} 1)", text)


def foam_list_body(path):
    """Return the text between the outermost '(' ... ')' of an OpenFOAM list file."""
    txt = Path(path).read_text()
    txt = re.sub(r"/\*.*?\*/", "", txt, flags=re.S)
    txt = re.sub(r"//[^\n]*", "", txt)
    txt = txt[txt.index("}") + 1:]                      # skip FoamFile header
    m = re.search(r"(\d+)\s*\(", txt)
    n = int(m.group(1))
    start = m.end()
    end = txt.rindex(")")
    return n, txt[start:end]


def read_points(p):
    n, body = foam_list_body(p)
    pts = np.array(re.findall(r"\(([^()]+)\)", body)[:n], dtype=object)
    return np.array([list(map(float, s.split())) for s in pts])


def read_faces(p):
    n, body = foam_list_body(p)
    faces = [list(map(int, s.split())) for s in re.findall(r"\d+\(([^()]+)\)", body)]
    assert len(faces) == n, (len(faces), n)
    return faces


def read_labels(p):
    n, body = foam_list_body(p)
    lab = np.array(body.split(), dtype=int)
    assert len(lab) == n
    return lab


def read_boundary(p):
    txt = re.sub(r"/\*.*?\*/", "", Path(p).read_text(), flags=re.S)
    txt = txt[txt.index("}") + 1:]
    out = []
    for name, blk in re.findall(r"(\w+)\s*\{([^{}]*)\}", txt):
        typ = re.search(r"type\s+(\w+);", blk).group(1)
        nf = int(re.search(r"nFaces\s+(\d+);", blk).group(1))
        sf = int(re.search(r"startFace\s+(\d+);", blk).group(1))
        out.append((name, typ, nf, sf))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("blockmeshdict")
    ap.add_argument("level", type=int)
    ap.add_argument("out")
    ap.add_argument("--workdir", default=None, help="OpenFOAM case dir (default: OUT dir/foam-L<level>)")
    a = ap.parse_args()
    out = Path(a.out).resolve()
    wd = Path(a.workdir) if a.workdir else out.parent / f"foam-L{a.level}"
    (wd / "system").mkdir(parents=True, exist_ok=True)
    (wd / "system/controlDict").write_text(CONTROL_DICT)
    (wd / "system/blockMeshDict").write_text(refine_blockmesh(Path(a.blockmeshdict).read_text(), a.level))
    subprocess.run(["bash", "-c", f"set +u; source {OF_BASHRC} >/dev/null 2>&1; "
                    f"blockMesh -case '{wd}' > '{wd}/log.blockMesh' 2>&1"], check=True)

    pm = wd / "constant/polyMesh"
    P = read_points(pm / "points")
    F = read_faces(pm / "faces")
    owner = read_labels(pm / "owner")
    bnd = read_boundary(pm / "boundary")
    z0 = P[:, 2].min()
    front = np.abs(P[:, 2] - z0) < 1e-12 * max(1.0, abs(P[:, 2]).max())
    new_id = -np.ones(len(P), dtype=int)
    new_id[front] = np.arange(front.sum())

    elems, markers = [], {}
    for name, typ, nf, sf in bnd:
        if typ == "empty":
            for f in F[sf:sf + nf]:
                if all(front[f]):
                    q = [int(new_id[v]) for v in f]
                    x, y = P[f, 0], P[f, 1]
                    area2 = np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y)
                    if area2 < 0:
                        q = q[::-1]
                    elems.append(q)
        else:
            edges = []
            for f in F[sf:sf + nf]:
                fr = [v for v in f if front[v]]
                assert len(fr) == 2, (name, f)
                edges.append(tuple(int(new_id[v]) for v in fr))
            markers[name] = edges
    n_cells = len(owner) and (owner.max() + 1)
    assert len(elems) == n_cells, f"{len(elems)} front faces vs {n_cells} cells"

    xy = P[front][:, :2]
    with open(out, "w") as fh:
        fh.write("NDIME= 2\n")
        fh.write(f"NELEM= {len(elems)}\n")
        for e, q in enumerate(elems):
            fh.write(f"9 {q[0]} {q[1]} {q[2]} {q[3]} {e}\n")
        fh.write(f"NPOIN= {len(xy)}\n")
        for i, (x, y) in enumerate(xy):
            fh.write(f"{x:.15e} {y:.15e} {i}\n")
        fh.write(f"NMARK= {len(markers)}\n")
        for name, edges in markers.items():
            fh.write(f"MARKER_TAG= {name}\nMARKER_ELEMS= {len(edges)}\n")
            for n0, n1 in edges:
                fh.write(f"3 {n0} {n1}\n")
    print(f"{out}: level {a.level}, {len(elems)} quads, {len(xy)} nodes, markers "
          + ", ".join(f"{k} {len(v)}" for k, v in markers.items()))


if __name__ == "__main__":
    main()
