#!/usr/bin/env python3
"""Create one cavityFlexibleBottom preCICE case (SU2 incompressible + CalculiX quasi-static).

    setup_case.py LEVEL OUTDIR [--dt 0.5] [--t-end 100] [--inner-iter 1000] [--inner-minval -10]
                  [--steady-iter 20000] [--write-every 20] [--depth 1.0] [--min-its 0]

The runs in this folder used --dt 100 --t-end 4000 --write-every 10 (and --dt 5/25/50 for the
time-step checks, with --t-end 250/1000/2000).

* fluid mesh: make_fluid_mesh.py on the solids4foam tutorial fluid blockMeshDict (case/fluid-blockMeshDict,
  solids4foam development @ 1b5ab09), every block division x LEVEL, written to OUTDIR/cavity_L<LEVEL>.su2
  (OpenFOAM blockMesh in OUTDIR/foam-L<LEVEL>); the plate patch is already named 'interface';
* solid: solid_mesh.py, (40 L) x (4 L) C3D8I, x 2-6, y -0.4-0, slab depth --depth;
  fluid and solid are refined together (solid S = fluid LEVEL), so the interface nodes coincide;
* stage 1 steady rigid plate (steady_rigid.cfg), stage 2 coupled from it for --t-end s with
  coupling window = fluid time step --dt.
"""
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE.parent / "case"
# solids4foam development @ 1b5ab09: tutorials/fluidSolidInteraction/cavityFlexibleBottom/system/fluid/blockMeshDict
BLOCKMESH = TEMPLATE / "fluid-blockMeshDict"


def subst(path, **kw):
    t = path.read_text()
    for k, v in kw.items():
        t = t.replace(f"@{k}@", str(v))
    path.write_text(t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("level", type=int, choices=(1, 2, 4, 8, 16))
    ap.add_argument("outdir")
    ap.add_argument("--dt", type=float, default=0.5)
    ap.add_argument("--t-end", type=float, default=100.0)
    ap.add_argument("--inner-iter", type=int, default=1000)
    ap.add_argument("--inner-minval", type=float, default=-10.0)
    ap.add_argument("--steady-iter", type=int, default=20000)
    ap.add_argument("--write-every", type=int, default=20)
    ap.add_argument("--depth", type=float, default=1.0)
    ap.add_argument("--min-its", type=int, default=0,
                    help="preCICE <min-iterations> per window (0 = not set). L16 needs 2: with 1-iteration "
                         "windows IQN never acts and a slow mode creeps below the 1e-6 tolerance (job 869)")
    a = ap.parse_args()
    L = a.level
    out = Path(a.outdir).resolve()
    if out.exists():
        sys.exit(f"{out} exists")
    shutil.copytree(TEMPLATE, out, ignore=shutil.ignore_patterns("fluid-blockMeshDict"))

    fl = out / "fluid-su2"
    mesh = out / f"cavity_L{L}.su2"
    subprocess.run([sys.executable, "-I", HERE / "make_fluid_mesh.py", BLOCKMESH, str(L), mesh], check=True)
    m = mesh.read_text()
    assert m.count("MARKER_TAG= interface\n") == 1, "plate marker must be 'interface'"
    shutil.copy(mesh, fl / "fluid.su2")
    subprocess.run([sys.executable, "-I", HERE / "make_inlet_profile.py", fl / "fluid.su2",
                    fl / "inlet_profile.dat"], check=True, stdout=subprocess.DEVNULL)
    subst(fl / "steady_rigid.cfg", MESH="fluid.su2", STEADY_ITER=a.steady_iter)
    subst(fl / "coupled.cfg", MESH="fluid.su2", DT=repr(a.dt), INNER_ITER=a.inner_iter,
          INNER_MINVAL=a.inner_minval, WRT=a.write_every)

    sd = out / "solid-calculix"
    r = subprocess.run([sys.executable, "-I", HERE / "solid_mesh.py", str(L), sd / "mesh.inp",
                        "--depth", repr(a.depth)], check=True, capture_output=True, text=True)

    subst(out / "precice-config.base.xml", DT=repr(a.dt), TEND=repr(a.t_end), WRT=a.write_every)
    if a.min_its:
        x = out / "precice-config.base.xml"
        s = x.read_text()
        key = '    <max-iterations value="50" />\n'
        assert s.count(key) == 1
        x.write_text(s.replace(key, f'    <min-iterations value="{a.min_its}" />\n' + key))

    for f in out.rglob("*"):
        if f.is_file() and f.suffix in (".cfg", ".inp", ".xml", ".yml") and re.search(r"@[A-Z_]+@", f.read_text()):
            sys.exit(f"unfilled placeholder in {f}")
    info = (f"fluid level {L} ({mesh.name}), solid level {L}: {r.stdout.strip().split(': ', 1)[1]}; "
            f"dt {a.dt} s, coupled t_end {a.t_end} s, inner iter {a.inner_iter}, inner minval {a.inner_minval}, "
            f"steady iter {a.steady_iter}, write every {a.write_every} windows"
            + (f", preCICE min-iterations {a.min_its}" if a.min_its else ""))
    (out / "CASE_INFO.txt").write_text(info + "\n")
    print(info + f"\ncase in {out}")


if __name__ == "__main__":
    main()
