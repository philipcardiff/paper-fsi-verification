# 3-D pressure-pulse tube, perspective figure (`figures/3dtube-3d.pdf`, `fig:3dtube_3d`)

Source run: the solution-gallery rerun of the paper's level-1 case (corrected
meshes, `robin_backward_mesh1` of the 3dTube mesh study: Robin-Neumann,
backward, Delta t = 2.5e-5 s, OpenFOAM-v2512, 8 ranks) to t = 5 ms,
reconstructed (see `../solution_gallery/README.md`, panel f).

Pipeline (ParaView 6.1.1 `pvbatch` for 1-2, matplotlib for 3):

1. `extract.py CASE.foam INIT.foam OUTDIR [--time 0.005] [--rho 1000]`
   writes, on the undeformed geometry (mm; committed in `data/`):
   fluid `symmetry-y`, `symmetry-x`, `inlet`, `outlet` patches with `p` (Pa,
   rho times the kinematic pressure) and `disp` (fluid mesh motion, mm), and
   the solid `outer-wall`, symmetry, inlet and outlet patches with `disp`
   (`D`, mm) and `ur` (radial displacement, mm). `INIT.foam` is a stub in a
   directory with links to the case's `constant/` and `system/` only, so that
   the reader returns the initial fluid mesh (it does not list the case's `0`
   directory).
2. `render.py DATADIR render.png --width 3600 --height 1400` warps wall and
   fluid mesh by ten times their displacement (`--scale`), colours the fluid
   faces by p (viridis, 0-1500 Pa, unlit) and the wall by u_r (Oranges,
   0-0.16 mm), and writes `render.json` with label anchor points. The quarter
   (x, y >= 0) is drawn mirrored about y = 0 and seen from +y, so the flow
   runs left to right; only the plane y = 0 is visible from any viewpoint
   that also shows the outer wall.
3. `fig.py [--render render.png] --out <repo>/figures/3dtube-3d.pdf` adds the
   labels, axis triad and both colour bars in the paper style (`../figstyle`).

Rebuild from the committed data:

    PV=/Applications/ParaView-6.1.1.app/Contents/bin/pvbatch
    $PV render.py data render.png --width 3600 --height 1400
    python3 fig.py --out ../../figures/3dtube-3d.pdf
