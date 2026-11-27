# Hessenthaler flap, 3-D perspective figure (`figures/hessenthaler-3d.pdf`, `fig:hessenthaler_3d`)

Complements the mid-plane views of `../hessenthaler_geometry_figure/`.

Source run: the partitioned L2 run (`rF1` in the run names), UCD Sonic
`/scratch/pcardiff/hv/P_rF1_LUST_b02_mwp`, reconstructed time `45` (s), copied
read-only to the Mac (`constant/`, `system/`, `45/`).

Pipeline (ParaView 6.1.1 `pvbatch` for 1-2, matplotlib for 3):

1. `extract.py CASE.foam OUTDIR [--time 45]` reads the case with ParaView's
   OpenFOAM reader (`CASE.foam` is an empty stub file in the case directory)
   and writes small surface files in mm to `OUTDIR` (committed in `data/`):
   `wall.vtp`, `inlets.vtp`, `outlet.vtp` (fluid patches at t = 45 s),
   `flap0.vtp` (undeformed solid surface), `flap.vtp` (warped by `pointD`,
   with `magD` in mm) and `slice_x0.vtp` (fluid cut at x = 0 on the moved
   mesh, with `magU` in m/s).
2. `render.py DATADIR render.png --width 3200 --height 1520 --zmax 185`
   renders the view (camera in `camera.json`; `--zmax` cuts the domain at
   z = 185 mm) and writes `render.json` with the image coordinates of the
   anchor points used for the labels. No text is rendered here.
3. `fig.py [--render render.png] --out <repo>/figures/hessenthaler-3d.pdf`
   adds labels, axis triad, gravity and the |u| colour bar in the paper style
   (`../figstyle`), at the full text width, saved at 400 dpi (the render has
   about 650 px per inch at final size).

Rebuild from the committed data:

    PV=/Applications/ParaView-6.1.1.app/Contents/bin/pvbatch
    $PV render.py data render.png --width 3200 --height 1520 --zmax 185
    python3 fig.py --out ../../figures/hessenthaler-3d.pdf

Colours: |u| viridis (unlit, so the colours match the bar), flap orange as in
the 2-D figure, inlets/outlet Okabe-Ito blue/vermillion. The tip of the
deformed flap is at y = 16.4 mm (consistent with the L2 value in the text).
