# Figure sources

Where each figure of the paper comes from. Captions carry no file paths; this
file and the README of each data directory do. All matplotlib figures use the
shared style (`../figstyle/paper.mplstyle`, loaded by `../figstyle/figstyle.py`,
which also registers the CMU Serif fonts of the local TeX installation), are
drawn at their printed width (the 427 pt text width, `figstyle.TEXTWIDTH`) and are included in the `.tex` at their
natural size, so the 8 pt text is not rescaled.

Python 3 with numpy and matplotlib (fontTools comes with matplotlib). Paths
below are relative to the repository root; `$TUT` is a solids4foam checkout's
`tutorials/fluidSolidInteraction`.

| Label | File (`figures/`) | Script | Inputs | Command |
|---|---|---|---|---|
| `fig:geometry_panel` | `benchmarks_geometry_panel.pdf` | `paper-data/figure_scripts/geometry_panel.py` | none (dimensions of the tutorials) | `python3 paper-data/figure_scripts/geometry_panel.py figures/benchmarks_geometry_panel.pdf` |
| `fig:ring_febio_ratio`, `fig:ring_febio_decay`, `fig:ring_febio_traction` | `verification/ringAddedMass_febio_R{1,3,4}.pdf` | `paper-data/febio_ring_added_mass/no_slip_resolved/scripts/ring_paper_figures.py` | stored results of that package and its parent's `reference/` and `tools/` | `python3 paper-data/febio_ring_added_mass/no_slip_resolved/scripts/ring_paper_figures.py figures/verification` |
| `fig:mms_fsi_case` | `mms_fsi_case.pdf` | `paper-data/mms_fsi/figure/mms_fsi_case_figure.py` | the exact evaluators built from `paper-data/mms_fsi/case/src/evalExact.C` (see `paper-data/mms_fsi/README.md`) | `python3 paper-data/mms_fsi/figure/mms_fsi_case_figure.py EVAL_FULL EVAL_FLAT figures/mms_fsi_case.pdf` |
| `fig:mms_fsi` | `mms_fsi_pressure_error.pdf` | `paper-data/mms_fsi/figure/mms_fsi_figure.py` | the n128 Aitken and QM run directories at t = 0.4 (not in the repository; see `paper-data/mms_fsi/README.md`) | `python3 paper-data/mms_fsi/figure/mms_fsi_figure.py AITKEN_RUN QM_RUN figures/mms_fsi_pressure_error.pdf` |
| `fig:collapsible_mesh` | `verification/collapsibleChannel_mesh_study.pdf` | `paper-data/figure_scripts/collapsible_mesh_figure.py` | `paper-data/collapsible_mesh_study/`, `$TUT/collapsibleChannel/verification/reference/collapsibleChannel_oomph_reference.csv` | `python3 paper-data/figure_scripts/collapsible_mesh_figure.py paper-data/collapsible_mesh_study paper-data/collapsible_mesh_study/oomph_extrapolated_dt160.csv $TUT/collapsibleChannel/verification/reference/collapsibleChannel_oomph_reference.csv figures/verification/collapsibleChannel_mesh_study.pdf 0.05 "model form"` |
| `fig:mok_refs` | `verification/mokLidDrivenCavity_study.pdf` | `paper-data/figure_scripts/mok_figure.py` | `$TUT/mokLidDrivenCavity/verification/reference`, `paper-data/mok_coupling_test/`, `paper-data/kratos_mok_lid_driven_cavity/` | `python3 paper-data/figure_scripts/mok_figure.py $TUT/mokLidDrivenCavity/verification/reference paper-data/mok_coupling_test figures/verification/mokLidDrivenCavity_study.pdf` (prints `paper-data/mok_coupling_test/statistics.csv`) |
| `fig:hessenthaler_geom` | `hessenthaler-geometry.pdf` | `paper-data/hessenthaler_geometry_figure/fig.py` | the data files in that directory | `python3 paper-data/hessenthaler_geometry_figure/fig.py --out figures/hessenthaler-geometry.pdf` |
| `fig:hronturek_geom` | `HronTurek3-geometry.pdf` | Inkscape, `figures-svg/HronTurekFsi3.svg` | — | Inkscape export of the page; the PDF was then cropped to the drawing (406 × 110 pt) |
| `fig:3dtube_geom`, `fig:cavity_geom_v`, `fig:blob_geom` | `3dTube-geometry.pdf`, `cavityFlexibleBottom-geometry.pdf`, `blobInTreacle-geometry.pdf` | Inkscape, `figures-svg/*.svg` | — | Inkscape export |

The Inkscape figures have their text converted to paths (Computer Modern from
TexText) and are included at `scale=1`.
