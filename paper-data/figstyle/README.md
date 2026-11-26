# Paper figure style

Every plotting script uses `paper.mplstyle`:

- serif Computer Modern, 8 pt text and 7 pt ticks and legends, at final size;
- the Okabe–Ito colour-blind-safe cycle;
- viridis or cividis for sequential fields, and a diverging map (e.g. `RdBu_r`) centred on zero for errors;
- vector PDF output, with rasterised field layers (`rasterized=True`) at 400 dpi.

Widths in the layout of `article.tex`:

- full text width: 15 cm = 427 pt = 5.93 in (`textwidth=15cm` in `article.tex`);
- half width: about 209 pt = 2.9 in.

`figstyle.py` holds these widths (`TEXTWIDTH`, `HALFWIDTH`) and `figstyle.use()` loads this style with the CMU Serif fonts of the local TeX installation.

Draw at the final width and include with `width=\textwidth` or the half width, so the fonts are not rescaled.

Other conventions:

- Label panels (a), (b), … at the top left, not with titles.
- Give units in the axis labels, e.g. `$u_y$ (mm)`.
- Keep file paths out of captions; put them in the paper-data READMEs.
