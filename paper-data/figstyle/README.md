# Paper figure style

Every plotting script uses `paper.mplstyle`:

- serif Computer Modern, 8 pt text and 7 pt ticks and legends, at final size;
- the Okabe–Ito colour-blind-safe cycle;
- viridis or cividis for sequential fields, and a diverging map (e.g. `RdBu_r`) centred on zero for errors;
- vector PDF output, with rasterised field layers (`rasterized=True`) at 400 dpi.

Widths in the sn-jnl layout:

- full text width: 372 pt = 5.16 in;
- half width: about 2.5 in.

Draw at the final width and include with `width=\textwidth` or the half width, so the fonts are not rescaled.

Other conventions:

- Label panels (a), (b), … at the top left, not with titles.
- Give units in the axis labels, e.g. `$u_y$ (mm)`.
- Keep file paths out of captions; put them in the paper-data READMEs.
