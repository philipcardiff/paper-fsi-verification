"""collapsibleChannel fluid-mesh study figure for the verification paper
(fig:collapsible_mesh): (a) wall-midpoint vertical displacement for fluid
levels F1-F4 against the oomph-lib reference at the same time step; (b, c)
zooms at the first trough and the second peak; (d) maximum error against the
reference, as a percentage of the reference peak, against the number of
fluid cells, with the reference uncertainty marked.

Usage: collapsible_mesh_figure.py HIST_DIR REF_CSV REF_CONVERGED_CSV OUT.pdf [UNCERTAINTY [LABEL]]
  HIST_DIR: F1.dat ... F4.dat (solidPointDisplacement_wallMid.dat of each level)
    and mesh_study.csv of the tutorial driver (fluid cell counts)
  REF_CSV: oomph-lib solution at the same time step (collapsibleChannel_oomph_dt160.csv)
  REF_CONVERGED_CSV: converged reference (collapsibleChannel_oomph_reference.csv),
    whose peak normalises the error, as in the tutorial driver
  UNCERTAINTY: uncertainty of REF_CSV in % of the peak (default 0.3, the
    tutorial reference at --refine 1)
  LABEL: label of the shaded band in (d) (default "reference uncertainty")
Prints the per-level maximum error."""
import bisect, csv, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "figstyle"))
import figstyle  # noqa: E402  (paper style: CM fonts, 8 pt at final size)
import matplotlib.pyplot as plt  # noqa: E402

figstyle.use()
C = figstyle.OKABE_ITO
hist_dir, ref_csv, conv_csv, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
# Fluid cells per level from the driver's study CSV (the levels are not a
# uniform ratio-2 family: 3200, 12800, 28800, 51200)
CELLS = {int(row["fluidLevel"]): int(row["fluidCells"])
         for row in csv.DictReader(open(hist_dir / "mesh_study.csv"))}
REF_UNCERTAINTY = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3   # % of the peak
BAND_LABEL = sys.argv[6] if len(sys.argv) > 6 else "reference uncertainty"


def ref(path):
    rows = list(csv.DictReader(l for l in open(path) if not l.startswith("#")))
    return [(float(r["time"]), float(r["wallMid"])) for r in rows]


def hist(path):
    v = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        v[round(float(f[0]), 10)] = float(f[2])
    return sorted(v.items())


def interp(c, t):
    ts = [p[0] for p in c]
    i = min(max(bisect.bisect_left(ts, t), 1), len(c) - 1)
    (a, u), (b, w) = c[i - 1], c[i]
    return u + (w - u) * (t - a) / (b - a)


r = ref(ref_csv)
scale = max(abs(v) for _, v in ref(conv_csv))
runs = {l: hist(hist_dir / f"F{l}.dat") for l in CELLS}
err = {}
for l, c in runs.items():
    err[l] = 100 * max(abs(u - interp(r, t)) for t, u in c if t <= r[-1][0]) / scale
    print(f"F{l}: {CELLS[l]} cells, max error {err[l]:.3f}% of the reference peak")

fig = plt.figure(figsize=(figstyle.TEXTWIDTH, 3.9))
gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1], hspace=0.45, wspace=0.55,
                      left=0.105, right=0.985, bottom=0.105, top=0.95)
ax = fig.add_subplot(gs[0, :])
# reference underneath (thick, grey-black); levels on top, colour and line style
sty = {1: (C[1], ":"), 2: (C[4], "-."), 3: (C[2], "--"), 4: (C[0], "-")}
curves = [(r, "0.55", 2.2, "-", "oomph-lib (same $\\Delta t$)")] + [
    (runs[l], sty[l][0], 1.0, sty[l][1], f"F{l}") for l in sorted(runs)]
for c, col, lw, ls, lab in curves:
    ax.plot(*zip(*c), color=col, lw=lw, ls=ls, label=lab)
ax.set_xlabel("$t$ (s)", labelpad=1); ax.set_ylabel("wall-midpoint $u_y$ (m)")
ax.set_xlim(0, r[-1][0]); ax.set_ylim(-0.3, 0.015)
figstyle.panel_label(ax, "(a)")
ax.legend(loc="lower center", ncol=5, borderaxespad=0.25, handlelength=2.4, columnspacing=1.4)
trough_t = min(r, key=lambda p: p[1])[0]
peak_t = max((p for p in r if p[0] > trough_t + 0.2), key=lambda p: p[1])[0]
w = 0.05
for k, (t0, title) in enumerate(((trough_t, "(b)"), (peak_t, "(c)"))):
    zx = fig.add_subplot(gs[1, k])
    for c, col, lw, ls, _ in curves:
        seg = [p for p in c if t0 - w <= p[0] <= t0 + w]
        zx.plot(*zip(*seg), color=col, lw=lw, ls=ls)
    ys = [p[1] for c, *_ in curves for p in c if t0 - w <= p[0] <= t0 + w]
    pad = 0.06 * (max(ys) - min(ys))
    zx.set_xlim(t0 - w, t0 + w); zx.set_ylim(min(ys) - pad, max(ys) + pad)
    zx.set_xlabel("$t$ (s)", labelpad=1); figstyle.panel_label(zx, title)
    if k == 0:
        zx.set_ylabel("$u_y$ (m)")
    ax.axvspan(t0 - w, t0 + w, color="0.85", lw=0, zorder=0)
bx = fig.add_subplot(gs[1, 2])
n = [CELLS[l] for l in sorted(CELLS)]
e = [err[l] for l in sorted(CELLS)]
bx.axhspan(0.01, REF_UNCERTAINTY, color="0.85", lw=0)
bx.axhline(REF_UNCERTAINTY, color="0.4", lw=0.6, ls="--")
bx.text(n[-1], REF_UNCERTAINTY * 0.9, BAND_LABEL, color="0.25", fontsize=7,
        ha="right", va="top")
bx.loglog(n, e, "o-", color="black", label="max error")
bx.loglog([n[0], n[2]], [e[0], e[0] * n[0] / n[2]], ":", color="0.4", lw=0.8, label="order 2")
for l, x, y in zip(sorted(CELLS), n, e):
    first = l == min(CELLS)
    bx.annotate(f"F{l}", (x, y), textcoords="offset points", xytext=(4, 1) if first else (-3, -3),
                fontsize=7, ha="left" if first else "right", va="bottom" if first else "top")
bx.set_xticks(n); bx.set_xticklabels([f"{v / 1000:g}" for v in n]); bx.minorticks_off()
bx.set_xlabel("fluid cells ($10^3$)", labelpad=1); bx.set_ylabel("max error (% of peak)", labelpad=1)
bx.set_ylim(0.02, 8)
bx.set_yticks([0.05, 0.1, 0.2, 0.5, 1, 2, 5])
bx.set_yticklabels(["0.05", "0.1", "0.2", "0.5", "1", "2", "5"])
figstyle.panel_label(bx, "(d)")
bx.legend(loc="lower left", bbox_to_anchor=(0.0, 0.2), borderaxespad=0.2, handlelength=1.6)
fig.savefig(out)
