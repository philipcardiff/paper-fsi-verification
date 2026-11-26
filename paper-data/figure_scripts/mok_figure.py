"""mokLidDrivenCavity figure for the verification paper (fig:mok_refs):
(a) membrane-midpoint histories over the last four lid periods: the
published group-A and group-B curves and the solids4foam group-B solutions
with the strict coupling test; (b) peak, mean and trough over t = 20-70 s
against the fluid cells per side, legacy and strict coupling test, with the
published group-B envelope, and the independent Kratos solution of the
group-B definition (openings fixed; ../kratos_mok_lid_driven_cavity/).

The Kratos tutorial example (reference "Kratos") lets its openings move with
the mesh (#36), so it is not a group-B solution: it is plotted but left out of
the group-B envelope, whatever the reference file says.

Usage: mok_figure.py REF_DIR RUN_DIR OUT.pdf [KRATOS_DIR]
  REF_DIR: the tutorial's verification/reference directory (published curves)
  RUN_DIR: ../mok_coupling_test (<variant>-mesh{N}.csv, columns t, u_y (m)).
    The series plotted are fixed in SERIES below, as in that directory's
    README: the legacy series is the paper's mesh study (legacy0 for 16^2,
    whose reuse-5 run stopped at t = 60.2 s), the strict points the completed
    levels with the strict coupling test (requireAllResidualMeasures)
Prints the statistics of every run as CSV on stdout."""
import csv, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "figstyle"))
import figstyle  # noqa: E402  (paper style: CM fonts, 8 pt at final size)
import matplotlib.pyplot as plt  # noqa: E402

figstyle.use()
C = figstyle.OKABE_ITO
ref_dir, run_dir, out = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
kratos_dir = Path(sys.argv[4]) if len(sys.argv) > 4 else \
    Path(__file__).resolve().parent.parent / "kratos_mok_lid_driven_cavity"
MOVING_OPENINGS = {"Kratos"}  # the Kratos tutorial example (#36)
cfg = json.loads((ref_dir / "mokLidDrivenCavity_verification_references.json").read_text())
T0, T1 = cfg["comparison_window_s"]
PERIOD = cfg["lid_period_s"]


def read(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        try:
            rows.append(tuple(float(v) for v in line.replace(",", " ").split()[:2]))
        except ValueError:
            continue
    return sorted(rows)


def interp(c, t):
    if t <= c[0][0]:
        return c[0][1]
    if t >= c[-1][0]:
        return c[-1][1]
    lo, hi = 0, len(c) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if c[mid][0] <= t else (lo, mid)
    (a, u), (b, v) = c[lo], c[hi]
    return u + (v - u) * (t - a) / (b - a)


def stats(c, start=T0, end=T1):
    """As periodic_stats in the tutorial driver: mean per-period peak and
    trough, and the time-mean over the window."""
    end = min(end, c[-1][0])
    peaks, troughs, s = [], [], start
    while s + PERIOD <= end + 1e-9:
        w = [u for t, u in c if s <= t < s + PERIOD]
        if len(w) >= 10:
            peaks.append(max(w)); troughs.append(min(w))
        s += PERIOD
    mean = sum(interp(c, start + (end - start) * i / 2000) for i in range(2001)) / 2001
    return dict(peak=sum(peaks) / len(peaks), trough=sum(troughs) / len(troughs), mean=mean)


refs = {k: dict(v) for k, v in cfg["references"].items()}
for k in MOVING_OPENINGS:
    refs[k]["group_b"] = False
curves = {k: read(ref_dir / v["file"]) for k, v in refs.items() if "file" in v}
group_b = {k: (stats(curves[k]) if k in curves else refs[k]["values"])
           for k, v in refs.items() if v["group_b"]}
env = {q: (min(s[q] for s in group_b.values() if q in s),
           max(s[q] for s in group_b.values() if q in s)) for q in ("peak", "trough", "mean")}

# (variant, n) -> file in RUN_DIR (see ../mok_coupling_test/README.md)
SERIES = {("legacy", 16): "legacy0-mesh16.csv", ("legacy", 32): "legacy-mesh32.csv",
          ("legacy", 64): "legacy-mesh64.csv", ("legacy", 96): "legacy-mesh96.csv",
          ("strict", 32): "strict-mesh32.csv", ("strict", 64): "strict-mesh64.csv"}
runs = {k: read(run_dir / f) for k, f in SERIES.items()}
w = csv.writer(sys.stdout)
w.writerow(["variant", "nx", "peak", "trough", "mean"])
table = {}
# Independent Kratos solution, openings fixed: the same combined sequence
kratos = {}
for nx, dt in ((16, "0.2"), (32, "0.1"), (64, "0.05"), (96, "0.025")):
    h = kratos_dir / "runs" / f"n{nx:03d}_dt{dt}_sy2_fixopen" / "history.csv"
    rows = list(csv.DictReader(h.open()))
    kratos[nx] = sorted((float(r["time"]), float(r["uy_mid"])) for r in rows)
for nx, c in sorted(kratos.items()):
    s = stats(c)
    w.writerow(["kratos_fixed_openings", nx] + [f"{s[q]:.5f}" for q in ("peak", "trough", "mean")])
for (variant, nx), c in sorted(runs.items()):
    s = stats(c); table[(variant, nx)] = s
    w.writerow([variant, nx] + [f"{s[q]:.5f}" for q in ("peak", "trough", "mean")])
for q, (lo, hi) in env.items():
    w.writerow(["groupB_envelope", q, f"{lo:.4f}", f"{hi:.4f}", ""])

fig, (ax, bx) = plt.subplots(1, 2, figsize=(figstyle.TEXTWIDTH, 3.45),
                             gridspec_kw=dict(width_ratios=[1.55, 1]))
t_show = (T1 - 4 * PERIOD, T1)
GA = dict(Mok2001="Mok (2001)", Wall1999="Wall (1999)", GerbeauVidrascu2003="Gerbeau & Vidrascu (2003)")
for (k, lab), ls in zip(GA.items(), ("-", "--", ":")):
    c = [p for p in curves[k] if t_show[0] <= p[0] <= t_show[1]]
    if c:
        ax.plot(*zip(*c), ls, color="0.6", lw=0.8, label=f"{lab} (A)")
c = [p for p in curves["Kassiotis2011"] if t_show[0] <= p[0] <= t_show[1]]
ax.plot(*zip(*c), "-.", color="0.45", lw=0.8, label="Kassiotis et al. (2011)")
c = [p for p in curves["Kratos"] if t_show[0] <= p[0] <= t_show[1]]
ax.plot(*zip(*c), "--", color=C[4], lw=0.9, label="Kratos example (moving openings)")
c = [p for p in curves["Valdes2007"] if t_show[0] <= p[0] <= t_show[1]]
ax.plot(*zip(*c), color=C[0], lw=1.0, label="Valdés (2007) (B)")
# solids4foam underneath (thick black), the Kratos 96^2 solution dashed on top:
# the two curves nearly coincide
finest = max(nx for v, nx in runs if v == "legacy")
c = [p for p in runs[("legacy", finest)] if t_show[0] <= p[0] <= t_show[1]]
ax.plot(*zip(*c), color="black", lw=1.6, label=f"solids4foam ${finest}^2$ (B)", zorder=3)
c = [p for p in kratos[96] if t_show[0] <= p[0] <= t_show[1]]
ax.plot(*zip(*c), color=C[2], lw=1.0, ls=(0, (3, 1.5)), label="Kratos $96^2$, openings fixed (B)",
        zorder=4)
ax.set_xlim(*t_show); ax.set_xlabel("$t$ (s)", labelpad=1); ax.set_ylabel("midpoint $u_y$ (m)")
figstyle.panel_label(ax, "(a)")
ax.legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.2), borderaxespad=0.0,
          columnspacing=1.0, handlelength=2.2)

for q, mk in (("peak", "o"), ("mean", "s"), ("trough", "v")):
    lo, hi = env[q]
    if hi - lo > 1e-4:
        bx.axhspan(lo, hi, color=C[0], alpha=0.15, lw=0)
    else:  # a single group-B value (the trough: Valdes only)
        bx.axhline(lo, color=C[0], alpha=0.5, lw=1.5)
    pts = sorted((nx, s[q]) for (v, nx), s in table.items() if v == "legacy")
    bx.plot([1 / nx for nx, _ in pts], [y for _, y in pts], marker=mk, color="black", mfc="black",
            label="solids4foam, legacy test" if q == "peak" else None, zorder=3)
    # strict test: larger open markers drawn over the legacy ones
    pts = sorted((nx, s[q]) for (v, nx), s in table.items() if v == "strict")
    bx.plot([1 / nx for nx, _ in pts], [y for _, y in pts], marker=mk, ls="none", ms=6.5, mew=1.0,
            color=C[1], mfc="none", label="solids4foam, strict test" if q == "peak" else None, zorder=4)
    pts = sorted((nx, stats(c)[q]) for nx, c in kratos.items())
    bx.plot([1 / nx for nx, _ in pts], [y for _, y in pts], marker=mk, ls="--", color=C[2],
            mfc="white", label="Kratos, openings fixed" if q == "peak" else None, zorder=2)
bx.set_xscale("log"); bx.set_xlabel("$1/n$ ($n$: fluid cells per side)", labelpad=1)
nxs = sorted({nx for _, nx in table})
bx.set_xticks([1 / n for n in nxs]); bx.set_xticklabels([f"1/{n}" for n in nxs]); bx.minorticks_off()
bx.set_xlim(1 / 115, 1 / 13.5)
bx.set_ylabel("$u_y$ (m): peak, mean, trough")
figstyle.panel_label(bx, "(b)")
for q in ("peak", "mean", "trough"):
    y = table[("legacy", max(nx for v, nx in table if v == "legacy"))][q]
    bx.annotate(q, (1 / max(nxs), y), textcoords="offset points", xytext=(-2, 5), fontsize=7, ha="left")
bx.plot([], [], color=C[0], alpha=0.3, lw=6, label="published group-B range")
bx.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), borderaxespad=0.0, ncol=1)
fig.subplots_adjust(left=0.1, right=0.985, top=0.95, bottom=0.36, wspace=0.36)
fig.savefig(out)
