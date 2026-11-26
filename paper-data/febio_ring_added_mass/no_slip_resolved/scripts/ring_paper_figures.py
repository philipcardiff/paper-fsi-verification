#!/usr/bin/env python3
"""Paper versions of the FEBio ringAddedMass figures, in the shared paper style
(paper-data/figstyle), drawn at the final printed width, as vector PDFs:

  R1  fig:ring_febio_ratio    wet/dry frequency-ratio errors: (a) FEBio no-slip
      against A and V; (b) solids4foam IQN-ILS (partitioned) against A, mesh 1, 2, 4
  R3  fig:ring_febio_decay    FEBio f1 displacement histories with exp(-alpha_V t)
  R4  fig:ring_febio_traction interface traction amplitudes; ring energy against
      traction work (strong level, potential-flow start)

The data are those of ring_figure_R1_partitioned.py (R1) and ring_figures.py
(R3, R4), which still write the package's PNG versions; nothing is recomputed
except what those scripts compute. Stored results only, no simulation.

Usage: python3 ring_paper_figures.py OUT_DIR [--pkg PKG_DIR]
  OUT_DIR: where ringAddedMass_febio_R{1,3,4}.pdf are written
           (the paper uses <repo>/figures/verification)
  PKG_DIR: this package (default: the parent of this script's directory)
"""
import argparse, csv, json, math, os, sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("out_dir")
ap.add_argument("--pkg", default=str(HERE.parent))
a = ap.parse_args()
PKG = Path(a.pkg).resolve()
PARENT = PKG.parent
REPO_DATA = PARENT.parent                      # <repo>/paper-data
sys.path.insert(0, str(REPO_DATA / "figstyle"))
sys.path.insert(0, str(PARENT / "tools"))
import figstyle                                  # noqa: E402
import analyze_ring as AR                        # noqa: E402
import matplotlib.pyplot as plt                  # noqa: E402

figstyle.use()
C = figstyle.OKABE_ITO
OUT = Path(a.out_dir); OUT.mkdir(parents=True, exist_ok=True)
W = figstyle.TEXTWIDTH

J = json.load(open(PARENT / "reference/s4f_ringAddedMass/verification/reference/"
                   "ringAddedMass_verification_references.json"))["exact"]["levels"]
LV = ["weak", "moderate", "strong"]
tab = {r["id"]: r for r in json.load(open(PKG / "results/ring_final_table.json"))}


def refs(l):
    v = J[l]; s = v["viscousOmega"]; z = v["viscousDampingRatio"]
    return dict(A=v["ratio"], alpha=z * s / math.sqrt(1 - z * z), omega=v["omega"])


# ---------------------------------------------------------------- R1
s4f = {}
for x in csv.DictReader(open(PKG / "reference/ringAddedMass_monolithic_mesh_ff150e4.csv")):
    n = x["name"]
    if n.startswith("wet_") and "iqnils" in n:
        s4f.setdefault(n.split("_")[1], {})[int(x["factor"])] = 100 * float(x["ratio_error"])
for l in LV:
    assert sorted(s4f[l]) == [1, 2, 4], l

fig, (ax, bx) = plt.subplots(1, 2, figsize=(W, 2.35), gridspec_kw=dict(width_ratios=[1.25, 1]))
x = np.arange(3); w = 0.19
bars = [("1", "errA_pct", C[0], None, "f1, vs A"), ("1", "errV_pct", "white", "////", "f1, vs V"),
        ("2", "errA_pct", C[1], None, "f2, vs A"), ("2", "errV_pct", "white", "\\\\\\\\", "f2, vs V")]
for k, (rid, key, fc, hatch, lab) in enumerate(bars):
    ec = C[0] if rid == "1" else C[1]
    ax.bar(x + (k - 1.5) * w, [tab[f"P{l[0]}{rid}"][key] for l in LV], w * 0.92, color=fc,
           edgecolor=ec, hatch=hatch, lw=0.6, label=lab)
for i, l in enumerate(LV):
    ax.plot([i - 2 * w, i + 2 * w], [100 * J[l]["viscousFrequencyShift"]] * 2, "k:", lw=1.0)
ax.plot([], [], "k:", lw=1.0, label="V relative to A")
ax.axhline(0, color="0.5", lw=0.5)
ax.set_xticks(x)
ax.set_xticklabels(["weak\n$m^*=0.101$", "moderate\n$m^*=1.008$", "strong\n$m^*=10.08$"])
ax.tick_params(axis="x", which="both", top=False, bottom=False)
ax.set_ylabel("wet/dry ratio error (%)")
ax.legend(loc="lower left", ncol=1, handlelength=1.6, borderaxespad=0.3)
figstyle.panel_label(ax, "(a)")

mk = {1: ("o", "-", C[0]), 2: ("s", "--", C[1]), 4: ("^", ":", C[2])}
ms = [J[l]["massRatio"] for l in LV]
for fac, (m, ls, c) in mk.items():
    bx.plot(ms, [s4f[l][fac] for l in LV], marker=m, ls=ls, color=c, mfc="white" if fac == 2 else c,
            label=f"mesh {fac}")
bx.set_xscale("log"); bx.axhline(0, color="0.5", lw=0.5)
bx.set_xlabel("mass ratio $m^*$")
bx.set_ylabel("ratio error vs A (%)")
bx.legend(loc="best")
figstyle.panel_label(bx, "(b)")
fig.tight_layout(w_pad=1.5)
fig.savefig(OUT / "ringAddedMass_febio_R1.pdf"); plt.close(fig)

# ---------------------------------------------------------------- R3
O = PKG / "results"
fig, axs = plt.subplots(1, 3, figsize=(W, 1.75), sharey=True)
for i, (l, ax) in enumerate(zip(LV, axs)):
    r = refs(l); T = 2 * math.pi / r["omega"]
    t, d, v = AR.series(str(O / f"runs/P{l[0]}1/probes.txt"))
    k0 = int(np.argmax(np.abs(d[: int(len(d) / 4)]))); a0 = np.abs(d[k0])
    ax.plot(t / T, d / a0, color=C[0], lw=0.5, label="FEBio f1")
    tt = np.linspace(t[k0], t[-1], 100)
    env = np.exp(-r["alpha"] * (tt - t[k0]))
    ax.plot(tt / T, env, "k--", lw=0.9, label=r"$\pm\exp(-\alpha_V t)$")
    ax.plot(tt / T, -env, "k--", lw=0.9)
    ax.set_xlabel("$t/T$")
    ax.set_ylim(-2.45, 2.05)
    ax.text(0.97, 0.93, f"{l}: fit ${tab[f'P{l[0]}1']['alpha_over_V']:.3f}\\,\\alpha_V$",
            transform=ax.transAxes, ha="right", va="top")
    figstyle.panel_label(ax, f"({'abc'[i]})")
axs[0].set_ylabel("$D(t)$ / first peak")
axs[0].legend(loc="lower left", borderaxespad=0.2, handlelength=1.6, labelspacing=0.2)
fig.tight_layout(w_pad=0.8)
fig.savefig(OUT / "ringAddedMass_febio_R3.pdf"); plt.close(fig)

# ---------------------------------------------------------------- R4
fig, (ax, bx) = plt.subplots(1, 2, figsize=(W, 2.35), gridspec_kw=dict(width_ratios=[1.3, 1]))
tj = {k: json.load(open(O / f"analysis/{k}.traction.json")) for k in ("Pw1", "Pm1", "Ps1", "Ds")}
names = ["Pw1", "Pm1", "Ps1", "Ds"]
lab = ["weak\n(rest)", "moderate\n(rest)", "strong\n(rest)", "strong\n(pot. flow)"]
x = np.arange(4); w = 0.26
for k, (key, c, h, lb) in enumerate([("pressure_traction_n2_amp_Pa", C[0], None, "pressure"),
                                      ("viscous_tangential_traction_n2_amp_Pa", C[1], "////", "viscous shear"),
                                      ("thin_layer_shear_estimate_Pa", C[2], "....", "thin-layer estimate")]):
    ax.bar(x + (k - 1) * w, [tj[n][key] for n in names], w * 0.92, color=c, hatch=h, edgecolor="white",
           lw=0.0, label=lb)
ax.set_yscale("log"); ax.set_ylim(1e-5, 1e3); ax.set_xticks(x); ax.set_xticklabels(lab)
ax.tick_params(axis="x", which="both", top=False, bottom=False)
ax.set_ylabel("$n=2$ traction amplitude (Pa)")
ax.legend(loc="upper left", borderaxespad=0.3, handlelength=1.4)
figstyle.panel_label(ax, "(a)")
z = np.load(O / "analysis/Ds.traction.npz"); e = np.load(O / "analysis/Ds.energy.npz")
Er = e["KEr"] + e["PEr"]
S = 1e4                                          # plot in units of 1e-4 J/m
bx.plot(z["t"], S * (Er - Er[0]), color=C[0], lw=1.4, label=r"ring energy change $\Delta E$")
bx.plot(z["t"], S * z["Wtr"], "--", color=C[1], lw=1.0, label="traction work $W$")
bx.plot(z["t"], S * z["Wv"] * 50, ":", color="k", lw=1.0, label=r"viscous part of $W$, $\times 50$")
bx.set_xlabel("$t$ (s)"); bx.set_ylabel(r"energy per unit depth ($10^{-4}$ J/m)")
bx.set_ylim(-0.5, 2.1)
bx.legend(loc="upper left", borderaxespad=0.3, handlelength=1.8)
figstyle.panel_label(bx, "(b)")
fig.tight_layout(w_pad=1.5)
fig.savefig(OUT / "ringAddedMass_febio_R4.pdf"); plt.close(fig)
print("wrote", *(OUT / f"ringAddedMass_febio_R{i}.pdf" for i in (1, 3, 4)))
