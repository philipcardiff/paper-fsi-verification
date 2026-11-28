#!/usr/bin/env python3
"""Point-A histories: axisymmetric levels against the 3-D levels.

Usage: plot_histories.py --axi L1=a.csv --axi L2=b.csv ... --td L1=x.csv ...
                         --out fig.png
CSV columns: time_s, radial_displacement_m, axial_displacement_m, ...
"""
import argparse
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def load(path):
    t, r, a = [], [], []
    for row in csv.DictReader(open(path)):
        t.append(float(row["time_s"]) * 1e3)
        r.append(float(row["radial_displacement_m"]) * 1e3)
        a.append(float(row["axial_displacement_m"]) * 1e3)
    return t, r, a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--axi", action="append", default=[])
    ap.add_argument("--td", action="append", default=[])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    blues = plt.cm.Blues
    greys = plt.cm.Oranges
    for i, spec in enumerate(args.axi):
        lab, path = spec.split("=", 1)
        t, r, a = load(path)
        c = blues(0.35 + 0.6 * i / max(1, len(args.axi) - 1))
        for ax, y in ((axes[0], r), (axes[2], a)):
            ax.plot(t, y, color=c, lw=1.2, label=f"axi {lab}")
        axes[1].plot(t, r, color=c, lw=1.2, label=f"axi {lab}")
    for i, spec in enumerate(args.td):
        lab, path = spec.split("=", 1)
        t, r, a = load(path)
        c = greys(0.4 + 0.5 * i / max(1, len(args.td) - 1))
        for ax, y in ((axes[0], r), (axes[2], a)):
            ax.plot(t, y, color=c, lw=1.0, ls="--", label=f"3-D {lab}")
        axes[1].plot(t, r, color=c, lw=1.0, ls="--", label=f"3-D {lab}")
    axes[0].set(xlabel="t (ms)", ylabel="u_r(A) (mm)", title="radial")
    axes[1].set(xlabel="t (ms)", ylabel="u_r(A) (mm)", title="radial peak",
                xlim=(6.9, 7.6), ylim=(0.1565, 0.1605))
    axes[2].set(xlabel="t (ms)", ylabel="u_z(A) (mm)", title="axial")
    for ax in axes:
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)


if __name__ == "__main__":
    main()
