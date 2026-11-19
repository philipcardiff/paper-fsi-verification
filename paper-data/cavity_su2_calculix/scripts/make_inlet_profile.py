#!/usr/bin/env python3
"""SU2 inlet profile for cavityFlexibleBottom: parabolic, peak UMAX = 1.5 m/s (mean 1 m/s),
channel y0 = 1 to y0 + H = 2 m (solids4foam tutorial: parabolicVelocity, maxValue 1.5, n (1 0 0)).

Writes the SU2 7.5 INC velocity-inlet profile format (x, y, T, |V|, nx, ny) at every node of
the 'inlet' marker of a 2-D .su2 mesh (SU2 is vertex-based: values at nodes).

Usage: make_inlet_profile.py MESH.su2 OUT.dat [--umax 1.5 --y0 1.0 --height 1.0]
"""
import argparse


def inlet_nodes(mesh):
    L = open(mesh).read().split("\n")
    pts, nodes = None, set()
    i = 0
    while i < len(L):
        s = L[i]
        if s.startswith("NPOIN="):
            n = int(s.split("=")[1])
            pts = [tuple(map(float, L[i + 1 + k].split()[:2])) for k in range(n)]
            i += n
        elif s.startswith("MARKER_TAG=") and s.split("=")[1].strip() == "inlet":
            n = int(L[i + 1].split("=")[1])
            for k in range(n):
                nodes.update(map(int, L[i + 2 + k].split()[1:3]))
            i += n + 1
        i += 1
    return sorted((pts[j] for j in nodes), key=lambda p: p[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("out")
    ap.add_argument("--umax", type=float, default=1.5)
    ap.add_argument("--y0", type=float, default=1.0)
    ap.add_argument("--height", type=float, default=1.0)
    a = ap.parse_args()
    H, U, y0 = a.height, a.umax, a.y0
    prof = lambda y: max(0.0, 4.0 * U * (y - y0) * (y0 + H - y) / H ** 2)
    P = inlet_nodes(a.mesh)
    with open(a.out, "w") as f:
        f.write(f"NMARK= 1\nMARKER_TAG= inlet\nNROW={len(P)}\nNCOL=6\n")
        f.write("# COORD-X COORD-Y TEMPERATURE VELOCITY NORMAL-X NORMAL-Y\n")
        for x, y in P:
            f.write(f"{x:.15e}\t{y:.15e}\t2.881500000000000e+02\t{prof(y):.15e}\t1.0\t0.0\n")
    ys = [p[1] for p in P]
    q = sum(0.5 * (prof(ys[k]) + prof(ys[k + 1])) * (ys[k + 1] - ys[k]) for k in range(len(ys) - 1))
    print(f"{a.out}: {len(P)} inlet nodes y [{ys[0]}, {ys[-1]}], Umax {U}, trapezoidal mean {q / H:.6f} m/s (exact {2 * U / 3:.6f})")


if __name__ == "__main__":
    main()
