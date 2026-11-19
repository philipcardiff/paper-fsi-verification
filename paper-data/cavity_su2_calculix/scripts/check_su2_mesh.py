#!/usr/bin/env python3
"""Sanity checks for a 2-D .su2 quad mesh: element areas, marker extents, cylinder radius."""
import sys
import numpy as np

def read_su2(p):
    L = open(p).read().split("\n"); i = 0; el = []; pts = []; mk = {}
    while i < len(L):
        s = L[i]
        if s.startswith("NELEM="):
            n = int(s.split("=")[1]); el = [list(map(int, L[i+1+k].split()[1:5])) for k in range(n)]; i += n
        elif s.startswith("NPOIN="):
            n = int(s.split("=")[1]); pts = np.array([list(map(float, L[i+1+k].split()[:2])) for k in range(n)]); i += n
        elif s.startswith("MARKER_TAG="):
            name = s.split("=")[1].strip(); n = int(L[i+1].split("=")[1])
            mk[name] = np.array([list(map(int, L[i+2+k].split()[1:3])) for k in range(n)]); i += n + 1
        i += 1
    return np.array(el), pts, mk

el, X, mk = read_su2(sys.argv[1])
x, y = X[el, 0], X[el, 1]
A = 0.5 * (np.sum(x * np.roll(y, -1, 1), 1) - np.sum(np.roll(x, -1, 1) * y, 1))
print(f"elements {len(el)}: area min {A.min():.3e} max {A.max():.3e} total {A.sum():.6f} m^2 (positive: {(A > 0).all()})")
print(f"expected fluid area 2.5*0.41 - pi*0.05^2 - plate {2.5*0.41 - np.pi*0.05**2 - (0.6-0.2489898)*0.02:.6f} (plate approx, cylinder arcs are polygonal)")
for name, e in mk.items():
    n = np.unique(e)
    print(f"{name:9s} {len(e):4d} edges  x [{X[n,0].min():.5f}, {X[n,0].max():.5f}]  y [{X[n,1].min():.5f}, {X[n,1].max():.5f}]")
r = np.hypot(X[np.unique(mk['cylinder']), 0] - 0.2, X[np.unique(mk['cylinder']), 1] - 0.2)
print(f"cylinder radius: min {r.min():.6f} max {r.max():.6f} (0.05)")
inl = np.unique(mk['inlet']); print(f"inlet nodes {len(inl)}, x = {np.unique(X[inl,0].round(9))}")
