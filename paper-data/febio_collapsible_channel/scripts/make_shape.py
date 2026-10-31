#!/usr/bin/env python3
"""Build the CC-F1 frozen-shape file (dx, dy of the interface nodes, final load step of the S1_320x8 static solid run) from Mon_bottom.txt (plain text).
usage: make_shape.py Mon_bottom.txt OUT.json [source-label]   (columns of Mon_bottom.txt: node x y ux uy)"""
import sys, json
rows = []
for l in open(sys.argv[1]):
    if l.startswith("*Time"): rows = []
    elif not l.startswith("*") and l.strip(): rows.append(l.split())
label = sys.argv[3] if len(sys.argv) > 3 else "runs/campaignCC/S1_320x8 final step"
json.dump({"dx": [float(r[3]) for r in rows], "dy": [float(r[4]) for r in rows], "source": label}, open(sys.argv[2], "w"))
