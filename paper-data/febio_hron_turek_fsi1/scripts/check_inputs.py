#!/usr/bin/env python3
"""Reads every .feb input of the package and tabulates the configuration actually written (material, loads, boundary conditions, solver tolerances, mesh counts, point-A node,
root nodes, pins), so the model definition is verified from the inputs themselves.  usage: check_inputs.py INPUT_DIR OUT.csv"""
import sys, os, re, glob, csv, hashlib
import numpy as np
D, OUT = sys.argv[1], sys.argv[2]; rows = []
g = lambda pat, t, d="": (re.search(pat, t, re.S).group(1) if re.search(pat, t, re.S) else d)
for f in sorted(glob.glob(D + "/**/*.feb", recursive=True)):
    t = open(f, encoding="latin-1").read(); nodes = {int(m.group(1)): [float(x) for x in m.group(2).split(",")] for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', t)}
    r = dict(file=os.path.basename(f), sha256=hashlib.sha256(open(f, "rb").read()).hexdigest()[:16], module=g(r'<Module type="(\w+)"', t), nodes=len(nodes),
             elements=len(re.findall(r"<elem id=", t)), element_type=g(r'<Elements type="(\w+)"', t))
    if r["module"] == "solid":
        r.update(material=g(r'<material id="1"[^>]*type="([^"]+)"', t), E=g(r"<E>([^<]+)</E>", t), v=g(r"<v>([^<]+)</v>", t), density=g(r"<density>([^<]+)</density>", t),
                 body_load_y_written=g(r'<body_load type="const">.*?<y lc="1">([^<]+)</y>', t), analysis=g(r"<analysis>(\w+)</analysis>", t), load_steps=g(r"<time_steps>(\d+)</time_steps>", t),
                 solver_tol=f"dtol={g('<dtol>([^<]+)', t)} etol={g('<etol>([^<]+)', t)} rtol={g('<rtol>([^<]+)', t)} qn={g('<qn_method type=.([^/]+)./>', t)}")
        ns = {m.group(1): [int(x) for x in m.group(2).split(",")] for m in re.finditer(r'<NodeSet name="(\w+)">([^<]+)</NodeSet>', t)}
        xs = [nodes[i][0] for i in ns["Root"]]; A = ns["PointA"][0]
        r.update(root_x_min=min(xs), root_x_max=max(xs), n_root_nodes=len(xs), pointA=",".join(f"{v:.6f}" for v in nodes[A]), uz_bc_on_all_nodes=("z_dof>1" in t), clamp_dofs="x,y fixed on Root")
    else:
        r.update(density=g(r"<density>([^<]+)</density>", t), K=g(r"<k>([^<]+)</k>", t), mu=g(r"<mu>([^<]+)</mu>", t), kappa=g(r"<kappa>([^<]+)</kappa>", t), rhoi=g(r"<rhoi>([^<]+)</rhoi>", t),
                 dt=g(r"<step_size>([^<]+)</step_size>", t), steps=g(r"<time_steps>(\d+)</time_steps>", t), inlet_velocity=g(r'type="fluid normal velocity"><velocity lc="1">([^<]+)<', t),
                 parabolic=("<parabolic>1<" in t), rim_pressure=("<prescribe_rim_pressure>1<" in t), outlet=("zero fluid dilatation" in t),
                 pins=";".join(re.findall(r'<bc name="RootPin\w*" node_set="\w+" type="prescribed fluid dilatation"><value lc="1">([^<]+)<', t)) or "none",
                 tolerances=f"vtol={g('<vtol>([^<]+)', t)} ftol={g('<ftol>([^<]+)', t)} etol={g('<etol>([^<]+)', t)} rtol={g('<rtol>([^<]+)', t)} max_refs={g('<max_refs>([^<]+)', t)}",
                 z_velocity_fixed_all_nodes=("wz_dof>1" in t))
    rows.append(r)
keys = []
for r in rows:
    for k in r:
        if k not in keys: keys.append(k)
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, keys); w.writeheader(); w.writerows(rows)
for r in rows: print({k: r.get(k) for k in ("file", "material", "E", "v", "density", "body_load_y_written", "root_x_min", "pointA", "K", "mu", "inlet_velocity", "pins") if r.get(k) is not None})
