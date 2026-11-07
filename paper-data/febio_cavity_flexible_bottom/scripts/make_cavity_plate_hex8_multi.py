#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""cavityFlexibleBottom structural check for the conformal hex8 solid used in the coupled model.
ONE FEBio input holding several INDEPENDENT (non-interacting, 1.5 m apart in y) StVK plates, plane strain, clamped ends, uniform follower pressure
1.315 Pa on the top face, hex8 (full integration), each with its own in-plane divisions (nx, ny). Compares u_y(A) at (4, 0) with the Q9 plane-strain
StVK reference (paper-data plate_laws.py: 20x2 -0.23887, 40x4 -0.24047; Richardson p=2 -> -0.24100 m).
usage: make_cavity_plate_hex8_multi.py --levels 40x4,40x8,80x8,160x16 --out DIR --tag TAG
"""
import argparse, os, json, hashlib
L0, L1, YB, YT, TZ = 2.0, 6.0, -0.4, 0.0, 0.1
E, NU, RHO, P = 500.0, 0.4, 1000.0, 1.315
ap = argparse.ArgumentParser(); ap.add_argument("--levels", required=True); ap.add_argument("--out", required=True); ap.add_argument("--tag", required=True); ap.add_argument("--steps", type=int, default=20); ap.add_argument("--pmap-from", default=None, help="C3/C2 FSI run dir: per-face mean interface pressure (final time) is applied as a face-constant follower pressure (math expression Pmap)")
a = ap.parse_args()
pface = None
if a.pmap_from:
    import re, json as _j
    meta = _j.load(open(os.path.join(a.pmap_from, 'model_meta.json'))); K = meta['k']
    feb = open(os.path.join(a.pmap_from, meta['tag'] + '.feb'), encoding='latin-1').read(); nd = {int(m.group(1)): float(m.group(2).split(',')[0]) for m in re.finditer(r'<node id="(\d+)">([^<]+)</node>', feb)}
    blocks, t = {}, None
    for ln in open(os.path.join(a.pmap_from, 'Mon_iface.txt')):
        if ln.startswith('*Time'): t = float(ln.split('=')[1]); blocks[t] = {}
        elif not ln.startswith('*') and ln.strip(): f = ln.split(); blocks[t][int(f[0])] = -K * float(f[5])
    last = blocks[max(blocks)]; ids = sorted(last, key=lambda i: nd[i]); pn = [last[i] for i in ids]; pface = [0.5 * (pn[k] + pn[k + 1]) for k in range(len(pn) - 1)]
nodes, els, bodies = [], [], []
for b, lv in enumerate(a.levels.split(",")):
    nx, ny = map(int, lv.split("x")); off = -1.5 * b; idx = {}
    def nid(i, j, k):
        if (i, j, k) not in idx:
            nodes.append((L0 + (L1 - L0) * i / nx, YB + (YT - YB) * j / ny + off, TZ * k)); idx[(i, j, k)] = len(nodes)
        return idx[(i, j, k)]
    el0 = len(els); top = []
    for ex in range(nx):
        for ey in range(ny):
            c = lambda aa, bb, kk: nid(ex + aa, ey + bb, kk)
            els.append([c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(0, 1, 0), c(0, 0, 1), c(1, 0, 1), c(1, 1, 1), c(0, 1, 1)])
            if ey == ny - 1: top.append([c(0, 1, 0), c(0, 1, 1), c(1, 1, 1), c(1, 1, 0)])   # normal +y
    ids = list(idx.items())
    bodies.append(dict(lv=lv, nx=nx, ny=ny, off=off, el=(el0 + 1, len(els)), top=top, left=[n for (i, j, k), n in ids if i == 0], right=[n for (i, j, k), n in ids if i == nx],
                       A=[n for (i, j, k), n in ids if i == nx // 2 and j == ny and k == 0][0], topline=[n for (i, j, k), n in ids if j == ny and k == 0]))
lst = lambda l: ",".join(map(str, l)); o = []; w = o.append
w('<?xml version="1.0" encoding="ISO-8859-1"?>'); w(f'<!-- cavity hex8 plate check {a.tag}: independent plates {a.levels} -->'); w('<febio_spec version="4.0">\n\t<Module type="solid"/>')
w(f'\t<Control><analysis>STATIC</analysis><time_steps>{a.steps}</time_steps><step_size>{1.0/a.steps:.12g}</step_size><plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range><plot_level>PLOT_MAJOR_ITRS</plot_level><output_level>OUTPUT_MAJOR_ITRS</output_level>')
w('\t\t<solver type="solid"><symmetric_stiffness>symmetric</symmetric_stiffness><max_refs>15</max_refs><dtol>0.001</dtol><etol>0.01</etol><rtol>0</rtol><lstol>0.9</lstol><min_residual>1e-20</min_residual><qn_method type="BFGS"><max_ups>0</max_ups></qn_method></solver></Control>')
w(f'\t<Material><material id="1" name="Plate" type="St.Venant-Kirchhoff"><density>{RHO:g}</density><E>{E:g}</E><v>{NU:g}</v></material></Material>')
w('\t<Mesh>\n\t\t<Nodes name="AllNodes">')
for n, (x, y, z) in enumerate(nodes, 1): w(f'\t\t\t<node id="{n}">{x:.15g},{y:.15g},{z:.15g}</node>')
w('\t\t</Nodes>')
for b, bd in enumerate(bodies):
    w(f'\t\t<Elements type="hex8" name="P{b}">')
    for i in range(bd["el"][0] - 1, bd["el"][1]): w(f'\t\t\t<elem id="{i+1}">{lst(els[i])}</elem>')
    w('\t\t</Elements>')
w(f'\t\t<NodeSet name="All">{lst(range(1, len(nodes)+1))}</NodeSet>')
for b, bd in enumerate(bodies):
    w(f'\t\t<NodeSet name="L{b}">{lst(bd["left"])}</NodeSet><NodeSet name="R{b}">{lst(bd["right"])}</NodeSet><NodeSet name="A{b}">{bd["A"]}</NodeSet><NodeSet name="T{b}">{lst(bd["topline"])}</NodeSet>')
    w(f'\t\t<Surface name="Top{b}">')
    for i, q in enumerate(bd["top"], 1): w(f'\t\t\t<quad4 id="{i}">{lst(q)}</quad4>')
    w('\t\t</Surface>')
w('\t</Mesh>\n\t<MeshDomains>' + "".join(f'<SolidDomain name="P{b}" mat="Plate"/>' for b in range(len(bodies))) + '</MeshDomains>')
if pface is not None:
    assert len(pface) == len(bodies[0]['top']), (len(pface), len(bodies[0]['top']))
    w('\t<MeshData>\n\t\t<SurfaceData name="Pmap" data_type="scalar" surface="Top0">')
    for k, v in enumerate(pface, 1): w(f'\t\t\t<face lid="{k}">{v:.12g}</face>')
    w('\t\t</SurfaceData>\n\t</MeshData>')
w('\t<Boundary>\n\t\t<bc name="uz" node_set="All" type="zero displacement"><x_dof>0</x_dof><y_dof>0</y_dof><z_dof>1</z_dof></bc>')
for b in range(len(bodies)):
    w(f'\t\t<bc name="cl{b}" node_set="L{b}" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc><bc name="cr{b}" node_set="R{b}" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>')
w('\t</Boundary>\n\t<Loads>')
for b in range(len(bodies)): w(f'\t\t<surface_load name="p{b}" surface="Top{b}" type="pressure"><pressure lc="1"' + (' type="math">Pmap</pressure>' if pface is not None else f'>{P}</pressure>') + f'<symmetric_stiffness>1</symmetric_stiffness><linear>0</linear><shell_bottom>0</shell_bottom></surface_load>')
w('\t</Loads>\n\t<LoadData><load_controller id="1" name="ramp" type="loadcurve"><interpolate>LINEAR</interpolate><points><pt>0,0</pt><pt>1,1</pt></points></load_controller></LoadData>')
w('\t<Output><plotfile type="febio"><var type="displacement"/></plotfile><logfile>')
for b, bd in enumerate(bodies):
    w(f'\t\t\t<node_data data="ux;uy" file="A{b}.txt" delim=" ">{bd["A"]}</node_data>')
w('\t\t</logfile></Output>\n</febio_spec>')
text = "\n".join(o) + "\n"; os.makedirs(a.out, exist_ok=True)
open(os.path.join(a.out, a.tag + ".feb"), "w", encoding="latin-1").write(text)
json.dump(dict(tag=a.tag, levels=a.levels, steps=a.steps, n_nodes=len(nodes), n_elem=len(els), A_nodes=[bd["A"] for bd in bodies], feb_sha256=hashlib.sha256(text.encode("latin-1")).hexdigest(),
               purpose="hex8 StVK plate-only structural check for the conformal FSI solid (independent plates in one run)", E=E, nu=NU, p=P), open(os.path.join(a.out, "model_meta.json"), "w"), indent=1)
print(a.tag, len(nodes), "nodes", len(els), "hex8", [bd["lv"] for bd in bodies])
