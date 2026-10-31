#!/usr/bin/env python3
"""Hron-Turek CFD1-type rigid-obstacle fluid model (cylinder AND flag rigid, FSI1 inlet/outlet formulation) for FEBio v4.13 (fluid module, fixed mesh).  Stage C of the FSI1 pilot.

Geometry/topology: the 24 blocks, 10 arcs and grading of solids4foam tutorials/fluidSolidInteraction/TurekHron2006/system/fluid/blockMeshDict (copy with sha256 in
tools/hronturek_data/, commit b200718); the fluid domain already excludes cylinder and flag, so the obstacle wall = cylinder + flag (both rigid).  Base cell counts x --factor
(rounded identically on shared edges); blocks meshed by transfinite (Coons) interpolation with arcs and simple grading (OpenFOAM expansion ratio = last/first cell), nodes merged by
coordinate, extruded one hex8 element (TZ = 0.015) in z with w_z = 0 (2-D).
Physics: rho = 1000, mu = rho nu = 1.0 Pa s (nu = 1e-3), weakly compressible p = -K e_f (K = --k), no-slip on cylinder/flag/top/bottom, inlet: built-in `fluid normal velocity` (parabolic map,
rim pressure; the nodal map is normalised so the trapezoid mean of the nodal speeds equals the load value 0.2 m/s => discrete flux exactly 0.2 x 0.41 g(t)), g(t) = (1-cos(pi t/2))/2 for
t <= 2 s then 1 (load curve exact at the step times), outlet: zero fluid dilatation (p = 0) on all outlet nodes, velocity natural.  Pseudo-transient to the steady state, dt constant.
Usage: make_hronturek_cfd_feb.py --out DIR --tag TAG --factor F [--k 1e6] [--dt 0.5] [--t-end 40] [--ramp 2] [--pin-root PA] [--relaxed] [--rhoi 0]
"""
import argparse, os, json, hashlib, math, re, sys

BMD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hronturek_data", "TurekHron2006_fluid_blockMeshDict.s4f_b200718")
TZ, HCH, LCH, MU, RHO, UBAR = 0.015, 0.41, 2.5, 1.0, 1000.0, 0.2
XROOT = 0.2 + math.sqrt(0.05**2 - 0.01**2)


def parse_bmd():
    txt = open(BMD).read()
    sc = float(re.search(r"scale\s+([0-9.eE+-]+);", txt).group(1))
    vs = txt[txt.index("vertices"):txt.index("blocks")]
    pts = [tuple(float(x) * sc for x in m.groups()) for m in re.finditer(r"\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", vs)]
    blocks = []
    for m in re.finditer(r"hex\s*\(([\d\s]+)\)\s*\((\d+)\s+(\d+)\s+(\d+)\)\s*simpleGrading\s*\(\s*([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s*\)", txt):
        ids = [int(x) for x in m.group(1).split()]
        blocks.append((tuple(ids[:4]), int(m.group(2)), int(m.group(3)), float(m.group(5)), float(m.group(6))))
    arcs = {}
    for m in re.finditer(r"arc\s+(\d+)\s+(\d+)\s*\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", txt):
        if abs(float(m.group(5))) < 1e-9:
            arcs[(int(m.group(1)), int(m.group(2)))] = (float(m.group(3)) * sc, float(m.group(4)) * sc)
    return pts, blocks, arcs


def circle(p, q, m):
    ax, ay = p; bx, by = q; cx, cy = m
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
    uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
    return ux, uy


def make_edge(V, ARCS, i, j):
    p, q = V[i][:2], V[j][:2]
    key = (i, j) if (i, j) in ARCS else (j, i) if (j, i) in ARCS else None
    if key is None:
        return lambda t: (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))
    m = ARCS[key]; c = circle(V[key[0]][:2], V[key[1]][:2], m); r = math.hypot(p[0] - c[0], p[1] - c[1])
    th = lambda P: math.atan2(P[1] - c[1], P[0] - c[0]); t0, t1, tm = th(p), th(q), th(m)
    wrap = lambda a: (a + math.pi) % (2 * math.pi) - math.pi
    d1, dm = wrap(t1 - t0), wrap(tm - t0)
    if not (min(0, d1) <= dm <= max(0, d1)):
        d1 = d1 - 2 * math.pi if d1 > 0 else d1 + 2 * math.pi
    return lambda t: (c[0] + r * math.cos(t0 + t * d1), c[1] + r * math.sin(t0 + t * d1))


def grade(n, r):
    if n == 1 or abs(r - 1) < 1e-12:
        return [i / n for i in range(n + 1)]
    q = r ** (1.0 / (n - 1))
    return [(1 - q**i) / (1 - q**n) for i in range(n + 1)]


def mesh(factor):
    V, BL, ARCS = parse_bmd()
    nodes, idx, quads, qblock = [], {}, [], []
    def nid(P):
        k = (round(P[0] * 1e8), round(P[1] * 1e8))
        if k not in idx:
            idx[k] = len(nodes); nodes.append(P)
        return idx[k]
    for bi, ((a, b, c, d), nx, ny, gx, gy) in enumerate(BL):
        nx, ny = max(1, int(round(nx * factor))), max(1, int(round(ny * factor)))
        B, T, L, R = make_edge(V, ARCS, a, b), make_edge(V, ARCS, d, c), make_edge(V, ARCS, a, d), make_edge(V, ARCS, b, c)
        P00, P10, P11, P01 = V[a][:2], V[b][:2], V[c][:2], V[d][:2]
        us, vs = grade(nx, gx), grade(ny, gy)
        def X(u, v):
            return tuple((1 - v) * B(u)[k] + v * T(u)[k] + (1 - u) * L(v)[k] + u * R(v)[k]
                         - ((1 - u) * (1 - v) * P00[k] + u * (1 - v) * P10[k] + (1 - u) * v * P01[k] + u * v * P11[k]) for k in (0, 1))
        grid = [[nid(X(us[i], vs[j])) for i in range(nx + 1)] for j in range(ny + 1)]
        for j in range(ny):
            for i in range(nx):
                quads.append((grid[j][i], grid[j][i + 1], grid[j + 1][i + 1], grid[j + 1][i])); qblock.append(bi)
    return nodes, quads, qblock


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--tag", required=True)
    ap.add_argument("--factor", type=float, default=1.0); ap.add_argument("--k", type=float, default=1e6)
    ap.add_argument("--dt", type=float, default=0.5); ap.add_argument("--t-end", type=float, default=40.0); ap.add_argument("--ramp", type=float, default=2.0)
    ap.add_argument("--rhoi", type=float, default=0.0); ap.add_argument("--pin-root", type=float, default=None, help="prescribe p [Pa] (dilatation -p/K) at the two flag-root junction nodes (both z layers)")
    ap.add_argument("--pin-root-values", default=None, help="two pressures [Pa] for the root junction nodes: first for the node at y = 0.19, second for y = 0.21 (declared rule: mean of the two wall neighbours in the unpinned run)")
    ap.add_argument("--relaxed", action="store_true", help="declared fallback tolerances (Womersley B1b): vtol=ftol=dtol=0, rtol=1e-3, etol=1e-4")
    ap.add_argument("--purpose", default="Hron-Turek CFD1 rigid-flag fluid-only (stage C)")
    a = ap.parse_args()
    nodes2, quads, qblock = mesh(a.factor); n2 = len(nodes2)
    # --- quality / topology checks
    def cross(p0, p1, p2): return (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0])
    minsj, maxar, area = 1e9, 0.0, 0.0
    for q in quads:
        P = [nodes2[i] for i in q]
        for k in range(4):
            p0, p1, p2 = P[k - 1], P[k], P[(k + 1) % 4]
            e1 = (p0[0] - p1[0], p0[1] - p1[1]); e2 = (p2[0] - p1[0], p2[1] - p1[1])
            sj = (e2[0] * e1[1] - e2[1] * e1[0]) / (math.hypot(*e1) * math.hypot(*e2))   # sin of corner angle, signed (CCW positive)
            # for a CCW quad the signed cross of (next - cur) x (prev - cur) is positive
            minsj = min(minsj, sj)
        ed = [math.dist(P[k], P[(k + 1) % 4]) for k in range(4)]; maxar = max(maxar, max(ed) / min(ed))
        area += 0.5 * ((P[2][0] - P[0][0]) * (P[3][1] - P[1][1]) - (P[3][0] - P[1][0]) * (P[2][1] - P[0][1]))
    assert minsj > 0.05, ("non-positive or degenerate quad corner, min sin(angle)", minsj)
    exact_area = LCH * HCH - math.pi * 0.05**2 - (0.6 - XROOT) * 0.02
    cnt = {}
    for q in quads:
        for e in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])):
            cnt.setdefault(frozenset(e), []).append(e)
    assert all(len(v) <= 2 for v in cnt.values()), "edge shared by more than 2 quads (merge failure)"
    bedges = [v[0] for v in cnt.values() if len(v) == 1]
    eps = 1e-7
    bnodes = {i for e in bedges for i in e}
    inlet = sorted(i for i in bnodes if abs(nodes2[i][0]) < eps)
    outlet = sorted(i for i in bnodes if abs(nodes2[i][0] - LCH) < eps)
    bottom = sorted(i for i in bnodes if abs(nodes2[i][1]) < eps)
    top = sorted(i for i in bnodes if abs(nodes2[i][1] - HCH) < eps)
    obst = sorted(bnodes - set(inlet) - set(outlet) - set(bottom) - set(top))
    # every boundary edge must lie on exactly one declared boundary
    def which(e):
        i, j = e
        if i in inlet and j in inlet: return "inlet"
        if i in outlet and j in outlet: return "outlet"
        if i in bottom and j in bottom: return "bottom"
        if i in top and j in top: return "top"
        if i in obst and j in obst: return "obstacle"
        return "UNCLASSIFIED"
    kinds = [which(e) for e in bedges]
    assert "UNCLASSIFIED" not in kinds, "boundary edge not on a declared boundary"
    in_free = [i for i in inlet if 1e-7 < nodes2[i][1] < HCH - 1e-7]
    in_corner = [i for i in inlet if i not in in_free]
    out_corner = [i for i in outlet if i in bottom or i in top]
    noslip = sorted(set(obst) | set(bottom) | set(top) | set(in_corner) | set(out_corner))
    root = [i for i in obst if abs(nodes2[i][0] - XROOT) < 1e-6 and (abs(nodes2[i][1] - 0.19) < 1e-6 or abs(nodes2[i][1] - 0.21) < 1e-6)]
    assert len(root) == 2, ("flag root junction nodes", root)
    # wall-adjacent element height on the obstacle (first-cell wall distance)
    hw = []
    for q in quads:
        for k in range(4):
            i, j = q[k], q[(k + 1) % 4]
            if which((i, j)) == "obstacle" and frozenset((i, j)) in cnt and len(cnt[frozenset((i, j))]) == 1:
                mid = ((nodes2[i][0] + nodes2[j][0]) / 2, (nodes2[i][1] + nodes2[j][1]) / 2)
                o1, o2 = q[(k + 3) % 4], q[(k + 2) % 4]
                om = ((nodes2[o1][0] + nodes2[o2][0]) / 2, (nodes2[o1][1] + nodes2[o2][1]) / 2)
                hw.append(math.dist(mid, om))
    hw.sort()
    N1 = lambda i: i + 1
    nodes = [(x, y, 0.0) for x, y in nodes2] + [(x, y, TZ) for x, y in nodes2]
    els = [[N1(q[0]), N1(q[1]), N1(q[2]), N1(q[3]), N1(q[0] + n2), N1(q[1] + n2), N1(q[2] + n2), N1(q[3] + n2)] for q in quads]
    both = lambda s: sorted([N1(i) for i in s] + [N1(i + n2) for i in s])
    lst = lambda l: ",".join(map(str, l))
    nsteps = int(round(a.t_end / a.dt))
    o = []; w = o.append
    w('<?xml version="1.0" encoding="ISO-8859-1"?>')
    w(f'<!-- Hron-Turek CFD1 rigid flag fluid-only: {a.tag}. Generated by tools/make_hronturek_cfd_feb.py v1 -->')
    w('<febio_spec version="4.0">\n\t<Module type="fluid"/>')
    w('\t<Control><analysis>DYNAMIC</analysis>')
    w(f'\t\t<time_steps>{nsteps}</time_steps><step_size>{a.dt:g}</step_size>')
    w('\t\t<plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range><plot_level>PLOT_MAJOR_ITRS</plot_level>')
    w('\t\t<output_level>OUTPUT_MAJOR_ITRS</output_level><plot_stride>20</plot_stride>')
    w('\t\t<solver type="fluid"><symmetric_stiffness>non-symmetric</symmetric_stiffness><equation_scheme>block</equation_scheme>')
    w('\t\t\t<lstol>0.9</lstol><lsmin>0.01</lsmin><lsiter>5</lsiter><max_refs>15</max_refs><check_zero_diagonal>0</check_zero_diagonal>')
    w('\t\t\t<reform_each_time_step>1</reform_each_time_step><diverge_reform>1</diverge_reform><min_residual>1e-30</min_residual><max_residual>1e+20</max_residual>')
    if a.relaxed:
        w(f'\t\t\t<vtol>0</vtol><ftol>0</ftol><dtol>0</dtol><etol>1e-4</etol><rtol>1e-3</rtol><rhoi>{a.rhoi:g}</rhoi><predictor>0</predictor><min_volume_ratio>0</min_volume_ratio>')
    else:
        w(f'\t\t\t<vtol>1e-6</vtol><ftol>1e-6</ftol><etol>0</etol><rtol>0</rtol><rhoi>{a.rhoi:g}</rhoi><predictor>0</predictor><min_volume_ratio>0</min_volume_ratio>')
    w('\t\t\t<qn_method type="full Newton"/></solver></Control>')
    w(f'\t<Material><material id="1" name="Fluid" type="fluid"><density>{RHO:g}</density><k>{a.k:g}</k>')
    w(f'\t\t<viscous type="Newtonian fluid"><kappa>0</kappa><mu>{MU:g}</mu></viscous></material></Material>')
    w('\t<Mesh>\n\t\t<Nodes name="AllNodes">')
    for n, (x, y, z) in enumerate(nodes, 1):
        w(f'\t\t\t<node id="{n}">{x:.15g},{y:.15g},{z:.15g}</node>')
    w('\t\t</Nodes>\n\t\t<Elements type="hex8" name="Fluid">')
    for i, e in enumerate(els, 1):
        w(f'\t\t\t<elem id="{i}">{lst(e)}</elem>')
    w('\t\t</Elements>')
    w(f'\t\t<NodeSet name="All">{lst(range(1, len(nodes) + 1))}</NodeSet>')
    w(f'\t\t<NodeSet name="NoSlip">{lst(both(noslip))}</NodeSet>')
    w('\t\t<Surface name="InletSurf">')
    inl = [e for e, kd in zip(bedges, kinds) if kd == "inlet"]
    for k, (p_, q_) in enumerate(inl, 1):
        w(f'\t\t\t<quad4 id="{k}">{N1(p_)},{N1(q_)},{N1(q_ + n2)},{N1(p_ + n2)}</quad4>')
    w('\t\t</Surface>')
    w(f'\t\t<NodeSet name="Outlet">{lst(both(outlet))}</NodeSet>')
    pinv = None
    if a.pin_root_values is not None:
        pa, pb = [float(x) for x in a.pin_root_values.split(",")]
        pinv = {min(root, key=lambda i: nodes2[i][1]): pa, max(root, key=lambda i: nodes2[i][1]): pb}
    elif a.pin_root is not None:
        pinv = {i: a.pin_root for i in root}
    if pinv is not None:
        for nm_, i_ in (("RootPinLow", min(root, key=lambda i: nodes2[i][1])), ("RootPinHigh", max(root, key=lambda i: nodes2[i][1]))):
            w(f'\t\t<NodeSet name="{nm_}">{lst(both([i_]))}</NodeSet>')
    w(f'\t\t<NodeSet name="Mon_all">{lst([N1(i) for i in range(n2)])}</NodeSet>')
    w('\t</Mesh>\n\t<MeshDomains><SolidDomain name="Fluid" mat="Fluid"/></MeshDomains>')
    w('\t<Boundary>')
    w('\t\t<bc name="wz" node_set="All" type="zero fluid velocity"><wx_dof>0</wx_dof><wy_dof>0</wy_dof><wz_dof>1</wz_dof></bc>')
    w('\t\t<bc name="noslip" node_set="NoSlip" type="zero fluid velocity"><wx_dof>1</wx_dof><wy_dof>1</wy_dof><wz_dof>0</wz_dof></bc>')
    w('\t\t<bc name="outlet_p0" node_set="Outlet" type="zero fluid dilatation"/>')
    if pinv is not None:
        for nm_, i_ in (("RootPinLow", min(root, key=lambda i: nodes2[i][1])), ("RootPinHigh", max(root, key=lambda i: nodes2[i][1]))):
            w(f'\t\t<bc name="{nm_}" node_set="{nm_}" type="prescribed fluid dilatation"><value lc="1">{-pinv[i_] / a.k:.12g}</value><relative>0</relative></bc>')
    w('\t</Boundary>')
    w(f'\t<Loads><surface_load name="InletNormalVelocity" surface="InletSurf" type="fluid normal velocity"><velocity lc="1">{-UBAR:g}</velocity>'
      '<prescribe_nodal_velocities>1</prescribe_nodal_velocities><parabolic>1</parabolic><prescribe_rim_pressure>1</prescribe_rim_pressure></surface_load></Loads>')
    pts = []
    kk = 0
    while kk * a.dt < a.ramp - 1e-9:
        t = kk * a.dt; pts.append((t, 0.5 * (1 - math.cos(math.pi * t / a.ramp)))); kk += 1
    pts.append((a.ramp, 1.0)); pts.append((a.t_end, 1.0))
    w('\t<LoadData><load_controller id="1" name="ramp" type="loadcurve"><interpolate>LINEAR</interpolate><extend>CONSTANT</extend><points>'
      + "".join(f"<pt>{t:.12g},{g:.15g}</pt>" for t, g in pts) + '</points></load_controller></LoadData>')
    w('\t<Output><plotfile type="febio"><var type="fluid pressure"/><var type="fluid velocity"/><var type="fluid dilatation"/></plotfile>')
    w('\t\t<logfile>')
    w(f'\t\t\t<node_data data="wx;wy;ef" file="Mon_all.txt" delim=" ">{lst([N1(i) for i in range(n2)])}</node_data>')
    w('\t\t</logfile></Output>\n</febio_spec>')
    text = "\n".join(o) + "\n"
    os.makedirs(a.out, exist_ok=True)
    open(os.path.join(a.out, f"{a.tag}.feb"), "w", encoding="latin-1").write(text)
    import numpy as np
    np.savez(os.path.join(a.out, "mesh2d.npz"), xy=np.array(nodes2), quads=np.array(quads), n2=n2, tz=TZ, inlet=np.array(inlet), outlet=np.array(outlet),
             bottom=np.array(bottom), top=np.array(top), obst=np.array(obst), root=np.array(root))
    meshtxt = "".join(o[o.index('\t<Mesh>\n\t\t<Nodes name="AllNodes">') if False else 0:0])
    mesh_sha = hashlib.sha256(("".join(f"{x:.15g},{y:.15g}\n" for x, y in nodes2) + "".join(f"{q}\n" for q in quads)).encode()).hexdigest()
    obst_edges = [e for e, kd in zip(bedges, kinds) if kd == "obstacle"]
    meta = dict(vars(a), n_nodes=len(nodes), n_elem=len(els), n_steps=nsteps, n_inlet_free=len(in_free), n_obstacle_nodes_layer=len(obst), n_obstacle_edges=len(obst_edges),
                area=area, area_exact=exact_area, area_rel_diff=(area - exact_area) / exact_area, min_sin_corner=minsj, max_edge_aspect=maxar,
                first_cell_wall_dist_min=hw[0], first_cell_wall_dist_median=hw[len(hw) // 2], first_cell_wall_dist_max=hw[-1], pin_nodes_layer=root,
                inlet_dy_max=max(abs(nodes2[e[0]][1] - nodes2[e[1]][1]) for e in inl), n_inlet_edges=len(inl),
                feb_sha256=hashlib.sha256(text.encode("latin-1")).hexdigest(), mesh_sha256=mesh_sha, bmd_sha256=hashlib.sha256(open(BMD, "rb").read()).hexdigest(),
                generator="tools/make_hronturek_cfd_feb.py v2 (v1 + --pin-root-values)", pin_values_applied=(None if pinv is None else {str(k): v for k, v in pinv.items()}), purpose=a.purpose)
    json.dump(meta, open(os.path.join(a.out, "model_meta.json"), "w"), indent=1, default=float)
    print(a.tag, f"cells={len(els)} nodes={len(nodes)} area_rel_diff={meta['area_rel_diff']:.2e} min_sin={minsj:.3f} max_aspect={maxar:.1f} wall_h(min/med/max)={hw[0]*1e3:.2f}/{hw[len(hw)//2]*1e3:.2f}/{hw[-1]*1e3:.2f} mm obst_edges={len(obst_edges)} inlet_edges={len(inl)}")


if __name__ == "__main__":
    main()
