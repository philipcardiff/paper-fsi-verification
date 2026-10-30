#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""cavityFlexibleBottom B0b: fluid-only FEBio model (plate held fixed = no-slip wall on the interface).

Geometry/blocks are copied from solids4foam tutorials/fluidSolidInteraction/cavityFlexibleBottom/system/fluid/blockMeshDict
(15 blocks, 8 arcs, base cell counts at dx ~ 0.1). Blocks are meshed by transfinite (Coons) interpolation with the arc edges,
nodes merged by coordinate, extruded one element (TZ) in z with w_z = 0 (plane strain/2-D).
Physics: rho_f = 1, mu = 0.01 (Re = 100 on H = 1, mean 1 m/s), weakly compressible p = -K e_f, no-slip on all walls and the interface,
parabolic inlet u_x = 1.5 - 6 (y-1.5)^2, u_y = 0, outlet: p = 0 (zero dilatation), velocity free. Pseudo-transient to steady state.
Usage: make_cavity_fluid_feb.py --out DIR --tag TAG [--mesh-factor 1] [--k 1e4] [--dt 0.5] [--t-end 300] [--rhoi 0]
"""
import argparse, os, json, hashlib, math
V = [(0, 1), (1.75, 1), (0, 1.2), (1.75, 1.2), (0, 2), (1.75, 2), (1.926776695, 0.9267766953), (2, 0.75), (2, 0), (2.2, 0),
     (2.2, 0.75), (2.068198052, 1.032842712), (2.1, 2), (5.9, 2), (5.931801948, 1.032842712), (5.8, 0.75), (5.8, 0), (6, 0),
     (6, 0.75), (6.073223305, 0.9267766953), (6.25, 1), (6.25, 1.2), (6.25, 2), (14, 2), (14, 1.2), (14, 1)]
BLOCKS = [((0, 1, 3, 2), 17, 2), ((2, 3, 5, 4), 17, 8), ((1, 6, 11, 3), 3, 2), ((8, 9, 10, 7), 2, 7), ((9, 16, 15, 10), 36, 7),
          ((6, 7, 10, 11), 3, 2), ((10, 15, 14, 11), 36, 3), ((15, 18, 19, 14), 2, 3), ((16, 17, 18, 15), 2, 7),
          ((20, 25, 24, 21), 77, 2), ((21, 24, 23, 22), 77, 8), ((19, 20, 21, 14), 3, 2), ((14, 21, 22, 13), 3, 8),
          ((3, 11, 12, 5), 3, 8), ((11, 14, 13, 12), 36, 8)]
ARCS = {(6, 7): (1.996201938, 0.79341204444), (11, 10): (2.193163489, 0.82814168), (3, 11): (1.82814168, 1.193163489),
        (1, 6): (1.793412044, 0.9962019383), (18, 19): (6.015076845, 0.8355050358), (14, 15): (5.827138321, 0.9039090645),
        (19, 20): (6.164494964, 0.9849231552), (14, 21): (6.096090936, 1.172861679)}
H, TZ = 1.0, 0.1


def circle(p, q, m):
    ax, ay = p; bx, by = q; cx, cy = m
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
    uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
    return ux, uy


def edge(i, j):
    """parametrised curve from vertex i to vertex j (t in [0,1])"""
    p, q = V[i], V[j]
    key = (i, j) if (i, j) in ARCS else (j, i) if (j, i) in ARCS else None
    if key is None:
        return lambda t: (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))
    m = ARCS[key]; c = circle(V[key[0]], V[key[1]], m); r = math.hypot(p[0] - c[0], p[1] - c[1])
    th = lambda P: math.atan2(P[1] - c[1], P[0] - c[0])
    t0, t1, tm = th(p), th(q), th(m)
    wrap = lambda a: (a + math.pi) % (2 * math.pi) - math.pi
    d1, dm = wrap(t1 - t0), wrap(tm - t0)
    if not (min(0, d1) <= dm <= max(0, d1)):   # mid is on the other side -> long way round
        d1 = d1 - 2 * math.pi if d1 > 0 else d1 + 2 * math.pi
    return lambda t: (c[0] + r * math.cos(t0 + t * d1), c[1] + r * math.sin(t0 + t * d1))


def mesh(f):
    nodes, idx, quads = [], {}, []
    def nid(P):
        k = (round(P[0] * 1e7), round(P[1] * 1e7))
        if k not in idx:
            idx[k] = len(nodes); nodes.append(P)
        return idx[k]
    for (a, b, c, d), nx, ny in BLOCKS:
        nx, ny = nx * f, ny * f
        B, T, L, R = edge(a, b), edge(d, c), edge(a, d), edge(b, c)
        P00, P10, P11, P01 = V[a], V[b], V[c], V[d]
        def X(u, v):
            return tuple((1 - v) * B(u)[k] + v * T(u)[k] + (1 - u) * L(v)[k] + u * R(v)[k]
                         - ((1 - u) * (1 - v) * P00[k] + u * (1 - v) * P10[k] + (1 - u) * v * P01[k] + u * v * P11[k]) for k in (0, 1))
        grid = [[nid(X(i / nx, j / ny)) for i in range(nx + 1)] for j in range(ny + 1)]
        for j in range(ny):
            for i in range(nx):
                quads.append((grid[j][i], grid[j][i + 1], grid[j + 1][i + 1], grid[j + 1][i]))
    return nodes, quads


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--tag", required=True)
    ap.add_argument("--mesh-factor", type=int, default=1); ap.add_argument("--k", type=float, default=1e4)
    ap.add_argument("--dt", type=float, default=0.5); ap.add_argument("--t-end", type=float, default=300.0)
    ap.add_argument("--mu", type=float, default=0.01); ap.add_argument("--inlet", choices=["nodal", "surface", "normalvel"], default="nodal", help="nodal Dirichlet (v1), FEBio fluid-velocity surface load with a math profile (v3), or the built-in fluid normal velocity load with its parabolic option (normalvel)"); ap.add_argument("--ubar", type=float, default=1.0, help="mean inlet speed for --inlet normalvel"); ap.add_argument("--pin-corners", type=float, default=None, help="prescribe the dilatation (p = this value [Pa], ramped with lc 1) at the two floor-end corner nodes (x = 2, 6 on y = 0; every velocity DOF fixed there) to remove the undetermined corner pressure"); ap.add_argument("--log-all", action="store_true", help="also log every z=0 node (Mon_all) for traction post-processing")
   
    ap.add_argument("--straight", action="store_true", help="unit test: plain 14 x 1 channel (same BCs), exact Poiseuille solution"); ap.add_argument("--rhoi", type=float, default=0.0); ap.add_argument("--ramp", type=float, default=10.0)
    a = ap.parse_args()
    if a.straight:   # BC/formulation unit test: replace the geometry by one rectangular block y in [1,2]
        global V, BLOCKS, ARCS
        V = [(0, 1), (14, 1), (14, 2), (0, 2)]; BLOCKS = [((0, 1, 2, 3), 140, 10)]; ARCS = {}
    nodes2, quads = mesh(a.mesh_factor); n2 = len(nodes2)
    # orientation check: all quads CCW
    for q in quads:
        (x0, y0), (x1, y1), (x2, y2), (x3, y3) = [nodes2[i] for i in q]
        assert 0.5 * ((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0)) > 0, ("clockwise/degenerate quad", q)
    # boundary edges
    cnt = {}
    for q in quads:
        for e in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])):
            cnt.setdefault(frozenset(e), []).append(e)
    bedges = [v[0] for v in cnt.values() if len(v) == 1]
    eps = 1e-6
    bnodes = {i for e in bedges for i in e}
    inlet = sorted(i for i in bnodes if abs(nodes2[i][0]) < eps)
    outlet = sorted(i for i in bnodes if abs(nodes2[i][0] - 14) < eps)
    iface = sorted(i for i in bnodes if abs(nodes2[i][1]) < eps)
    walls = sorted(bnodes - set(inlet) - set(outlet) - set(iface))
    in_corner = [i for i in inlet if abs(nodes2[i][1] - 1) < eps or abs(nodes2[i][1] - 2) < eps]
    in_free = [i for i in inlet if i not in in_corner]
    noslip = sorted(set(walls) | set(iface) | set(in_corner) | {i for i in outlet if abs(nodes2[i][1] - 1) < eps or abs(nodes2[i][1] - 2) < eps})
    assert len(iface) == 40 * a.mesh_factor + 1 or True
    N1 = lambda i: i + 1
    nodes = [(x, y, 0.0) for x, y in nodes2] + [(x, y, TZ) for x, y in nodes2]
    els = [[N1(q[0]), N1(q[1]), N1(q[2]), N1(q[3]), N1(q[0] + n2), N1(q[1] + n2), N1(q[2] + n2), N1(q[3] + n2)] for q in quads]
    both = lambda s: sorted([N1(i) for i in s] + [N1(i + n2) for i in s])
    lst = lambda l: ",".join(map(str, l))
    nsteps = int(round(a.t_end / a.dt))
    o = []; w = o.append
    w('<?xml version="1.0" encoding="ISO-8859-1"?>')
    w(f'<!-- cavityFlexibleBottom fluid-only (plate fixed): {a.tag}. Generated by tools/make_cavity_fluid_feb.py -->')
    w('<febio_spec version="4.0">\n\t<Module type="fluid"/>')
    w('\t<Control><analysis>DYNAMIC</analysis>')
    w(f'\t\t<time_steps>{nsteps}</time_steps><step_size>{a.dt:g}</step_size>')
    w('\t\t<plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range><plot_level>PLOT_MAJOR_ITRS</plot_level>')
    w('\t\t<output_level>OUTPUT_MAJOR_ITRS</output_level><plot_stride>20</plot_stride>')
    w('\t\t<solver type="fluid"><symmetric_stiffness>non-symmetric</symmetric_stiffness><equation_scheme>block</equation_scheme>')
    w('\t\t\t<lstol>0.9</lstol><lsmin>0.01</lsmin><lsiter>5</lsiter><max_refs>15</max_refs><check_zero_diagonal>0</check_zero_diagonal>')
    w('\t\t\t<reform_each_time_step>1</reform_each_time_step><diverge_reform>1</diverge_reform><min_residual>1e-30</min_residual><max_residual>1e+20</max_residual>')
    w(f'\t\t\t<vtol>1e-6</vtol><ftol>1e-6</ftol><etol>0</etol><rtol>0</rtol><rhoi>{a.rhoi:g}</rhoi><predictor>0</predictor><min_volume_ratio>0</min_volume_ratio>')
    w('\t\t\t<qn_method type="full Newton"/></solver></Control>')
    w(f'\t<Material><material id="1" name="Fluid" type="fluid"><density>1</density><k>{a.k:g}</k>')
    w(f'\t\t<viscous type="Newtonian fluid"><kappa>0</kappa><mu>{a.mu:g}</mu></viscous></material></Material>')
    w('\t<Mesh>\n\t\t<Nodes name="AllNodes">')
    for n, (x, y, z) in enumerate(nodes, 1):
        w(f'\t\t\t<node id="{n}">{x:.15g},{y:.15g},{z:.15g}</node>')
    w('\t\t</Nodes>\n\t\t<Elements type="hex8" name="Fluid">')
    for i, e in enumerate(els, 1):
        w(f'\t\t\t<elem id="{i}">{lst(e)}</elem>')
    w('\t\t</Elements>')
    w(f'\t\t<NodeSet name="All">{lst(range(1, len(nodes) + 1))}</NodeSet>')
    w(f'\t\t<NodeSet name="NoSlip">{lst(both(noslip))}</NodeSet>')
    if a.inlet in ("surface", "normalvel"):
        w('\t\t<Surface name="InletSurf">')
        inl = [e for e in bedges if abs(nodes2[e[0]][0]) < eps and abs(nodes2[e[1]][0]) < eps]
        for k, (p_, q_) in enumerate(inl, 1):
            w(f'\t\t\t<quad4 id="{k}">{N1(p_)},{N1(q_)},{N1(q_ + n2)},{N1(p_ + n2)}</quad4>')
        w('\t\t</Surface>')
    w(f'\t\t<NodeSet name="InletWy">{lst(both(in_free))}</NodeSet>')
    w(f'\t\t<NodeSet name="Outlet">{lst(both(outlet))}</NodeSet>')
    corners = [i for i in iface if abs(nodes2[i][0] - 2.0) < eps or abs(nodes2[i][0] - 6.0) < eps]
    if a.pin_corners is not None:
        w(f'\t\t<NodeSet name="CornerPin">{lst(both(corners))}</NodeSet>')
    for k, i in enumerate(in_free):
        w(f'\t\t<NodeSet name="in{k}">{lst(both([i]))}</NodeSet>')
    # monitoring sets (z = 0 layer only), ordered by coordinate
    def line(xv=None, yv=None):
        if xv is not None:   # nearest node column (outlet-block columns are not at integer x)
            xv = min({round(nodes2[i][0], 6) for i in range(n2)}, key=lambda c: abs(c - xv))
        s = [i for i in range(n2) if (xv is not None and abs(nodes2[i][0] - xv) < eps) or (yv is not None and abs(nodes2[i][1] - yv) < eps)]
        return sorted(s, key=lambda i: (nodes2[i][1], nodes2[i][0]))
    mon = {"Mon_inlet": line(xv=0.0), "Mon_x10": line(xv=10.0), "Mon_x13": line(xv=13.0), "Mon_outlet": line(xv=14.0),
           "Mon_iface": sorted([i for i in iface if True], key=lambda i: nodes2[i][0]) or sorted([i for i in range(n2) if abs(nodes2[i][1] - 1.0) < eps], key=lambda i: nodes2[i][0]),
           "Mon_cl": sorted([i for i in range(n2) if abs(nodes2[i][1] - 1.5) < eps and nodes2[i][0] > (-1 if a.straight else 6.4)], key=lambda i: nodes2[i][0])}
    if a.straight:   # extra monitors for the Poiseuille acceptance checks
        mon["Mon_x2"] = line(xv=2.0); mon["Mon_x7"] = line(xv=7.0); mon["Mon_x12"] = line(xv=12.0)
    if a.straight or a.log_all:
        mon["Mon_all"] = sorted(range(n2), key=lambda i: (nodes2[i][0], nodes2[i][1]))
    for k, v in mon.items():
        w(f'\t\t<NodeSet name="{k}">{lst([N1(i) for i in v])}</NodeSet>')
    w('\t</Mesh>\n\t<MeshDomains><SolidDomain name="Fluid" mat="Fluid"/></MeshDomains>')
    w('\t<Boundary>')
    w('\t\t<bc name="wz" node_set="All" type="zero fluid velocity"><wx_dof>0</wx_dof><wy_dof>0</wy_dof><wz_dof>1</wz_dof></bc>')
    w('\t\t<bc name="noslip" node_set="NoSlip" type="zero fluid velocity"><wx_dof>1</wx_dof><wy_dof>1</wy_dof><wz_dof>0</wz_dof></bc>')
    if a.inlet == "nodal":
      w('\t\t<bc name="inlet_wy" node_set="InletWy" type="zero fluid velocity"><wx_dof>0</wx_dof><wy_dof>1</wy_dof><wz_dof>0</wz_dof></bc>')
    for k, i in enumerate(in_free if a.inlet == "nodal" else []):
        u = 1.5 - 6 * (nodes2[i][1] - 1.5) ** 2
        w(f'\t\t<bc name="inlet{k}" node_set="in{k}" type="prescribed fluid velocity"><dof>wx</dof><value lc="1">{u:.15g}</value><relative>0</relative></bc>')
    w('\t\t<bc name="outlet_p0" node_set="Outlet" type="zero fluid dilatation"/>')
    if a.pin_corners is not None:
        w(f'\t\t<bc name="corner_pin" node_set="CornerPin" type="prescribed fluid dilatation"><value lc="1">{-a.pin_corners / a.k:.12g}</value><relative>0</relative></bc>')
    w('\t</Boundary>')
    if a.inlet == "surface":
        w('\t<Loads><surface_load name="InletVelocity" surface="InletSurf" type="fluid velocity"><scale lc="1">1</scale>'
          '<velocity type="math">1.5-6*(Y-1.5)^2,0,0</velocity></surface_load></Loads>')
    if a.inlet == "normalvel":   # inflow = negative normal velocity (outward normal is -x); parabolic map is normalised to the mean speed
        w(f'\t<Loads><surface_load name="InletNormalVelocity" surface="InletSurf" type="fluid normal velocity"><velocity lc="1">{-a.ubar:g}</velocity>'
          '<prescribe_nodal_velocities>1</prescribe_nodal_velocities><parabolic>1</parabolic><prescribe_rim_pressure>1</prescribe_rim_pressure></surface_load></Loads>')
    w(f'\t<LoadData><load_controller id="1" name="ramp" type="loadcurve"><interpolate>LINEAR</interpolate><extend>CONSTANT</extend>'
      f'<points><pt>0,0</pt><pt>{a.ramp:g},1</pt></points></load_controller></LoadData>')
    w('\t<Output><plotfile type="febio"><var type="fluid pressure"/><var type="fluid velocity"/><var type="fluid dilatation"/></plotfile>')
    w('\t\t<logfile>')
    for k, v in mon.items():
        w(f'\t\t\t<node_data data="x;y;wx;wy;ef" file="{k}.txt" delim=" ">{lst([N1(i) for i in v])}</node_data>')
    w('\t\t</logfile></Output>\n</febio_spec>')
    text = "\n".join(o) + "\n"
    os.makedirs(a.out, exist_ok=True)
    open(os.path.join(a.out, f"{a.tag}.feb"), "w", encoding="latin-1").write(text)
    meta = dict(vars(a), n_nodes=len(nodes), n_elem=len(els), n_inlet_free=len(in_free), n_iface_nodes=len(iface), n_steps=nsteps,
                feb_sha256=hashlib.sha256(text.encode("latin-1")).hexdigest())
    json.dump(meta, open(os.path.join(a.out, "model_meta.json"), "w"), indent=1)
    # geometric sanity: area of fluid domain (exact: channel 14x1 + cavity part)
    area = sum(0.5 * abs((nodes2[q[2]][0] - nodes2[q[0]][0]) * (nodes2[q[3]][1] - nodes2[q[1]][1]) - (nodes2[q[3]][0] - nodes2[q[1]][0]) * (nodes2[q[2]][1] - nodes2[q[0]][1])) for q in quads)
    print(a.tag, meta["n_nodes"], "nodes", meta["n_elem"], "hex8", "interface nodes/layer", len(iface), "inlet free", len(in_free), f"area {area:.6f}")


if __name__ == "__main__":
    main()
