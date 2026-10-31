#!/usr/bin/env python3
"""collapsibleChannel (solids4foam tutorial, commit b200718) FEBio v4.13 generator: three modes.  v1 (2026-10-31, COLLAPSIBLE-CHANNEL pilot).
Geometry (CI, tutorials/fluidSolidInteraction/collapsibleChannel/system/*/blockMeshDict): channel y in [0,1], x in [0,25]; rigid upstream 0<x<5, flexible TOP wall 5<x<15 occupying
y in [1, 1.05] (thickness h = 0.05), rigid downstream 15<x<25.  Fluid (tutorial): 40 + 80 + 80 columns x 16 uniform rows = 3,200 cells; wall 80 x 8.  One-element slab in z (TZ = 0.1), plane strain: u_z = w_z = 0.
Modes
  solid : STATIC solid-only wall (SVK, E = 4.55e8, nu = 0.3, plane strain), clamped end faces (all nodes at x = 5, 15 fixed), follower `pressure` p_ext = 200 Pa on the TOP face, 20 equal load
          increments (the tutorial 'static' study).  --nsc x --ns hex8 (or hex20) cells.  Log: bottom-face (y = 1, interface) nodes ux, uy (final step is the result).
  fluid : frozen-geometry fluid-only (fluid module): the fluid mesh is the (optionally deformed) channel; no-slip on walls and interface; pressure-driven ends as the tutorial (prescribed fluid dilatation
          e_f = -p/K, p = 6 inlet / 0 outlet kinematic = Pa at rho = 1, wy = 0 at the end faces, wx natural); exact-Poiseuille initial condition (per-column/row constants: the v4.13 math-IC crash).
  fsi   : coupled monolithic ALE fluid-FSI: conformal hex8 wall (nsc = fluid wall columns) sharing the interface nodes, as the cavity/ring models; follower p_ext on the top face with the tutorial
          ramp 200(1-cos(pi t/0.25))/2 for t <= 0.25 s (tabulated every 0.00125 s as the tutorial's constant/solid/pExt), fluid-FSI traction on the interface, mesh regulariser neo-Hookean E_mesh.
Physics: rho_f = 1, nu = 0.02 (mu = 0.02), K = --k; wall rho_s = 1.  Shape (fluid mode): --shape FILE.json with 'dx','dy' lists per wall interface node (from a solid run), columns follow the interface x, rows scale in y.
usage: make_cc_feb.py --mode {solid,fluid,fsi} --out DIR --tag TAG [options]  (see --help)
"""
import argparse, os, json, hashlib, math, sys

LUP, LW, LDN, HCH, HW, TZ = 5.0, 10.0, 10.0, 1.0, 0.05, 0.1
E_S, NU_S, RHO_S, RHO_F, MU_F = 4.55e8, 0.3, 1.0, 1.0, 0.02
P_IN, P_OUT, P_EXT, T_RAMP = 6.0, 0.0, 200.0, 0.25


def pext_table(tend):
    pts = [(k * 0.00125, 0.5 * (1 - math.cos(math.pi * k * 0.00125 / T_RAMP))) for k in range(0, 201)]
    pts.append((tend, 1.0)) if tend > T_RAMP else None
    return pts


def pois_p(x):   # kinematic = Pa at rho = 1: p = 12 nu (25 - x)
    return 12 * MU_F / RHO_F * (LUP + LW + LDN - x) * RHO_F


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["solid", "fluid", "fsi"], required=True); ap.add_argument("--out", required=True); ap.add_argument("--tag", required=True)
    ap.add_argument("--nu", type=int, default=40); ap.add_argument("--nw", type=int, default=80); ap.add_argument("--nd", type=int, default=80); ap.add_argument("--nr", type=int, default=16)
    ap.add_argument("--nsc", type=int, default=None, help="solid columns (default nw)"); ap.add_argument("--ns", type=int, default=8, help="solid rows through the thickness")
    ap.add_argument("--elem", choices=["hex8", "hex20"], default="hex8", help="solid element (solid mode only)")
    ap.add_argument("--k", type=float, default=1e7); ap.add_argument("--dt", type=float, default=0.025); ap.add_argument("--t-end", type=float, default=3.5)
    ap.add_argument("--rhoi", type=float, default=0.9); ap.add_argument("--e-mesh", type=float, default=1.0); ap.add_argument("--max-refs", type=int, default=25)
    ap.add_argument("--dtol", type=float, default=1e-5); ap.add_argument("--vtol", type=float, default=1e-5); ap.add_argument("--ftol", type=float, default=1e-5)
    ap.add_argument("--etol", type=float, default=0.0); ap.add_argument("--rtol", type=float, default=0.0); ap.add_argument("--min-res", type=float, default=1e-20)
    ap.add_argument("--steps", type=int, default=20, help="solid mode load steps"); ap.add_argument("--shape", default=None, help="fluid mode: JSON with dx, dy per wall interface node")
    ap.add_argument("--pins", default=None, help="fluid mode: comma list of junction pin values 'lo,hi' [Pa] (x = 5 and x = 15 interface-corner nodes)")
    ap.add_argument("--plot-stride", type=int, default=0); ap.add_argument("--rows-logged", type=int, default=3, help="fluid rows (from the wall) logged for every x")
    ap.add_argument("--purpose", default="")
    a = ap.parse_args()
    nsc = a.nsc if a.nsc else a.nw
    lst = lambda l: ",".join(map(str, l)); eps = 1e-9
    o = []; w = o.append
    tag = a.tag
    ctrl_fsi = (f'\t<Control><analysis>DYNAMIC</analysis><time_steps>{int(round(a.t_end / a.dt))}</time_steps><step_size>{a.dt:.12g}</step_size>'
                f'<plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range><plot_level>PLOT_MAJOR_ITRS</plot_level><output_level>OUTPUT_MAJOR_ITRS</output_level>'
                f'<plot_stride>{a.plot_stride or 40}</plot_stride>\n')
    solver_fsi = ('\t\t<solver type="fluid-FSI"><symmetric_stiffness>non-symmetric</symmetric_stiffness><equation_scheme>staggered</equation_scheme><equation_order>default</equation_order><optimize_bw>0</optimize_bw>'
                  f'<lstol>0.9</lstol><lsmin>0.01</lsmin><lsiter>5</lsiter><max_refs>{a.max_refs}</max_refs><check_zero_diagonal>0</check_zero_diagonal><reform_each_time_step>1</reform_each_time_step><diverge_reform>1</diverge_reform>'
                  f'<min_residual>{a.min_res:g}</min_residual><max_residual>1e+20</max_residual><dtol>{a.dtol:g}</dtol><vtol>{a.vtol:g}</vtol><ftol>{a.ftol:g}</ftol><etol>{a.etol:g}</etol><rtol>{a.rtol:g}</rtol>'
                  f'<rhoi>{a.rhoi:g}</rhoi><predictor>0</predictor><min_volume_ratio>0</min_volume_ratio><order>2</order><qn_method type="full Newton"/></solver></Control>')
    meta = dict(vars(a), generator="tools/make_cc_feb.py v1", source_commit="b200718f4f8390c90f9c7e3ded8b74f3e13941cc", TZ=TZ)

    # ------------------------------------------------------------------ solid mode
    if a.mode == "solid":
        nx, ny = nsc, a.ns
        nodes = []; idx = {}
        order20 = a.elem == "hex20"
        def nid(i, j, k):   # grid index on the (2nx+1) x (2ny+1) x 3 lattice for hex20, (nx+1)x(ny+1)x2 for hex8
            if (i, j, k) not in idx:
                if order20: x = LUP + LW * i / (2 * nx); y = HCH + HW * j / (2 * ny); z = TZ * k / 2
                else: x = LUP + LW * i / nx; y = HCH + HW * j / ny; z = TZ * k
                nodes.append((x, y, z)); idx[(i, j, k)] = len(nodes)
            return idx[(i, j, k)]
        els = []; topf = []; botf = []
        for ex in range(nx):
            for ey in range(ny):
                if not order20:
                    c = lambda da, db, k: nid(ex + da, ey + db, k)
                    els.append([c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(0, 1, 0), c(0, 0, 1), c(1, 0, 1), c(1, 1, 1), c(0, 1, 1)])
                    if ey == ny - 1: topf.append([c(0, 1, 0), c(0, 1, 1), c(1, 1, 1), c(1, 1, 0)])
                else:
                    i0, j0 = 2 * ex, 2 * ey; c = lambda da, db, k: nid(i0 + da, j0 + db, k)
                    corners = [c(0, 0, 0), c(2, 0, 0), c(2, 2, 0), c(0, 2, 0), c(0, 0, 2), c(2, 0, 2), c(2, 2, 2), c(0, 2, 2)]
                    mb = [c(1, 0, 0), c(2, 1, 0), c(1, 2, 0), c(0, 1, 0)]; mt = [c(1, 0, 2), c(2, 1, 2), c(1, 2, 2), c(0, 1, 2)]
                    vert = [c(0, 0, 1), c(2, 0, 1), c(2, 2, 1), c(0, 2, 1)]
                    els.append(corners + mb + mt + vert)
                    if ey == ny - 1: topf.append([c(0, 2, 0), c(0, 2, 2), c(2, 2, 2), c(2, 2, 0), c(0, 2, 1), c(1, 2, 2), c(2, 2, 1), c(1, 2, 0)])
        N = nodes
        left = [n for n, (x, y, z) in enumerate(N, 1) if abs(x - LUP) < eps]; right = [n for n, (x, y, z) in enumerate(N, 1) if abs(x - LUP - LW) < eps]
        bot = sorted([n for n, (x, y, z) in enumerate(N, 1) if abs(y - HCH) < eps and abs(z) < eps], key=lambda n: N[n - 1][0])
        w('<?xml version="1.0" encoding="ISO-8859-1"?>'); w(f'<!-- collapsibleChannel static wall, {tag}. Generated by tools/make_cc_feb.py v1 -->')
        w('<febio_spec version="4.0">\n\t<Module type="solid"/>\n\t<Control><analysis>STATIC</analysis>')
        w(f'\t\t<time_steps>{a.steps}</time_steps><step_size>{1.0 / a.steps:.12g}</step_size><plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range><plot_level>PLOT_MAJOR_ITRS</plot_level><output_level>OUTPUT_MAJOR_ITRS</output_level>')
        w('\t\t<solver type="solid"><symmetric_stiffness>non-symmetric</symmetric_stiffness><max_refs>25</max_refs><dtol>1e-9</dtol><etol>1e-14</etol><rtol>0</rtol><lstol>0.9</lstol><min_residual>1e-30</min_residual><qn_method type="full Newton"/></solver>\n\t</Control>')
        w(f'\t<Material><material id="1" name="Wall" type="St.Venant-Kirchhoff"><density>{RHO_S:g}</density><E>{E_S:g}</E><v>{NU_S:g}</v></material></Material>')
        w('\t<Mesh>\n\t\t<Nodes name="AllNodes">')
        for n, (x, y, z) in enumerate(N, 1): w(f'\t\t\t<node id="{n}">{x:.15g},{y:.15g},{z:.15g}</node>')
        w(f'\t\t</Nodes>\n\t\t<Elements type="{a.elem}" name="Wall">')
        for i, e in enumerate(els, 1): w(f'\t\t\t<elem id="{i}">{lst(e)}</elem>')
        w('\t\t</Elements>')
        w(f'\t\t<NodeSet name="All">{lst(range(1, len(N) + 1))}</NodeSet><NodeSet name="ClampL">{lst(left)}</NodeSet><NodeSet name="ClampR">{lst(right)}</NodeSet><NodeSet name="Bottom">{lst(bot)}</NodeSet>')
        w('\t\t<Surface name="Top">')
        for k, q in enumerate(topf, 1): w(f'\t\t\t<{"quad4" if not order20 else "quad8"} id="{k}">{lst(q)}</{"quad4" if not order20 else "quad8"}>')
        w('\t\t</Surface>\n\t</Mesh>\n\t<MeshDomains><SolidDomain name="Wall" mat="Wall"/></MeshDomains>\n\t<Boundary>')
        w('\t\t<bc name="uz" node_set="All" type="zero displacement"><x_dof>0</x_dof><y_dof>0</y_dof><z_dof>1</z_dof></bc>')
        w('\t\t<bc name="clampL" node_set="ClampL" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>')
        w('\t\t<bc name="clampR" node_set="ClampR" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>\n\t</Boundary>')
        w(f'\t<Loads><surface_load name="pext" surface="Top" type="pressure"><pressure lc="1">{P_EXT:g}</pressure><symmetric_stiffness>0</symmetric_stiffness><linear>0</linear><shell_bottom>0</shell_bottom></surface_load></Loads>')
        w('\t<LoadData><load_controller id="1" name="ramp" type="loadcurve"><interpolate>LINEAR</interpolate><points><pt>0,0</pt><pt>1,1</pt></points></load_controller></LoadData>')
        w('\t<Output><plotfile type="febio"><var type="displacement"/><var type="stress"/></plotfile>\n\t\t<logfile>')
        w(f'\t\t\t<node_data data="x;y;ux;uy" file="Mon_bottom.txt" delim=" ">{lst(bot)}</node_data>')
        w('\t\t</logfile></Output>\n</febio_spec>')
        text = "\n".join(o) + "\n"
        meta.update(n_nodes=len(N), n_elem=len(els), aspect=(LW / nx) / (HW / ny), bottom_nodes=bot, n_top_faces=len(topf))
    else:
        # ------------------------------------------------------------------ fluid / fsi: shared fluid mesh
        nu_, nw_, nd_, nr = a.nu, a.nw, a.nd, a.nr
        xs = [LUP * i / nu_ for i in range(nu_ + 1)] + [LUP + LW * i / nw_ for i in range(1, nw_ + 1)] + [LUP + LW + LDN * i / nd_ for i in range(1, nd_ + 1)]
        NX = len(xs) - 1; i0w = nu_; i1w = nu_ + nw_
        F = lambda i, j: i * (nr + 1) + j
        n2 = (NX + 1) * (nr + 1)
        dxs = [0.0] * (NX + 1); dys = [0.0] * (NX + 1)
        if a.shape:
            sh = json.load(open(a.shape)); assert len(sh["dx"]) == nw_ + 1, (len(sh["dx"]), nw_)
            for q in range(nw_ + 1): dxs[i0w + q] = sh["dx"][q]; dys[i0w + q] = sh["dy"][q]
            meta["shape_sha256"] = hashlib.sha256(open(a.shape, "rb").read()).hexdigest()
        nodes2 = [None] * n2
        for i in range(NX + 1):
            for j in range(nr + 1): nodes2[F(i, j)] = (xs[i] + dxs[i], (j / nr) * (HCH + dys[i]))
        quads = [[F(i, j), F(i + 1, j), F(i + 1, j + 1), F(i, j + 1)] for i in range(NX) for j in range(nr)]
        for q in quads:
            (x0, y0), (x1, y1), (x2, y2), (x3, y3) = [nodes2[m] for m in q]; assert 0.5 * ((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0)) > 0, "negative cell"
        nodes = [(x, y, 0.0) for x, y in nodes2] + [(x, y, TZ) for x, y in nodes2]
        N1 = lambda i: i + 1
        both = lambda s: sorted([N1(i) for i in s] + [N1(i + n2) for i in s])
        fl_els = [[N1(q[0]), N1(q[1]), N1(q[2]), N1(q[3]), N1(q[0] + n2), N1(q[1] + n2), N1(q[2] + n2), N1(q[3] + n2)] for q in quads]
        inlet = [F(0, j) for j in range(nr + 1)]; outlet = [F(NX, j) for j in range(nr + 1)]
        bottom_w = [F(i, 0) for i in range(NX + 1)]
        topw = [F(i, nr) for i in range(NX + 1) if i <= i0w or i >= i1w]       # rigid top wall incl. the two junction nodes
        iface = [F(i, nr) for i in range(i0w, i1w + 1)]
        junc = [F(i0w, nr), F(i1w, nr)]
        inlet_c = [F(0, 0), F(0, nr)]; outlet_c = [F(NX, 0), F(NX, nr)]
        noslip_f = set(bottom_w) | set(topw) | set(iface) | set(inlet_c) | set(outlet_c)
        rowkeep = a.rows_logged
        mon_rows = [F(i, nr - r) for i in range(NX + 1) for r in range(rowkeep)] + [F(i, 0) for i in range(NX + 1)]
        mon_rows = sorted(set(mon_rows))
        # exact Poiseuille initial condition (columns for dilatation, rows for velocity)
        sets_ef = {}; sets_v = {}
        for i in range(NX + 1): sets_ef[f"ef{i}"] = (-pois_p(xs[i]) / a.k, both([F(i, j) for j in range(nr + 1)]))
        for j in range(nr + 1): sets_v[f"wv{j}"] = (6.0 * (j / nr) * (1 - j / nr), both([F(i, j) for i in range(NX + 1)]))
        meta.update(n_fluid_elem=len(fl_els), NX=NX, nodes_per_layer=n2, interface_nodes_layer=iface, junction_nodes_layer=junc, interface_x=[xs[i] for i in range(i0w, i1w + 1)])
        w('<?xml version="1.0" encoding="ISO-8859-1"?>')
        w(f'<!-- collapsibleChannel {a.mode}, {tag}. Generated by tools/make_cc_feb.py v1 -->')
        if a.mode == "fluid":
            w('<febio_spec version="4.0">\n\t<Module type="fluid"/>')
            nst = int(round(a.t_end / a.dt))
            w(f'\t<Control><analysis>DYNAMIC</analysis><time_steps>{nst}</time_steps><step_size>{a.dt:.12g}</step_size><plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range><plot_level>PLOT_MAJOR_ITRS</plot_level><output_level>OUTPUT_MAJOR_ITRS</output_level><plot_stride>{a.plot_stride or 1000}</plot_stride>')
            w('\t\t<solver type="fluid"><symmetric_stiffness>non-symmetric</symmetric_stiffness><equation_scheme>staggered</equation_scheme><equation_order>default</equation_order><optimize_bw>0</optimize_bw>')
            w(f'\t\t\t<lstol>0.9</lstol><lsmin>0.01</lsmin><lsiter>5</lsiter><max_refs>15</max_refs><check_zero_diagonal>0</check_zero_diagonal><reform_each_time_step>1</reform_each_time_step><diverge_reform>1</diverge_reform><min_residual>1e-30</min_residual><max_residual>1e+20</max_residual>')
            w(f'\t\t\t<vtol>1e-6</vtol><ftol>1e-6</ftol><etol>0</etol><rtol>0</rtol><rhoi>{a.rhoi:g}</rhoi><predictor>0</predictor><min_volume_ratio>0</min_volume_ratio><qn_method type="full Newton"/></solver></Control>')
            w(f'\t<Material><material id="1" name="Fluid" type="fluid"><density>{RHO_F:g}</density><k>{a.k:g}</k><viscous type="Newtonian fluid"><kappa>0</kappa><mu>{MU_F:g}</mu></viscous></material></Material>')
        else:
            w('<febio_spec version="4.0">\n\t<Module type="fluid-FSI"/>'); w(ctrl_fsi.rstrip("\n")); w(solver_fsi)
            w('\t<Material>')
            w(f'\t\t<material id="1" name="Wall" type="St.Venant-Kirchhoff"><density>{RHO_S:g}</density><E>{E_S:g}</E><v>{NU_S:g}</v></material>')
            w('\t\t<material id="2" name="Fluid" type="fluid-FSI">')
            w(f'\t\t\t<solid type="neo-Hookean"><density>0</density><E>{a.e_mesh:g}</E><v>0</v></solid>')
            w(f'\t\t\t<fluid type="fluid"><density>{RHO_F:g}</density><k>{a.k:g}</k><viscous type="Newtonian fluid"><kappa>0</kappa><mu>{MU_F:g}</mu></viscous></fluid></material>\n\t</Material>')
        w('\t<Mesh>\n\t\t<Nodes name="AllNodes">')
        # solid nodes (fsi): bottom row shared with fluid interface nodes (requires nsc == nw)
        pl_els = []; top = []; sol_nodes = []
        if a.mode == "fsi":
            assert nsc == nw_, "conformal interface needs solid columns == fluid wall columns"
            pidx = {}; ns = a.ns
            def pn(i, j, k):
                if j == 0: return N1(iface[i] + (n2 if k else 0))
                if (i, j, k) not in pidx:
                    nodes.append((LUP + LW * i / nsc, HCH + HW * j / ns, TZ * k)); pidx[(i, j, k)] = len(nodes)
                return pidx[(i, j, k)]
            for ex in range(nsc):
                for ey in range(ns):
                    c = lambda aa, bb, kk: pn(ex + aa, ey + bb, kk)
                    pl_els.append([c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(0, 1, 0), c(0, 0, 1), c(1, 0, 1), c(1, 1, 1), c(0, 1, 1)])
                    if ey == ns - 1: top.append([c(0, 1, 0), c(0, 1, 1), c(1, 1, 1), c(1, 1, 0)])
            sol_nodes = sorted(set(n for e in pl_els for n in e))
        for n, (x, y, z) in enumerate(nodes, 1): w(f'\t\t\t<node id="{n}">{x:.15g},{y:.15g},{z:.15g}</node>')
        w('\t\t</Nodes>')
        eid = 0
        if pl_els:
            w('\t\t<Elements type="hex8" name="Wall">')
            for e in pl_els: eid += 1; w(f'\t\t\t<elem id="{eid}">{lst(e)}</elem>')
            w('\t\t</Elements>')
        w('\t\t<Elements type="hex8" name="Fluid">')
        for e in fl_els: eid += 1; w(f'\t\t\t<elem id="{eid}">{lst(e)}</elem>')
        w('\t\t</Elements>')
        w(f'\t\t<NodeSet name="AllN">{lst(range(1, len(nodes) + 1))}</NodeSet><NodeSet name="FluidAll">{lst(range(1, 2 * n2 + 1))}</NodeSet>')
        w(f'\t\t<NodeSet name="NoSlipW">{lst(sorted(both(noslip_f)))}</NodeSet>')
        endnodes = sorted(set(inlet) | set(outlet)) ; endfree = [i for i in endnodes if i not in set(inlet_c + outlet_c)]
        w(f'\t\t<NodeSet name="EndWy">{lst(both(endfree))}</NodeSet><NodeSet name="InletP">{lst(both(inlet))}</NodeSet><NodeSet name="OutletP">{lst(both(outlet))}</NodeSet>')
        for k_, (v_, ids_) in sets_ef.items(): w(f'\t\t<NodeSet name="{k_}">{lst(ids_)}</NodeSet>')
        for k_, (v_, ids_) in sets_v.items(): w(f'\t\t<NodeSet name="{k_}">{lst(ids_)}</NodeSet>')
        if a.pins:
            w(f'\t\t<NodeSet name="PinLo">{lst(both([junc[0]]))}</NodeSet><NodeSet name="PinHi">{lst(both([junc[1]]))}</NodeSet>')
        if a.mode == "fsi":
            fixmesh = sorted(set(both(set(bottom_w) | set(topw) | set(inlet) | set(outlet))) - set(both(iface)))
            fixmesh = sorted(set(both(set(bottom_w) | set([i for i in topw if i not in junc]) | set(inlet) | set(outlet))))
            left = [n for n in sol_nodes if abs(nodes[n - 1][0] - LUP) < eps]; right = [n for n in sol_nodes if abs(nodes[n - 1][0] - LUP - LW) < eps]
            w(f'\t\t<NodeSet name="FixMesh">{lst(fixmesh)}</NodeSet><NodeSet name="ClampL">{lst(left)}</NodeSet><NodeSet name="ClampR">{lst(right)}</NodeSet>')
            w('\t\t<Surface name="Top">')
            for k, q in enumerate(top, 1): w(f'\t\t\t<quad4 id="{k}">{lst(q)}</quad4>')
            w('\t\t</Surface>\n\t\t<Surface name="Interface">')
            for k, (p_, q_) in enumerate(zip(iface[:-1], iface[1:]), 1): w(f'\t\t\t<quad4 id="{k}">{N1(p_)},{N1(q_)},{N1(q_ + n2)},{N1(p_ + n2)}</quad4>')
            w('\t\t</Surface>')
        mons = {"Mon_iface": ("ux;uy;wx;wy;ef" if a.mode == "fsi" else "wx;wy;ef", [N1(i) for i in iface]),
                "Mon_rows": ("wx;wy;ef", [N1(i) for i in mon_rows]),
                "Mon_inlet": ("wx;wy;ef", [N1(i) for i in inlet]), "Mon_outlet": ("wx;wy;ef", [N1(i) for i in outlet])}
        for k_, (dat, ids_) in mons.items(): w(f'\t\t<NodeSet name="{k_}">{lst(ids_)}</NodeSet>')
        w('\t</Mesh>')
        w('\t<MeshDomains>' + ('<SolidDomain name="Wall" mat="Wall"/><SolidDomain name="Fluid" mat="Fluid"/>' if a.mode == "fsi" else '<SolidDomain name="Fluid" mat="Fluid"/>') + '</MeshDomains>')
        w('\t<Boundary>')
        if a.mode == "fsi":
            w('\t\t<bc name="uz" node_set="AllN" type="zero displacement"><x_dof>0</x_dof><y_dof>0</y_dof><z_dof>1</z_dof></bc>')
            w('\t\t<bc name="clampL" node_set="ClampL" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>')
            w('\t\t<bc name="clampR" node_set="ClampR" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>')
            w('\t\t<bc name="fixmesh" node_set="FixMesh" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>')
        w('\t\t<bc name="wz" node_set="FluidAll" type="zero fluid velocity"><wx_dof>0</wx_dof><wy_dof>0</wy_dof><wz_dof>1</wz_dof></bc>')
        w('\t\t<bc name="noslip" node_set="NoSlipW" type="zero fluid velocity"><wx_dof>1</wx_dof><wy_dof>1</wy_dof><wz_dof>0</wz_dof></bc>')
        w('\t\t<bc name="endwy" node_set="EndWy" type="zero fluid velocity"><wx_dof>0</wx_dof><wy_dof>1</wy_dof><wz_dof>0</wz_dof></bc>')
        w(f'\t\t<bc name="p_in" node_set="InletP" type="prescribed fluid dilatation"><value lc="2">{-P_IN / a.k:.15g}</value><relative>0</relative></bc>')
        w(f'\t\t<bc name="p_out" node_set="OutletP" type="prescribed fluid dilatation"><value lc="2">{-P_OUT / a.k:.15g}</value><relative>0</relative></bc>')
        if a.pins:
            plo, phi = [float(v) for v in a.pins.split(",")]
            w(f'\t\t<bc name="pin_lo" node_set="PinLo" type="prescribed fluid dilatation"><value lc="2">{-plo / a.k:.15g}</value><relative>0</relative></bc>')
            w(f'\t\t<bc name="pin_hi" node_set="PinHi" type="prescribed fluid dilatation"><value lc="2">{-phi / a.k:.15g}</value><relative>0</relative></bc>')
        w('\t</Boundary>')
        if a.mode == "fsi":
            w('\t<Loads>')
            w('\t\t<surface_load name="FSItraction" surface="Interface" type="fluid-FSI traction"><shell_bottom>0</shell_bottom></surface_load>')
            w(f'\t\t<surface_load name="pext" surface="Top" type="pressure"><pressure lc="1">{P_EXT:g}</pressure><symmetric_stiffness>0</symmetric_stiffness><linear>0</linear><shell_bottom>0</shell_bottom></surface_load>')
            w('\t</Loads>')
        w('\t<Initial>')
        for k_, (v_, ids_) in sets_v.items(): w(f'\t\t<ic type="initial fluid velocity" node_set="{k_}"><value>{v_:.17g},0,0</value></ic>')
        for k_, (v_, ids_) in sets_ef.items(): w(f'\t\t<ic type="initial fluid dilatation" node_set="{k_}"><value>{v_:.17g}</value></ic>')
        w('\t</Initial>')
        tend = a.t_end
        pts = pext_table(tend)
        w('\t<LoadData>')
        w('\t\t<load_controller id="1" name="pext" type="loadcurve"><interpolate>LINEAR</interpolate><extend>CONSTANT</extend><points>' + "".join(f"<pt>{t:.10g},{g:.15g}</pt>" for t, g in pts) + '</points></load_controller>')
        w('\t\t<load_controller id="2" name="const" type="loadcurve"><interpolate>LINEAR</interpolate><extend>CONSTANT</extend><points><pt>0,1</pt><pt>1,1</pt></points></load_controller>')
        w('\t</LoadData>')
        w('\t<Output><plotfile type="febio"><var type="fluid pressure"/>' + ('<var type="displacement"/><var type="relative fluid velocity"/>' if a.mode == "fsi" else '<var type="fluid velocity"/>') + '<var type="fluid dilatation"/></plotfile>\n\t\t<logfile>')
        for k_, (dat, ids_) in mons.items(): w(f'\t\t\t<node_data data="{dat}" file="{k_}.txt" delim=" ">{lst(ids_)}</node_data>')
        w('\t\t</logfile></Output>\n</febio_spec>')
        text = "\n".join(o) + "\n"
        meta.update(n_nodes=len(nodes), n_solid_elem=len(pl_els), n_top_faces=len(top))
        os.makedirs(a.out, exist_ok=True)
        json.dump({"xs": xs, "nr": nr, "n2": n2, "i0w": i0w, "i1w": i1w, "nodes2": nodes2, "quads": quads, "iface": iface, "inlet": inlet, "outlet": outlet, "mon_rows": mon_rows, "junc": junc},
                  open(os.path.join(a.out, "mesh2d.json"), "w"))
    os.makedirs(a.out, exist_ok=True)
    open(os.path.join(a.out, tag + ".feb"), "w", encoding="latin-1").write(text)
    meta["feb_sha256"] = hashlib.sha256(text.encode("latin-1")).hexdigest()
    json.dump(meta, open(os.path.join(a.out, "model_meta.json"), "w"), indent=1, default=float)
    print(tag, a.mode, meta.get("n_nodes"), "nodes", meta.get("n_elem", meta.get("n_fluid_elem")), "elems", meta["feb_sha256"][:12])


if __name__ == "__main__":
    main()
