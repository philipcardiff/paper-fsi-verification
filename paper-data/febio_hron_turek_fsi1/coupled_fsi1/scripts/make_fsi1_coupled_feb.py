#!/usr/bin/env python3
"""Hron-Turek FSI1 COUPLED model for FEBio v4.13 (fluid-FSI module): the stage-C CFD1 fluid mesh (24 blocks, same generator functions) with the flag as a deformable
St.Venant-Kirchhoff solid (rho 1000, nu 0.4, E 1.4e6) meshed conformally with the fluid mesh on its wetted surface (top, bottom, tip), clamped on its root face, rigid cylinder
(u fixed), `fluid-FSI traction` on the flag surface, no-slip (w = 0) on cylinder, flag, top and bottom walls, the stage-C inlet and outlet formulation, pseudo-transient to the steady
state. Fluid mesh motion by the massless regularising solid of the mixture (E_mesh). Same physical problem as the Featflow FSI1 reference (u_x(A) 2.2705e-5 m, u_y(A) 8.2088e-4 m,
A = (0.6, 0.2)); the 2-D problem is modelled as a one-element-thick slab (u_z = 0, w_z = 0).
usage: make_fsi1_coupled_feb.py --out DIR --tag TAG --factor F [--k 1e6] [--dt 0.5] [--t-end 40] [--ramp 2] [--pin-root-values a,b] [--e-mesh 1e-3] [--rhoi 0] [--no-fsi-traction]"""
import argparse, os, json, hashlib, math, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_hronturek_cfd_feb as C

TZ, HCH, LCH, MU, RHO, UBAR, XROOT = C.TZ, C.HCH, C.LCH, C.MU, C.RHO, C.UBAR, C.XROOT
E_FLAG, NU_FLAG, RHO_FLAG, XTIP = 1.4e6, 0.4, 1000.0, 0.6


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--tag", required=True)
    ap.add_argument("--factor", type=float, default=1.0); ap.add_argument("--k", type=float, default=1e6)
    ap.add_argument("--dt", type=float, default=0.5); ap.add_argument("--t-end", type=float, default=40.0); ap.add_argument("--ramp", type=float, default=2.0)
    ap.add_argument("--rhoi", type=float, default=0.0); ap.add_argument("--e-mesh", type=float, default=1e-3)
    ap.add_argument("--pin-root-values", default="0,0"); ap.add_argument("--no-fsi-traction", action="store_true")
    ap.add_argument("--max-refs", type=int, default=15)
    ap.add_argument("--relaxed", action="store_true", help="declared fallback: disable the displacement/velocity/dilatation relative criteria (they have a zero initial norm and sit on an arithmetic floor), keep residual (rtol 1e-10) and energy (etol 1e-8)")
    a = ap.parse_args()
    nodes2, quads, qblock = C.mesh(a.factor); n2 = len(nodes2)
    cnt = {}
    for q in quads:
        for e in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])):
            cnt.setdefault(frozenset(e), []).append(e)
    assert all(len(v) <= 2 for v in cnt.values())
    bedges = [v[0] for v in cnt.values() if len(v) == 1]
    eps = 1e-7
    bnodes = {i for e in bedges for i in e}
    inlet = sorted(i for i in bnodes if abs(nodes2[i][0]) < eps)
    outlet = sorted(i for i in bnodes if abs(nodes2[i][0] - LCH) < eps)
    bottom = sorted(i for i in bnodes if abs(nodes2[i][1]) < eps)
    top = sorted(i for i in bnodes if abs(nodes2[i][1] - HCH) < eps)
    obst = sorted(bnodes - set(inlet) - set(outlet) - set(bottom) - set(top))
    in_free = [i for i in inlet if 1e-7 < nodes2[i][1] < HCH - 1e-7]
    in_corner = [i for i in inlet if i not in in_free]
    out_corner = [i for i in outlet if i in bottom or i in top]
    noslip = sorted(set(obst) | set(bottom) | set(top) | set(in_corner) | set(out_corner))
    key = lambda x, y: (round(x * 1e8), round(y * 1e8))
    _h = {}
    for i in range(n2):
        _h.setdefault((round(nodes2[i][0] * 1e5), round(nodes2[i][1] * 1e5)), []).append(i)
    class _F:
        def get(self, k):
            x, y = k[0] / 1e8, k[1] / 1e8; cx, cy = round(x * 1e5), round(y * 1e5)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for i in _h.get((cx + dx, cy + dy), []):
                        if abs(nodes2[i][0] - x) < 1e-7 and abs(nodes2[i][1] - y) < 1e-7:
                            return i
            return None
    fkey = _F()
    # ---- flag grid from the fluid wetted nodes
    flag_obst = [i for i in obst if nodes2[i][0] >= XROOT - 1e-9 and 0.19 - 1e-9 <= nodes2[i][1] <= 0.21 + 1e-9]
    topx = sorted({round(nodes2[i][0], 12) for i in flag_obst if abs(nodes2[i][1] - 0.21) < 1e-9})
    botx = sorted({round(nodes2[i][0], 12) for i in flag_obst if abs(nodes2[i][1] - 0.19) < 1e-9})
    assert len(topx) == len(botx) and all(abs(p - q) < 1e-9 for p, q in zip(topx, botx)), "top/bottom wetted node positions differ"
    tipy = sorted(nodes2[i][1] for i in flag_obst if abs(nodes2[i][0] - XTIP) < 1e-9)
    nx, ny = len(topx) - 1, len(tipy) - 1
    assert abs(topx[0] - XROOT) < 1e-8 and abs(topx[-1] - XTIP) < 1e-9 and abs(tipy[0] - 0.19) < 1e-9 and abs(tipy[-1] - 0.21) < 1e-9
    n_fluid3 = 2 * n2
    nodes = [(x, y, 0.0) for x, y in nodes2] + [(x, y, TZ) for x, y in nodes2]
    fl_id = lambda i, k: i + 1 + (n2 if k else 0)
    gid = {}
    def flag_node(ix, iy, k):
        x, y = topx[ix], tipy[iy]
        f = fkey.get(key(x, y))
        on_wet = (iy == 0 or iy == ny or ix == nx)
        if on_wet and f is not None:
            return fl_id(f, k)
        if on_wet and f is None and not (ix == 0):   # wetted but no fluid node: mesh mismatch
            raise AssertionError(("flag surface node without fluid partner", ix, iy, x, y))
        if (ix, iy, k) not in gid:
            nodes.append((x, y, TZ if k else 0.0)); gid[(ix, iy, k)] = len(nodes)
        return gid[(ix, iy, k)]
    # root column: junction nodes (iy == 0, ny) exist in the fluid mesh; interior root nodes are new
    flag_els = []
    for ix in range(nx):
        for iy in range(ny):
            c = lambda dx, dy, k: flag_node(ix + dx, iy + dy, k)
            flag_els.append([c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(0, 1, 0), c(0, 0, 1), c(1, 0, 1), c(1, 1, 1), c(0, 1, 1)])
    root_nodes = sorted({flag_node(0, iy, k) for iy in range(ny + 1) for k in (0, 1)})
    A = flag_node(nx, ny // 2, 0)
    assert ny % 2 == 0 and abs(tipy[ny // 2] - 0.2) < 1e-9, "no tip node at y = 0.2"
    # flag wetted surface quads (bottom row, top row, tip column), nodes shared with the fluid mesh
    surf = []
    for ix in range(nx):
        for iy, kind in ((0, "b"), (ny, "t")):
            p, q = flag_node(ix, iy, 0), flag_node(ix + 1, iy, 0); surf.append((p, q, flag_node(ix + 1, iy, 1), flag_node(ix, iy, 1)))
    for iy in range(ny):
        p, q = flag_node(nx, iy, 0), flag_node(nx, iy + 1, 0); surf.append((p, q, flag_node(nx, iy + 1, 1), flag_node(nx, iy, 1)))
    # fluid elements
    els = [[fl_id(q[0], 0), fl_id(q[1], 0), fl_id(q[2], 0), fl_id(q[3], 0), fl_id(q[0], 1), fl_id(q[1], 1), fl_id(q[2], 1), fl_id(q[3], 1)] for q in quads]
    both = lambda s: sorted([fl_id(i, 0) for i in s] + [fl_id(i, 1) for i in s])
    lst = lambda l: ",".join(map(str, l))
    flagnodes_in_fluid = {fl_id(i, k) for i in flag_obst for k in (0, 1)}
    cyl_only = [i for i in obst if i not in set(flag_obst)]
    fixu_nodes = sorted(set(both(inlet)) | set(both(outlet)) | set(both(bottom)) | set(both(top)) | set(both(cyl_only)) | set(root_nodes))
    root_f = [i for i in obst if abs(nodes2[i][0] - XROOT) < 1e-6 and (abs(nodes2[i][1] - 0.19) < 1e-6 or abs(nodes2[i][1] - 0.21) < 1e-6)]
    assert len(root_f) == 2
    pa, pb = [float(x) for x in a.pin_root_values.split(",")]
    pinv = {min(root_f, key=lambda i: nodes2[i][1]): pa, max(root_f, key=lambda i: nodes2[i][1]): pb}
    nsteps = int(round(a.t_end / a.dt))
    o = []; w = o.append
    w('<?xml version="1.0" encoding="ISO-8859-1"?>')
    w(f'<!-- Hron-Turek FSI1 coupled (fluid-FSI): {a.tag}. Generated by tools/make_fsi1_coupled_feb.py -->')
    w('<febio_spec version="4.0">\n\t<Module type="fluid-FSI"/>')
    w('\t<Control><analysis>DYNAMIC</analysis>')
    w(f'\t\t<time_steps>{nsteps}</time_steps><step_size>{a.dt:g}</step_size>')
    w('\t\t<plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range><plot_level>PLOT_MAJOR_ITRS</plot_level>')
    w('\t\t<output_level>OUTPUT_MAJOR_ITRS</output_level><plot_stride>20</plot_stride>')
    w('\t\t<solver type="fluid-FSI"><symmetric_stiffness>non-symmetric</symmetric_stiffness><equation_scheme>staggered</equation_scheme><equation_order>default</equation_order><optimize_bw>0</optimize_bw>')
    w(f'\t\t\t<lstol>0.9</lstol><lsmin>0.01</lsmin><lsiter>5</lsiter><max_refs>{a.max_refs}</max_refs><check_zero_diagonal>0</check_zero_diagonal>')
    w('\t\t\t<reform_each_time_step>1</reform_each_time_step><diverge_reform>1</diverge_reform><min_residual>1e-30</min_residual><max_residual>1e+20</max_residual>')
    tol = '<dtol>0</dtol><vtol>0</vtol><ftol>0</ftol><etol>1e-8</etol><rtol>1e-10</rtol>' if a.relaxed else '<dtol>1e-6</dtol><vtol>1e-6</vtol><ftol>1e-6</ftol><etol>1e-8</etol><rtol>1e-10</rtol>'
    w(f'\t\t\t{tol}<rhoi>{a.rhoi:g}</rhoi><predictor>0</predictor><min_volume_ratio>0</min_volume_ratio><order>2</order>')
    w('\t\t\t<qn_method type="full Newton"/></solver></Control>')
    w('\t<Material>')
    w(f'\t\t<material id="1" name="Flag" type="St.Venant-Kirchhoff"><density>{RHO_FLAG:g}</density><E>{E_FLAG:g}</E><v>{NU_FLAG:g}</v></material>')
    w('\t\t<material id="2" name="Fluid" type="fluid-FSI">')
    w(f'\t\t\t<solid type="neo-Hookean"><density>0</density><E>{a.e_mesh:g}</E><v>0</v></solid>')
    w(f'\t\t\t<fluid type="fluid"><density>{RHO:g}</density><k>{a.k:g}</k><viscous type="Newtonian fluid"><kappa>0</kappa><mu>{MU:g}</mu></viscous></fluid>')
    w('\t\t</material>\n\t</Material>')
    w('\t<Mesh>\n\t\t<Nodes name="AllNodes">')
    for n, (x, y, z) in enumerate(nodes, 1):
        w(f'\t\t\t<node id="{n}">{x:.15g},{y:.15g},{z:.15g}</node>')
    w('\t\t</Nodes>')
    w('\t\t<Elements type="hex8" name="Flag">')
    eid = 0
    for e in flag_els:
        eid += 1; w(f'\t\t\t<elem id="{eid}">{lst(e)}</elem>')
    w('\t\t</Elements>\n\t\t<Elements type="hex8" name="Fluid">')
    for e in els:
        eid += 1; w(f'\t\t\t<elem id="{eid}">{lst(e)}</elem>')
    w('\t\t</Elements>')
    w(f'\t\t<NodeSet name="AllNodes_z">{lst(range(1, len(nodes) + 1))}</NodeSet>')
    w(f'\t\t<NodeSet name="FluidAll">{lst(both(range(n2)))}</NodeSet>')
    w(f'\t\t<NodeSet name="FixU">{lst(fixu_nodes)}</NodeSet>')
    w(f'\t\t<NodeSet name="NoSlip">{lst(both(noslip))}</NodeSet>')
    w(f'\t\t<NodeSet name="Outlet">{lst(both(outlet))}</NodeSet>')
    for nm_, i_ in (("RootPinLow", min(pinv, key=lambda i: nodes2[i][1])), ("RootPinHigh", max(pinv, key=lambda i: nodes2[i][1]))):
        w(f'\t\t<NodeSet name="{nm_}">{lst(both([i_]))}</NodeSet>')
    w(f'\t\t<NodeSet name="PointA">{A}</NodeSet>')
    w(f'\t\t<NodeSet name="RootFace">{lst(root_nodes)}</NodeSet>')
    w('\t\t<Surface name="InletSurf">')
    inl = [e for e in bedges if e[0] in set(inlet) and e[1] in set(inlet)]
    for k, (p_, q_) in enumerate(inl, 1):
        w(f'\t\t\t<quad4 id="{k}">{fl_id(p_,0)},{fl_id(q_,0)},{fl_id(q_,1)},{fl_id(p_,1)}</quad4>')
    w('\t\t</Surface>\n\t\t<Surface name="FlagSurf">')
    for k, s in enumerate(surf, 1):
        w(f'\t\t\t<quad4 id="{k}">{lst(s)}</quad4>')
    w('\t\t</Surface>')
    w('\t</Mesh>\n\t<MeshDomains><SolidDomain name="Flag" mat="Flag"/><SolidDomain name="Fluid" mat="Fluid"/></MeshDomains>')
    w('\t<Boundary>')
    w('\t\t<bc name="uz" node_set="AllNodes_z" type="zero displacement"><x_dof>0</x_dof><y_dof>0</y_dof><z_dof>1</z_dof></bc>')
    w('\t\t<bc name="uxy" node_set="FixU" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>')
    w('\t\t<bc name="wz" node_set="FluidAll" type="zero fluid velocity"><wx_dof>0</wx_dof><wy_dof>0</wy_dof><wz_dof>1</wz_dof></bc>')
    w('\t\t<bc name="noslip" node_set="NoSlip" type="zero fluid velocity"><wx_dof>1</wx_dof><wy_dof>1</wy_dof><wz_dof>0</wz_dof></bc>')
    w('\t\t<bc name="outlet_p0" node_set="Outlet" type="zero fluid dilatation"/>')
    for nm_, i_ in (("RootPinLow", min(pinv, key=lambda i: nodes2[i][1])), ("RootPinHigh", max(pinv, key=lambda i: nodes2[i][1]))):
        w(f'\t\t<bc name="{nm_}" node_set="{nm_}" type="prescribed fluid dilatation"><value lc="1">{-pinv[i_] / a.k:.12g}</value><relative>0</relative></bc>')
    w('\t</Boundary>')
    w('\t<Loads>')
    w(f'\t\t<surface_load name="InletNormalVelocity" surface="InletSurf" type="fluid normal velocity"><velocity lc="1">{-UBAR:g}</velocity><prescribe_nodal_velocities>1</prescribe_nodal_velocities><parabolic>1</parabolic><prescribe_rim_pressure>1</prescribe_rim_pressure></surface_load>')
    if not a.no_fsi_traction:
        w('\t\t<surface_load name="FSItraction" surface="FlagSurf" type="fluid-FSI traction"><shell_bottom>0</shell_bottom></surface_load>')
    w('\t</Loads>')
    pts = []; kk = 0
    while kk * a.dt < a.ramp - 1e-9:
        t = kk * a.dt; pts.append((t, 0.5 * (1 - math.cos(math.pi * t / a.ramp)))); kk += 1
    pts.append((a.ramp, 1.0)); pts.append((a.t_end, 1.0))
    w('\t<LoadData><load_controller id="1" name="ramp" type="loadcurve"><interpolate>LINEAR</interpolate><extend>CONSTANT</extend><points>' + "".join(f"<pt>{t:.12g},{g:.15g}</pt>" for t, g in pts) + '</points></load_controller></LoadData>')
    w('\t<Output><plotfile type="febio"><var type="displacement"/><var type="fluid pressure"/><var type="fluid velocity"/></plotfile>')
    w('\t\t<logfile>')
    w(f'\t\t\t<node_data data="x;y;ux;uy" file="pointA.txt" delim=" ">{A}</node_data>')
    w(f'\t\t\t<node_data data="Rx;Ry" file="react_root.txt" delim=" ">{lst(root_nodes)}</node_data>')
    w('\t\t</logfile></Output>\n</febio_spec>')
    text = "\n".join(o) + "\n"
    os.makedirs(a.out, exist_ok=True)
    open(os.path.join(a.out, f"{a.tag}.feb"), "w", encoding="latin-1").write(text)
    meta = dict(vars(a), n_nodes=len(nodes), n_fluid_el=len(els), n_flag_el=len(flag_els), flag_nx=nx, flag_ny=ny, n_root_nodes=len(root_nodes), point_A_node=A,
                feb_sha256=hashlib.sha256(text.encode("latin-1")).hexdigest(), n_steps=nsteps, tag=a.tag, generator="tools/make_fsi1_coupled_feb.py")
    json.dump(meta, open(os.path.join(a.out, "model_meta.json"), "w"), indent=1, default=float)
    print(a.tag, f"fluid_cells={len(els)} flag_cells={len(flag_els)} ({nx}x{ny}) nodes={len(nodes)} A={A} root_nodes={len(root_nodes)} fixU={len(fixu_nodes)}")


if __name__ == "__main__":
    main()
