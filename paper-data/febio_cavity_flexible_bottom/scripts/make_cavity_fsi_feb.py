#!/usr/bin/env python3
import os, os as _os; HERE = _os.path.dirname(_os.path.abspath(__file__)); PKG = _os.path.dirname(HERE)   # package root; run directories of retained HPC outputs: env CAVITY_RUNS (default runs/cavityNext)
"""cavityFlexibleBottom coupled FEBio fluid-FSI model (monolithic ALE), conformal hex8 fluid + hex8 StVK plate.
Fluid mesh = the validated 15-block cavity (tools/make_cavity_fluid_feb.py); plate 4 x 0.4 m, (40f) x (4f) hex8 cells whose top-row nodes coincide with the fluid
interface nodes (shared nodes, both z layers). Plane strain (u_z = w_z = 0 everywhere). StVK plate E = 500, nu = 0.4, rho = 1000; fluid rho = 1, mu = 0.01, weakly compressible K;
fluid-FSI mixture material with a massless neo-Hookean mesh regulariser (E_mesh, nu = 0). No-slip: w = 0 on walls and interface nodes. Inlet: built-in parabolic fluid normal velocity (+ rim pressure).
Clamped plate ends (u = 0, x = 2 and 6 over the full thickness). Fluid mesh displacement fixed on walls/inlet/outlet. Corner pin (optional): prescribed dilatation at the 4 floor-end corner nodes.
usage: make_cavity_fsi_feb.py --out DIR --tag TAG [--mesh-factor 1] [--k 1e6] [--dt 2.5] [--t-end 1250] [--ramp 25] [--rhoi 0] [--e-mesh 1e-3] [--pin-corners 1.15] [--log-all]
"""
import argparse, os, json, hashlib, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_cavity_fluid_feb as G

ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--tag", required=True)
ap.add_argument("--mesh-factor", type=int, default=1); ap.add_argument("--k", type=float, default=1e6); ap.add_argument("--dt", type=float, default=2.5)
ap.add_argument("--t-end", type=float, default=1250.0); ap.add_argument("--ramp", type=float, default=25.0); ap.add_argument("--rhoi", type=float, default=0.0)
ap.add_argument("--e-mesh", type=float, default=1e-3); ap.add_argument("--mu", type=float, default=0.01); ap.add_argument("--pin-corners", type=float, default=None)
ap.add_argument("--dtol", type=float, default=1e-3); ap.add_argument("--vtol", type=float, default=1e-6); ap.add_argument("--ftol", type=float, default=1e-6)
ap.add_argument("--etol", type=float, default=1e-4); ap.add_argument("--rtol", type=float, default=1e-6); ap.add_argument("--max-refs", type=int, default=25)
ap.add_argument("--plot-stride", type=int, default=40); ap.add_argument("--log-all", action="store_true", help="log all fluid z=0 nodes (large)")
a = ap.parse_args()
f = a.mesh_factor; TZ = G.TZ; E_S, NU_S, RHO_S = 500.0, 0.4, 1000.0
nodes2, quads = G.mesh(f); n2 = len(nodes2); eps = 1e-6
for q in quads:
    (x0, y0), (x1, y1), (x2, y2), (x3, y3) = [nodes2[i] for i in q]; assert 0.5 * ((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0)) > 0
cnt = {}
for q in quads:
    for e in ((q[0], q[1]), (q[1], q[2]), (q[2], q[3]), (q[3], q[0])): cnt.setdefault(frozenset(e), []).append(e)
bedges = [v[0] for v in cnt.values() if len(v) == 1]; bnodes = {i for e in bedges for i in e}
inlet = sorted(i for i in bnodes if abs(nodes2[i][0]) < eps); outlet = sorted(i for i in bnodes if abs(nodes2[i][0] - 14) < eps)
iface = sorted((i for i in bnodes if abs(nodes2[i][1]) < eps), key=lambda i: nodes2[i][0])
walls = sorted(bnodes - set(inlet) - set(outlet) - set(iface)); inlet_corner = [i for i in inlet if abs(nodes2[i][1] - 1) < eps or abs(nodes2[i][1] - 2) < eps]
outlet_corner = [i for i in outlet if abs(nodes2[i][1] - 1) < eps or abs(nodes2[i][1] - 2) < eps]
N1 = lambda i: i + 1; lst = lambda l: ",".join(map(str, l))
# --- nodes: fluid z=0 layer (0..n2-1), fluid z=TZ layer (n2..2n2-1), then plate nodes (sharing the top row with the fluid interface nodes)
nodes = [(x, y, 0.0) for x, y in nodes2] + [(x, y, TZ) for x, y in nodes2]
fl_els = [[N1(q[0]), N1(q[1]), N1(q[2]), N1(q[3]), N1(q[0] + n2), N1(q[1] + n2), N1(q[2] + n2), N1(q[3] + n2)] for q in quads]
nxp, nyp = 40 * f, 4 * f; assert len(iface) == nxp + 1, (len(iface), nxp)
xs_if = [nodes2[i][0] for i in iface]; assert all(abs(xs_if[k] - (2.0 + 4.0 * k / nxp)) < 1e-6 for k in range(nxp + 1)), "interface nodes must coincide with plate columns"
pidx = {}
def pn(i, j, k):
    if j == nyp:   # top row: shared with the fluid interface nodes
        return N1(iface[i] + (n2 if k else 0))
    if (i, j, k) not in pidx:
        nodes.append((2.0 + 4.0 * i / nxp, -0.4 + 0.4 * j / nyp, TZ * k)); pidx[(i, j, k)] = len(nodes)
    return pidx[(i, j, k)]
pl_els, top = [], []
for ex in range(nxp):
    for ey in range(nyp):
        c = lambda aa, bb, kk: pn(ex + aa, ey + bb, kk)
        pl_els.append([c(0, 0, 0), c(1, 0, 0), c(1, 1, 0), c(0, 1, 0), c(0, 0, 1), c(1, 0, 1), c(1, 1, 1), c(0, 1, 1)])
        if ey == nyp - 1: top.append([c(0, 1, 0), c(0, 1, 1), c(1, 1, 1), c(1, 1, 0)])   # normal +y (out of the solid, into the fluid)
plate_nodes = sorted(set(n for e in pl_els for n in e)); fluid_nodes = list(range(1, 2 * n2 + 1))
both = lambda s: sorted([N1(i) for i in s] + [N1(i + n2) for i in s])
left = sorted(n for n in plate_nodes if abs(nodes[n - 1][0] - 2.0) < eps); right = sorted(n for n in plate_nodes if abs(nodes[n - 1][0] - 6.0) < eps)
Aid = [n for n in plate_nodes if abs(nodes[n - 1][0] - 4.0) < eps and abs(nodes[n - 1][1]) < eps and abs(nodes[n - 1][2]) < eps]; assert len(Aid) == 1
bottomA = [n for n in plate_nodes if abs(nodes[n - 1][0] - 4.0) < eps and abs(nodes[n - 1][1] + 0.4) < eps and abs(nodes[n - 1][2]) < eps]
# fluid interface surface = the plate top faces (conformal): the fluid-FSI traction surface
nsteps = int(round(a.t_end / a.dt)); o = []; w = o.append
w('<?xml version="1.0" encoding="ISO-8859-1"?>'); w(f'<!-- cavityFlexibleBottom coupled FSI: {a.tag}. Generated by tools/make_cavity_fsi_feb.py -->'); w('<febio_spec version="4.0">\n\t<Module type="fluid-FSI"/>')
w('\t<Control><analysis>DYNAMIC</analysis>'); w(f'\t\t<time_steps>{nsteps}</time_steps><step_size>{a.dt:g}</step_size>')
w(f'\t\t<plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range><plot_level>PLOT_MAJOR_ITRS</plot_level><output_level>OUTPUT_MAJOR_ITRS</output_level><plot_stride>{a.plot_stride}</plot_stride>')
w('\t\t<solver type="fluid-FSI"><symmetric_stiffness>non-symmetric</symmetric_stiffness><equation_scheme>staggered</equation_scheme><equation_order>default</equation_order><optimize_bw>0</optimize_bw>')
w(f'\t\t\t<lstol>0.9</lstol><lsmin>0.01</lsmin><lsiter>5</lsiter><max_refs>{a.max_refs}</max_refs><check_zero_diagonal>0</check_zero_diagonal><reform_each_time_step>1</reform_each_time_step><diverge_reform>1</diverge_reform>')
w(f'\t\t\t<min_residual>1e-30</min_residual><max_residual>1e+20</max_residual><dtol>{a.dtol:g}</dtol><vtol>{a.vtol:g}</vtol><ftol>{a.ftol:g}</ftol><etol>{a.etol:g}</etol><rtol>{a.rtol:g}</rtol>')
w(f'\t\t\t<rhoi>{a.rhoi:g}</rhoi><predictor>0</predictor><min_volume_ratio>0</min_volume_ratio><order>2</order><qn_method type="full Newton"/></solver></Control>')
w('\t<Material>'); w(f'\t\t<material id="1" name="Plate" type="St.Venant-Kirchhoff"><density>{RHO_S:g}</density><E>{E_S:g}</E><v>{NU_S:g}</v></material>')
w('\t\t<material id="2" name="Fluid" type="fluid-FSI">'); w(f'\t\t\t<solid type="neo-Hookean"><density>0</density><E>{a.e_mesh:g}</E><v>0</v></solid>')
w(f'\t\t\t<fluid type="fluid"><density>1</density><k>{a.k:g}</k><viscous type="Newtonian fluid"><kappa>0</kappa><mu>{a.mu:g}</mu></viscous></fluid></material>\n\t</Material>')
w('\t<Mesh>\n\t\t<Nodes name="AllNodes">')
for n, (x, y, z) in enumerate(nodes, 1): w(f'\t\t\t<node id="{n}">{x:.15g},{y:.15g},{z:.15g}</node>')
w('\t\t</Nodes>\n\t\t<Elements type="hex8" name="Plate">')
for i, e in enumerate(pl_els, 1): w(f'\t\t\t<elem id="{i}">{lst(e)}</elem>')
w('\t\t</Elements>\n\t\t<Elements type="hex8" name="Fluid">')
for i, e in enumerate(fl_els, len(pl_els) + 1): w(f'\t\t\t<elem id="{i}">{lst(e)}</elem>')
w('\t\t</Elements>')
ifn = both(iface); wall_nodes = both(walls) ; noslip = sorted(set(both(walls)) | set(ifn) | set(both(inlet_corner)) | set(both(outlet_corner)))
w(f'\t\t<NodeSet name="AllNodes_z">{lst(range(1, len(nodes) + 1))}</NodeSet>')
w(f'\t\t<NodeSet name="FluidAll">{lst(fluid_nodes)}</NodeSet>'); w(f'\t\t<NodeSet name="NoSlipW">{lst(noslip)}</NodeSet>')
fixmesh = sorted(set(both(walls)) | set(both(inlet)) | set(both(outlet)))   # fluid mesh displacement fixed on walls, inlet, outlet (interface nodes belong to the plate)
w(f'\t\t<NodeSet name="FixMesh">{lst(fixmesh)}</NodeSet>'); w(f'\t\t<NodeSet name="Outlet">{lst(both(outlet))}</NodeSet>')
w(f'\t\t<NodeSet name="ClampL">{lst(left)}</NodeSet><NodeSet name="ClampR">{lst(right)}</NodeSet><NodeSet name="PointA">{Aid[0]}</NodeSet><NodeSet name="PointA_bottom">{lst(bottomA)}</NodeSet>')
corners = [i for i in iface if abs(nodes2[i][0] - 2.0) < eps or abs(nodes2[i][0] - 6.0) < eps]
if a.pin_corners is not None: w(f'\t\t<NodeSet name="CornerPin">{lst(both(corners))}</NodeSet>')
w(f'\t\t<NodeSet name="InterfaceZ0">{lst([N1(i) for i in iface])}</NodeSet>')
w('\t\t<Surface name="InletSurf">')
inl = [e for e in bedges if abs(nodes2[e[0]][0]) < eps and abs(nodes2[e[1]][0]) < eps]
for k, (p_, q_) in enumerate(inl, 1): w(f'\t\t\t<quad4 id="{k}">{N1(p_)},{N1(q_)},{N1(q_ + n2)},{N1(p_ + n2)}</quad4>')
w('\t\t</Surface>\n\t\t<Surface name="Interface">')
for k, q in enumerate(top, 1): w(f'\t\t\t<quad4 id="{k}">{lst(q)}</quad4>')
w('\t\t</Surface>')
mon = {"Mon_inlet": sorted([i for i in inlet], key=lambda i: nodes2[i][1]), "Mon_outlet": sorted(outlet, key=lambda i: nodes2[i][1])}
for k, v in mon.items(): w(f'\t\t<NodeSet name="{k}">{lst([N1(i) for i in v])}</NodeSet>')
w('\t</Mesh>\n\t<MeshDomains><SolidDomain name="Plate" mat="Plate"/><SolidDomain name="Fluid" mat="Fluid"/></MeshDomains>')
w('\t<Boundary>')
w('\t\t<bc name="uz" node_set="AllNodes_z" type="zero displacement"><x_dof>0</x_dof><y_dof>0</y_dof><z_dof>1</z_dof></bc>')
w('\t\t<bc name="clampL" node_set="ClampL" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>'); w('\t\t<bc name="clampR" node_set="ClampR" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>')
w('\t\t<bc name="fixmesh" node_set="FixMesh" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>')
w('\t\t<bc name="wz" node_set="FluidAll" type="zero fluid velocity"><wx_dof>0</wx_dof><wy_dof>0</wy_dof><wz_dof>1</wz_dof></bc>')
w('\t\t<bc name="noslip" node_set="NoSlipW" type="zero fluid velocity"><wx_dof>1</wx_dof><wy_dof>1</wy_dof><wz_dof>0</wz_dof></bc>')
w('\t\t<bc name="outlet_p0" node_set="Outlet" type="zero fluid dilatation"/>')
if a.pin_corners is not None: w(f'\t\t<bc name="corner_pin" node_set="CornerPin" type="prescribed fluid dilatation"><value lc="1">{-a.pin_corners / a.k:.12g}</value><relative>0</relative></bc>')
w('\t</Boundary>\n\t<Loads>')
w('\t\t<surface_load name="FSItraction" surface="Interface" type="fluid-FSI traction"><shell_bottom>0</shell_bottom></surface_load>')
w('\t\t<surface_load name="InletNormalVelocity" surface="InletSurf" type="fluid normal velocity"><velocity lc="1">-1</velocity><prescribe_nodal_velocities>1</prescribe_nodal_velocities><parabolic>1</parabolic><prescribe_rim_pressure>1</prescribe_rim_pressure></surface_load>')
w('\t</Loads>')
w(f'\t<LoadData><load_controller id="1" name="ramp" type="loadcurve"><interpolate>LINEAR</interpolate><extend>CONSTANT</extend><points><pt>0,0</pt><pt>{a.ramp:g},1</pt></points></load_controller></LoadData>')
w('\t<Output><plotfile type="febio"><var type="displacement"/><var type="fluid pressure"/><var type="relative fluid velocity"/><var type="fluid dilatation"/><var type="reaction forces"/></plotfile>\n\t\t<logfile>')
w(f'\t\t\t<node_data data="x;y;ux;uy" file="pointA.txt" delim=" ">{Aid[0]}</node_data>')
w(f'\t\t\t<node_data data="ux;uy" file="pointA_bottom.txt" delim=" ">{lst(bottomA)}</node_data>')
w(f'\t\t\t<node_data data="ux;uy;wx;wy;ef" file="Mon_iface.txt" delim=" ">{lst([N1(i) for i in iface])}</node_data>')
if a.log_all: w(f'\t\t\t<node_data data="ux;uy;wx;wy;ef" file="Mon_all.txt" delim=" ">{lst([N1(i) for i in range(n2)])}</node_data>')
for k in mon: w(f'\t\t\t<node_data data="wx;wy;ef" file="{k}.txt" delim=" ">{lst([N1(i) for i in mon[k]])}</node_data>')
w('\t\t</logfile></Output>\n</febio_spec>')
text = "\n".join(o) + "\n"; os.makedirs(a.out, exist_ok=True); open(os.path.join(a.out, a.tag + ".feb"), "w", encoding="latin-1").write(text)
meta = dict(vars(a), n_nodes=len(nodes), n_fluid_elem=len(fl_els), n_solid_elem=len(pl_els), n_interface_facets=len(top), point_A_node=Aid[0], thickness_z=TZ, feb_sha256=hashlib.sha256(text.encode("latin-1")).hexdigest(), purpose="coupled FSI")
json.dump(meta, open(os.path.join(a.out, "model_meta.json"), "w"), indent=1)
print(a.tag, len(nodes), "nodes", len(fl_els), "fluid hex8", len(pl_els), "plate hex8", len(top), "interface facets; point A node", Aid[0])
