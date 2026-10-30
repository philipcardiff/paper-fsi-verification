#!/usr/bin/env python3
"""Generate a FEBio 4.x (fluid-FSI) input for the ringAddedMass benchmark.

Problem (solids4foam tutorial `ringAddedMass`, development 1576edcc):
  plane-strain elastic ring Ri=0.95, Ro=1.05 m (mean radius a=1, h=0.1),
  E=1e6 Pa, nu=0.3, rho_s=1000; fluid in the annulus Ro<r<b=1.5 (rigid wall),
  rho_f in {14,140,1400}, nu_f=1e-7; n=2 (ovalling) mode, V0=1.5e-3 m/s,
  amplitude ~5e-4*a; quarter model with symmetry planes; 100 steps/wet period.

FEBio representation (documented departures from the OpenFOAM tutorial are
flagged DEPARTURE):
  * 3D slab one element thick in z (thickness 0.1 m); u_z = 0 and fluid w_z = 0
    on all nodes (plane strain).
  * ring: `isotropic elastic` (St Venant-Kirchhoff; identical linearisation),
    density 1000.
  * fluid: `fluid-FSI` mixture; fluid is WEAKLY COMPRESSIBLE (p = -K e^f).
    DEPARTURE from the incompressible reference: K is a swept parameter.
    Mesh-regularising solid: neo-Hookean, nu=0, E = E_mesh, density 0.
  * conformal interface (shared nodes) + `fluid-FSI traction` surface load.
  * no-slip: w = 0 on the interface and on the rigid wall (the reference has no-slip walls).
  * symmetry planes x=0 and y=0: mesh normal displacement fixed, fluid normal
    relative velocity fixed (frictionless tangentially).
  * IC (DEPARTURE from the impulsive start of the tutorial): ring velocity is the
    thin-ring n=2 mode velocity of the tutorial; the fluid starts with the exact
    potential-flow velocity field that matches the ring's radial velocity, so no
    acoustic pulse is launched. `--ic impulsive` reproduces the tutorial
    (fluid at rest) as a diagnostic only.
  * time integration: fluid-FSI solver generalized-alpha; `rhoi` is a parameter.
"""
from __future__ import annotations
import argparse, hashlib, json, math, os, sys

RI, RO, B, H = 0.95, 1.05, 1.5, 0.1
E_RING, NU_RING, RHO_S = 1.0e6, 0.3, 1000.0
NU_F = 1.0e-7
V0 = 1.5e-3
OMEGA_REF = {"dry": 2.55504257, 14.0: 2.43561524, 140.0: 1.80451344, 1400.0: 0.76865460}


def mode_expressions(v0: float, b: float, ro: float):
    """Cartesian expressions in X,Y for ring velocity and potential-flow fluid velocity."""
    q = "(X^2+Y^2)"
    r = "sqrt(X^2+Y^2)"
    c2 = f"((X^2-Y^2)/{q})"      # cos 2theta
    s2 = f"((2*X*Y)/{q})"        # sin 2theta
    # ring: Vr = V0 cos2t ; Vt = V0 sin2t (-1/2 + 3 (r-a)/(2a)), a = 1
    vr = f"({v0!r}*{c2})"
    vt = f"({v0!r}*{s2}*(-0.5+1.5*({r}-1.0)))"
    ring_x = f"{vr}*X/{r}-{vt}*Y/{r}"
    ring_y = f"{vr}*Y/{r}+{vt}*X/{r}"
    # fluid potential phi=(A r^2 + Bc r^-2) cos2t, dphi/dr=0 at r=b, dphi/dr=V0 at r=Ro
    A = v0 / (2.0 * (ro - b ** 4 * ro ** -3))
    Bc = A * b ** 4
    # r^3 written as q*r: the FEBio expression parser crashed (SIGSEGV, v4.13) on `sqrt(...)^3`
    r3 = f"({q}*{r})"
    # coefficients are negative: parenthesise them (the parser mishandles `--x` and `+-x`)
    Ap, Bp = f"({A!r})", f"({Bc!r})"
    fr = f"(2*({Ap}*{r}-{Bp}/{r3})*{c2})"
    ft = f"(-2*({Ap}*{r}+{Bp}/{r3})*{s2})"
    fl_x = f"{fr}*X/{r}-{ft}*Y/{r}"
    fl_y = f"{fr}*Y/{r}+{ft}*X/{r}"
    return (ring_x, ring_y), (fl_x, fl_y), (A, Bc)


def fluid_velocity(X, Y, mode="potential", wall_noslip=True):
    """Absolute fluid velocity (vx, vy) imposed as the initial relative velocity w at a free node (mesh velocity is 0 there).

    mode 'potential'  : exact irrotational n=2 flow matching the ring's radial velocity (v1 input).
    mode 'compatible' : potential flow + stationary solenoidal correction (stream function g(r) sin2t, cubic Hermite,
                        g=0 at r=Ro and r=b so no normal velocity is added) that makes the tangential velocity equal the
                        ring's at r=Ro and zero at r=b (when the wall is no-slip). Exactly divergence free.
    """
    (_, _), (_, _), (A, Bc) = mode_expressions(V0, B, RO)
    r = math.hypot(X, Y); th = math.atan2(Y, X)
    vr = 2 * (A * r - Bc / r ** 3) * math.cos(2 * th)
    vt = -2 * (A * r + Bc / r ** 3) * math.sin(2 * th)
    if mode == "compatible":
        pot_t = lambda rr: -2 * (A * rr + Bc / rr ** 3)               # coefficient of sin2t
        ring_t = V0 * (-0.5 + 1.5 * (RO - 1.0))                       # ring surface tangential velocity coefficient
        g1 = pot_t(RO) - ring_t                                       # g'(Ro)
        g2 = pot_t(B) if wall_noslip else 0.0                         # g'(b)
        L = B - RO; s = (r - RO) / L
        g = L * (g1 * (s - 2 * s ** 2 + s ** 3) + g2 * (s ** 3 - s ** 2))
        gp = g1 * (1 - 4 * s + 3 * s ** 2) + g2 * (3 * s ** 2 - 2 * s)
        vr += 2 * g / r * math.cos(2 * th)
        vt += -gp * math.sin(2 * th)
    return (vr * math.cos(th) - vt * math.sin(th), vr * math.sin(th) + vt * math.cos(th))


def evaluate(expr: str, X: float, Y: float) -> float:
    return eval(expr.replace("^", "**"), {"sqrt": math.sqrt}, {"X": X, "Y": Y})


def build(args):
    f = args.mesh_factor
    nrs, nt, nrf = 4 * f, 48 * f, 12 * f
    wet = args.rho_f is not None
    nodes, ring_el, fl_el = [], [], []
    idx = {}

    def add(key, x, y, z):
        if key in idx:
            return idx[key]
        nodes.append((x, y, z)); idx[key] = len(nodes); return idx[key]

    dth = (math.pi / 2) / nt
    for k in (0, 1):
        z = k * H
        for j in range(nt + 1):
            th = j * dth
            for i in range(nrs + 1):
                r = RI + (RO - RI) * i / nrs
                add(("s", i, j, k), r * math.cos(th), r * math.sin(th), z)
    if wet:
        for k in (0, 1):
            z = k * H
            for j in range(nt + 1):
                th = j * dth
                for i in range(nrf + 1):
                    r = RO + (B - RO) * i / nrf
                    key = ("s", nrs, j, k) if i == 0 else ("f", i, j, k)
                    add(key, r * math.cos(th), r * math.sin(th), z)

    def hex8(kind, i, j):
        g = lambda ii, jj, kk: idx[(kind, ii, jj, kk)]
        if kind == "f" and i == 0:   # first fluid layer uses shared interface nodes
            g0 = lambda jj, kk: idx[("s", nrs, jj, kk)]
        else:
            g0 = lambda jj, kk: g(i, jj, kk)
        n1, n2, n3, n4 = g0(j, 0), g(i + 1, j, 0), g(i + 1, j + 1, 0), g0(j + 1, 0)
        n5, n6, n7, n8 = g0(j, 1), g(i + 1, j, 1), g(i + 1, j + 1, 1), g0(j + 1, 1)
        return (n1, n2, n3, n4, n5, n6, n7, n8)

    for j in range(nt):
        for i in range(nrs):
            ring_el.append(hex8("s", i, j))
    if wet:
        for j in range(nt):
            for i in range(nrf):
                fl_el.append(hex8("f", i, j))

    ring_nodes = sorted({n for e in ring_el for n in e})
    fl_nodes = sorted({n for e in fl_el for n in e}) if wet else []
    iface_nodes = sorted(idx[("s", nrs, j, k)] for j in range(nt + 1) for k in (0, 1)) if wet else []
    wall_nodes = sorted(idx[("f", nrf, j, k)] for j in range(nt + 1) for k in (0, 1)) if wet else []
    symx = sorted(n for n, (x, y, z) in enumerate(nodes, 1) if abs(x) < 1e-12)
    symy = sorted(n for n, (x, y, z) in enumerate(nodes, 1) if abs(y) < 1e-12)
    symx_f = sorted(set(symx) & set(fl_nodes)); symy_f = sorted(set(symy) & set(fl_nodes))
    iface_faces = []
    if wet:
        for j in range(nt):   # interface = outer ring faces (i = nrs), normal sense handled by FEBio
            iface_faces.append((idx[("s", nrs, j, 0)], idx[("s", nrs, j + 1, 0)],
                                idx[("s", nrs, j + 1, 1)], idx[("s", nrs, j, 1)]))
    return dict(nodes=nodes, ring_el=ring_el, fl_el=fl_el, ring_nodes=ring_nodes, fl_nodes=fl_nodes,
                iface=iface_nodes, wall=wall_nodes, symx=symx, symy=symy, symx_f=symx_f, symy_f=symy_f,
                faces=iface_faces, idx=idx, dims=(nrs, nt, nrf))


def lst(ids):
    return ",".join(str(i) for i in ids)


def make(args):
    if not hasattr(args, "wall"): args.wall = "noslip"
    if not hasattr(args, "interface"): args.interface = "noslip"
    if not hasattr(args, "ic_field"): args.ic_field = "potential"
    m = build(args)
    wet = args.rho_f is not None
    nrs, nt, nrf = m["dims"]
    key = args.rho_f if wet else "dry"
    omega = OMEGA_REF[key]
    period = 2 * math.pi / omega
    dt = period / args.steps_per_period
    nsteps = int(round(args.periods * args.steps_per_period))
    (rx, ry), (fx, fy), (Acoef, Bcoef) = mode_expressions(V0, B, RO)

    allnodes = list(range(1, len(m["nodes"]) + 1))
    fixux = sorted(set(m["symx"]) | set(m["wall"]))
    fixuy = sorted(set(m["symy"]) | set(m["wall"]))
    wall_w = set(m["wall"]) if args.wall == "noslip" else set()     # wall nodes with w = 0 (no-slip); slip: natural BC only
    # interface nodes: no-slip fixes w=0; 'slip' leaves w free so only the natural condition w.n = 0 (weak impermeability,
    # Theory Manual eq. 3.6.3) and t_tau = 0 apply. The initial w at interface nodes stays 0 in both cases, and
    # fluid_free (the IC set) is deliberately unchanged, so the two inputs differ only in the FixWX/FixWY node sets.
    iface_w = set(m["iface"]) if args.interface == "noslip" else set()
    fixwx = sorted(set(m["symx_f"]) | iface_w | wall_w)
    fixwy = sorted(set(m["symy_f"]) | iface_w | wall_w)
    fluid_free = sorted(set(m["fl_nodes"]) - set(m["iface"]) - wall_w)
    probes = []
    ix = m["idx"]
    probes = [ix[("s", nrs, 0, 0)], ix[("s", nrs, nt, 0)], ix[("s", 0, 0, 0)], ix[("s", 0, nt, 0)]]

    o = []
    w = o.append
    w('<?xml version="1.0" encoding="ISO-8859-1"?>')
    w(f'<!-- ringAddedMass FEBio model: {args.tag}. Generated by tools/make_ring_feb.py. See PHASE_B_RING_PILOT.md -->')
    w('<febio_spec version="4.0">')
    w('\t<Module type="fluid-FSI"/>')
    w('\t<Control>')
    w('\t\t<analysis>DYNAMIC</analysis>')
    w(f'\t\t<time_steps>{nsteps}</time_steps>')
    w(f'\t\t<step_size>{dt:.12g}</step_size>')
    w('\t\t<plot_zero_state>0</plot_zero_state>')
    w('\t\t<plot_range>0,-1</plot_range>')
    w('\t\t<plot_level>PLOT_MAJOR_ITRS</plot_level>')
    w('\t\t<output_level>OUTPUT_MAJOR_ITRS</output_level>')
    w(f'\t\t<plot_stride>{args.plot_stride}</plot_stride>')
    w('\t\t<solver type="fluid-FSI">')
    w('\t\t\t<symmetric_stiffness>non-symmetric</symmetric_stiffness>')
    w('\t\t\t<equation_scheme>staggered</equation_scheme>')
    w('\t\t\t<equation_order>default</equation_order>')
    w('\t\t\t<optimize_bw>0</optimize_bw>')
    w('\t\t\t<lstol>0.9</lstol><lsmin>0.01</lsmin><lsiter>5</lsiter>')
    w(f'\t\t\t<max_refs>{args.max_refs}</max_refs>')
    w('\t\t\t<check_zero_diagonal>0</check_zero_diagonal>')
    w('\t\t\t<reform_each_time_step>1</reform_each_time_step>')
    w('\t\t\t<diverge_reform>1</diverge_reform>')
    w('\t\t\t<min_residual>1e-30</min_residual><max_residual>1e+20</max_residual>')
    w(f'\t\t\t<dtol>{args.dtol:g}</dtol><vtol>{args.vtol:g}</vtol><ftol>{args.ftol:g}</ftol>')
    w(f'\t\t\t<etol>{args.etol:g}</etol><rtol>{args.rtol:g}</rtol>')
    w(f'\t\t\t<rhoi>{args.rhoi:g}</rhoi>')
    w('\t\t\t<predictor>0</predictor><min_volume_ratio>0</min_volume_ratio><order>2</order>')
    w('\t\t\t<qn_method type="full Newton"/>')
    w('\t\t</solver>')
    w('\t</Control>')
    w('\t<Material>')
    w(f'\t\t<material id="1" name="Ring" type="isotropic elastic">')
    w(f'\t\t\t<density>{RHO_S:g}</density><E>{E_RING:g}</E><v>{NU_RING:g}</v>')
    w('\t\t</material>')
    if wet:
        mu = args.rho_f * (args.nu_f if args.nu_f is not None else NU_F)
        w('\t\t<material id="2" name="Fluid" type="fluid-FSI">')
        w(f'\t\t\t<solid type="neo-Hookean"><density>0</density><E>{args.e_mesh:g}</E><v>0</v></solid>')
        w(f'\t\t\t<fluid type="fluid"><density>{args.rho_f:g}</density><k>{args.k:g}</k>')
        w(f'\t\t\t\t<viscous type="Newtonian fluid"><kappa>0</kappa><mu>{mu:.12g}</mu></viscous></fluid>')
        w('\t\t</material>')
    w('\t</Material>')
    w('\t<Mesh>')
    w('\t\t<Nodes name="AllNodes">')
    for n, (x, y, z) in enumerate(m["nodes"], 1):
        w(f'\t\t\t<node id="{n}">{x:.15g},{y:.15g},{z:.15g}</node>')
    w('\t\t</Nodes>')
    eid = 0
    w('\t\t<Elements type="hex8" name="Ring">')
    for e in m["ring_el"]:
        eid += 1; w(f'\t\t\t<elem id="{eid}">{lst(e)}</elem>')
    w('\t\t</Elements>')
    if wet:
        w('\t\t<Elements type="hex8" name="Fluid">')
        for e in m["fl_el"]:
            eid += 1; w(f'\t\t\t<elem id="{eid}">{lst(e)}</elem>')
        w('\t\t</Elements>')
    sets = {"AllNodes_z": allnodes, "FixUX": fixux, "FixUY": fixuy, "RingNodes": m["ring_nodes"], "Probes": probes}
    if wet:
        sets.update({"FluidAll": m["fl_nodes"], "FixWX": fixwx, "FixWY": fixwy, "FluidFree": fluid_free})
    for name, ids in sets.items():
        w(f'\t\t<NodeSet name="{name}">{lst(ids)}</NodeSet>')
    pernode = wet and args.ic == "potential" and args.ic_mode == "pernode"
    if args.ic_field == "compatible" and args.ic_mode == "math":
        raise SystemExit("--ic-field compatible requires --ic-mode pernode")
    if pernode:   # v4.13 workaround: one constant-vector IC per node (see PHASE_B root-cause note)
        for n in fluid_free:
            w(f'\t\t<NodeSet name="fv{n}">{n}</NodeSet>')
    if wet:
        w('\t\t<Surface name="Interface">')
        for i, q in enumerate(m["faces"], 1):
            w(f'\t\t\t<quad4 id="{i}">{lst(q)}</quad4>')
        w('\t\t</Surface>')
    w('\t</Mesh>')
    w('\t<MeshDomains>')
    w('\t\t<SolidDomain name="Ring" mat="Ring"/>')
    if wet:
        w('\t\t<SolidDomain name="Fluid" mat="Fluid"/>')
    w('\t</MeshDomains>')
    w('\t<Boundary>')
    w('\t\t<bc name="uz" node_set="AllNodes_z" type="zero displacement"><x_dof>0</x_dof><y_dof>0</y_dof><z_dof>1</z_dof></bc>')
    w('\t\t<bc name="ux" node_set="FixUX" type="zero displacement"><x_dof>1</x_dof><y_dof>0</y_dof><z_dof>0</z_dof></bc>')
    w('\t\t<bc name="uy" node_set="FixUY" type="zero displacement"><x_dof>0</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>')
    if wet:
        w('\t\t<bc name="wz" node_set="FluidAll" type="zero fluid velocity"><wx_dof>0</wx_dof><wy_dof>0</wy_dof><wz_dof>1</wz_dof></bc>')
        w('\t\t<bc name="wx" node_set="FixWX" type="zero fluid velocity"><wx_dof>1</wx_dof><wy_dof>0</wy_dof><wz_dof>0</wz_dof></bc>')
        w('\t\t<bc name="wy" node_set="FixWY" type="zero fluid velocity"><wx_dof>0</wx_dof><wy_dof>1</wy_dof><wz_dof>0</wz_dof></bc>')
    w('\t</Boundary>')
    if wet:
        w('\t<Loads>')
        w('\t\t<surface_load name="FSItraction" surface="Interface" type="fluid-FSI traction"><shell_bottom>0</shell_bottom></surface_load>')
        w('\t</Loads>')
    w('\t<Initial>')
    w(f'\t\t<ic type="velocity" node_set="RingNodes"><value type="math">{rx},{ry},0</value></ic>')
    if wet and args.ic == "potential" and args.ic_mode == "math":
        w(f'\t\t<ic type="initial fluid velocity" node_set="FluidFree"><value type="math">{fx},{fy},0</value></ic>')
    if pernode:
        for n in fluid_free:
            X, Y, _ = m["nodes"][n - 1]
            vx_, vy_ = fluid_velocity(X, Y, args.ic_field, args.wall == "noslip")
            w(f'\t\t<ic type="initial fluid velocity" node_set="fv{n}"><value>{vx_:.17g},{vy_:.17g},0</value></ic>')
    w('\t</Initial>')
    w('\t<Output>')
    w('\t\t<plotfile type="febio"><var type="displacement"/><var type="velocity"/>')
    if wet:
        w('\t\t\t<var type="fluid pressure"/><var type="relative fluid velocity"/><var type="fluid volume ratio"/>')
    w('\t\t</plotfile>')
    w('\t\t<logfile>')
    w(f'\t\t\t<node_data data="x;y;ux;uy;vx;vy" file="probes.txt" delim=" ">{lst(probes)}</node_data>')
    w(f'\t\t\t<node_data data="ux;uy;vx;vy" file="ring_nodes.txt" delim=" ">{lst(m["ring_nodes"])}</node_data>')
    if wet:
        w(f'\t\t\t<node_data data="ux;uy;vx;vy;nfvx;nfvy" file="fluid_nodes.txt" delim=" ">{lst(m["fl_nodes"])}</node_data>')
        w('\t\t\t<element_data data="fp;fJ" elem_set="Fluid" file="fluid_elems.txt" delim=" "/>')
    w('\t\t</logfile>')
    w('\t</Output>')
    w('</febio_spec>')
    text = "\n".join(o) + "\n"

    os.makedirs(args.out, exist_ok=True)
    feb = os.path.join(args.out, f"{args.tag}.feb")
    with open(feb, "w", encoding="latin-1") as fh:
        fh.write(text)
    meta = dict(vars(args)); meta.update(
        omega_ref=omega, period_ref=period, dt=dt, n_steps=nsteps, n_nodes=len(m["nodes"]),
        n_ring_el=len(m["ring_el"]), n_fluid_el=len(m["fl_el"]), dims_nrs_nt_nrf=m["dims"],
        probe_nodes=dict(ring_outer_theta0=probes[0], ring_outer_theta90=probes[1],
                         ring_inner_theta0=probes[2], ring_inner_theta90=probes[3]),
        potential_coeffs=dict(A=Acoef, B=Bcoef), feb_sha256=hashlib.sha256(text.encode("latin-1")).hexdigest(),
        ring_velocity_expr=[rx, ry], fluid_velocity_expr=[fx, fy])
    with open(os.path.join(args.out, "model_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    return feb, meta


def selfcheck():
    """IC expressions vs direct polar evaluation, and kinematic consistency at the interface."""
    (rx, ry), (fx, fy), (A, Bc) = mode_expressions(V0, B, RO)
    worst = 0.0
    for th in (0.1, 0.5, 0.9, 1.3):
        for r in (0.95, 1.0, 1.05):
            vr = V0 * math.cos(2 * th); vt = V0 * math.sin(2 * th) * (-0.5 + 1.5 * (r - 1.0))
            vx = vr * math.cos(th) - vt * math.sin(th); vy = vr * math.sin(th) + vt * math.cos(th)
            X, Y = r * math.cos(th), r * math.sin(th)
            worst = max(worst, abs(evaluate(rx, X, Y) - vx), abs(evaluate(ry, X, Y) - vy))
        for r in (RO, 1.2, B):
            vr = 2 * (A * r - Bc / r ** 3) * math.cos(2 * th); vt = -2 * (A * r + Bc / r ** 3) * math.sin(2 * th)
            vx = vr * math.cos(th) - vt * math.sin(th); vy = vr * math.sin(th) + vt * math.cos(th)
            X, Y = r * math.cos(th), r * math.sin(th)
            worst = max(worst, abs(evaluate(fx, X, Y) - vx), abs(evaluate(fy, X, Y) - vy))
    # radial velocity continuity at the interface and zero normal velocity at the wall
    th = 0.7
    vr_if = 2 * (A * RO - Bc / RO ** 3) * math.cos(2 * th)
    vr_wall = 2 * (A * B - Bc / B ** 3) * math.cos(2 * th)
    return worst, vr_if - V0 * math.cos(2 * th), vr_wall


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", required=False)
    p.add_argument("--tag", default="ring")
    p.add_argument("--rho-f", type=float, default=None, help="fluid density; omit for the dry ring")
    p.add_argument("--nu-f", type=float, default=None, help="kinematic viscosity (default 1e-7)")
    p.add_argument("--k", type=float, default=1e8, help="fluid bulk modulus K [Pa]")
    p.add_argument("--e-mesh", type=float, default=1.0, help="mesh-regularising solid E [Pa]")
    p.add_argument("--rhoi", type=float, default=1.0)
    p.add_argument("--mesh-factor", type=int, default=1)
    p.add_argument("--steps-per-period", type=int, default=100)
    p.add_argument("--periods", type=float, default=4.0)
    p.add_argument("--ic", choices=["potential", "impulsive"], default="potential")
    p.add_argument("--ic-field", choices=["potential", "compatible"], default="potential",
                   help="potential: v1 field (incompatible with no-slip); compatible: potential + solenoidal no-slip correction (v2)")
    p.add_argument("--interface", choices=["noslip", "slip"], default="noslip",
                   help="fluid-ring interface: w=0 (no-slip, primary) or natural w.n=0/frictionless (slip diagnostic)")
    p.add_argument("--wall", choices=["noslip", "slip"], default="noslip", help="rigid outer wall: w=0 (no-slip) or natural (free-slip)")
    p.add_argument("--ic-mode", choices=["pernode", "math"], default="pernode",
                   help="pernode: constant IC per node (works on pristine v4.13); math: expression IC (needs the one-line patch)")
    p.add_argument("--dtol", type=float, default=1e-6); p.add_argument("--vtol", type=float, default=1e-6)
    p.add_argument("--ftol", type=float, default=1e-6); p.add_argument("--etol", type=float, default=1e-8)
    p.add_argument("--rtol", type=float, default=1e-10); p.add_argument("--max-refs", type=int, default=25)
    p.add_argument("--plot-stride", type=int, default=10)
    p.add_argument("--v0", type=float, default=1.5e-3, help="initial velocity amplitude V0 [m/s] (default 1.5e-3 = baseline; amplitude diagnostic only)")
    p.add_argument("--selfcheck", action="store_true")
    a = p.parse_args(argv)
    global V0
    V0 = a.v0
    if a.selfcheck:
        w, d1, d2 = selfcheck()
        print(f"IC expression max |expr - polar formula| = {w:.3e}; interface radial mismatch = {d1:.3e}; wall radial velocity = {d2:.3e}")
        return 0 if (w < 1e-12 and abs(d1) < 1e-12 and abs(d2) < 1e-12) else 1
    if not a.out:
        p.error("--out required")
    feb, meta = make(a)
    print(feb, json.dumps({k: meta[k] for k in ("n_nodes", "n_ring_el", "n_fluid_el", "dt", "n_steps", "feb_sha256")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
