#!/usr/bin/env python3
"""radialTorsionalRing (TORSION member) FEBio v4.13 input generator (RTR agent).  Stock FEBio only; nothing here patches FEBio.

PROVENANCE.  Geometry, materials, drives and mesh family from the user's own draft proposal (solids4foam PR #637 head 30ae00eb, tutorials/fluidSolidInteraction/radialTorsionalRing:
README Tables 1-2, system/{solid,fluid}/blockMeshDict, constant/{solid,fluid}, 0/*, verification/scripts/rt_reference.py; manuscript PR #25 head 97f9efef).  Copies of the reference
scripts are vendored (sha256 in RTR_EQUIVALENCE_AUDIT.md).  Plane-strain annulus modelled as a z-slab one hex8 layer thick (dz = 0.1 m, as the tutorial), full 360 degrees.
  solid ring  A=0.5 <= R <= B=1   : FEBio `Mooney-Rivlin` c1 = mu/2, c2 = 0, k = K, pressure_model NIKE3D  (W = mu/2 (Ibar1-3) + K/4 (J^2-1-2 ln J); mu = E/(2(1+nu)) = 7692.3077,
                                    K = E/(3(1-2 nu)) = 16666.667; this is the proposal's `neoHookeanElastic` law, NOT FEBio's `neo-Hookean`), rho_s = 1000, large strain (total Lagrangian)
  fluid       B <= r <= c=2       : `fluid-FSI` mixture (coupled) or `fluid` (fixed-geometry): rho_f = 1000, mu_f = 300, weakly compressible K_f (default 1e7, declared), laminar, convective
  mesh        M1 nth=64, nrs=8, nrf=16 (512 + 1024 cells); M2 128/16/32 (2048 + 4096); M3 256/32/64; uniform radial spacing; nodes exactly on the circles (polygonal faces, as the tutorial's
              blockMesh arcs at the vertices); interface conformal (shared nodes).
  modes       law1     2x2 hex8 block, all boundary nodes prescribed to a homogeneous F (plane strain), centre column free: Cauchy stress vs the hand formula (RT-A law identity)
              rotcheck one hex8 under the exact rotation expressions used for the coupled run (RT-A syntax/time check of the `math` BC with the C-infinity ramp)
              twist    RT-A structural test: static annulus, inner circle rotated by alpha_A(t) = amp S(t) sin(3t), outer circle clamped (alpha_B = 0, radially clamped), plane strain
              fluid    RT-B: fixed-geometry fluid annulus 1<=r<=2, wall circle moves tangentially with the tabulated omega_w(t) = omega_B(t) of the proposal's coupled torsion reference,
                       outer circle natural (zero traction) with p = 0 (`zero fluid dilatation`), K_f weak compressibility
              fsi      RT-C: coupled ring + fluid, inner solid circle prescribed rotation alpha_A(t) = 0.1 S(t) sin(3t) (radially clamped = rigid rotation), no-slip conformal interface, fluid
                       outer circle natural + p = 0, `fluid-FSI traction` on the interface, regularising skeleton neo-Hookean (E_mesh, nu 0, density 0)
Rotation BC (FEBio `math` expression of X, Y, t; FEBio's math parser string-substitution `name=expr;` is used; C-infinity ramp S(tau) = f/(f+g), f = exp(-1/tau), g = exp(-1/(1-tau)),
tau clamped to [0.005, 0.995]: S differs from the exact C-infinity step by < 1e-86):  u_x = (cos a - 1) X - sin a Y,  u_y = sin a X + (cos a - 1) Y.
"""
import argparse, hashlib, json, math, os, sys
import numpy as np

A_IN, B_OUT, C_OUT, DZ = 0.5, 1.0, 2.0, 0.1
E_S, NU_S, RHO_S = 2.0e4, 0.3, 1000.0
MU_S = E_S / (2 * (1 + NU_S)); K_S = E_S / (3 * (1 - 2 * NU_S))
RHO_F, MU_F = 1000.0, 300.0
MESH = {"M0": (16, 2, 4), "M1": (64, 8, 16), "M2": (128, 16, 32), "M3": (256, 32, 64)}


def ramp_S_expr():
    TC = "min(max(t,0.005),0.995)"
    return f"exp(-1/{TC})/(exp(-1/{TC})+exp(-1/(1-{TC})))"


def rot_exprs(amp):
    """(ux_expr, uy_expr) rigid rotation by a(t) = amp S(t) sin(3 t) about the z axis (reference position X, Y)."""
    sub = f"SSS={ramp_S_expr()};AAA={amp:.12g}*SSS*sin(3*t);"
    return (sub + "(cos(AAA)-1)*X-sin(AAA)*Y", sub + "sin(AAA)*X+(cos(AAA)-1)*Y")


def lst(ids):
    return ",".join(str(i) for i in ids)


class Grid:
    """Structured polar hex8 grids sharing the interface circle.  Node ids 1-based; z-layers k = 0, 1."""
    def __init__(self, nth, nrs, nrf, solid=True, fluid=True):
        self.nth, self.nrs, self.nrf = nth, nrs, nrf
        self.nodes, self.idx = [], {}
        self.solid_el, self.fluid_el = [], []
        dth = 2 * math.pi / nth

        def add(key, x, y, z):
            if key in self.idx:
                return self.idx[key]
            self.nodes.append((x, y, z)); self.idx[key] = len(self.nodes); return self.idx[key]
        if solid:
            for k in (0, 1):
                for j in range(nth):
                    for i in range(nrs + 1):
                        r = A_IN + (B_OUT - A_IN) * i / nrs
                        add(("s", i, j, k), r * math.cos(j * dth), r * math.sin(j * dth), k * DZ)
        if fluid:
            for k in (0, 1):
                for j in range(nth):
                    for i in range(nrf + 1):
                        r = B_OUT + (C_OUT - B_OUT) * i / nrf
                        key = ("s", nrs, j, k) if (i == 0 and solid) else ("f", i, j, k)
                        add(key, r * math.cos(j * dth), r * math.sin(j * dth), k * DZ)
        self.has_solid, self.has_fluid = solid, fluid

        def hex8(kind, i, j, k0=("s", "f")):
            jn = (j + 1) % nth

            def g(ii, jj, kk):
                if kind == "f" and ii == 0 and solid:
                    return self.idx[("s", nrs, jj, kk)]
                return self.idx[(kind, ii, jj, kk)]
            return (g(i, j, 0), g(i + 1, j, 0), g(i + 1, jn, 0), g(i, jn, 0), g(i, j, 1), g(i + 1, j, 1), g(i + 1, jn, 1), g(i, jn, 1))
        if solid:
            for j in range(nth):
                for i in range(nrs):
                    self.solid_el.append(hex8("s", i, j))
        if fluid:
            for j in range(nth):
                for i in range(nrf):
                    self.fluid_el.append(hex8("f", i, j))

    def ring(self, kind, i, layers=(0, 1)):
        return sorted(self.idx[(kind, i, j, k)] for j in range(self.nth) for k in layers)

    def iface_ring(self, layers=(0, 1)):
        return sorted(self.idx[("s", self.nrs, j, k)] for j in range(self.nth) for k in layers) if self.has_solid else self.ring("f", 0, layers)

    def iface_faces(self):
        n = self.nth; kind, i = ("s", self.nrs) if self.has_solid else ("f", 0)
        return [(self.idx[(kind, i, j, 0)], self.idx[(kind, i, (j + 1) % n, 0)], self.idx[(kind, i, (j + 1) % n, 1)], self.idx[(kind, i, j, 1)]) for j in range(n)]


def spec_header(tag, module, note):
    return ['<?xml version="1.0" encoding="ISO-8859-1"?>', f'<!-- {note} Generated by tools/make_rtr_feb.py (RTR agent) tag={tag} -->',
            '<febio_spec version="4.0">', f'\t<Module type="{module}"/>']


def build(a):
    o = []; w = o.append
    mode = a.mode
    nth, nrs, nrf = MESH[a.mesh] if a.mesh else (a.nth, a.nrs, a.nrf)
    meta = dict(vars(a)); meta.update(nth=nth, nrs=nrs, nrf=nrf)
    if mode == "law1":
        # 2 x 2 x 1 elements (unit square, one layer), plane strain F33 = 1; ALL boundary nodes prescribed to a homogeneous F, the single interior column of nodes free:
        # the affine field is then the exact discrete equilibrium (v1 of this test left the free faces unconstrained: a design error, see RTR_PILOT_RESULTS.md)
        Fm = np.array([[1.1, 0.4, 0.0], [0.1, 0.9, 0.0], [0.0, 0.0, 1.0]])
        nid = lambda ix, iy, iz: 1 + ix + 3 * iy + 9 * iz
        nodes = [(0.5 * ix, 0.5 * iy, float(iz)) for iz in (0, 1) for iy in (0, 1, 2) for ix in (0, 1, 2)]
        els = [(nid(ix, iy, 0), nid(ix + 1, iy, 0), nid(ix + 1, iy + 1, 0), nid(ix, iy + 1, 0), nid(ix, iy, 1), nid(ix + 1, iy, 1), nid(ix + 1, iy + 1, 1), nid(ix, iy + 1, 1)) for iy in (0, 1) for ix in (0, 1)]
        alln = list(range(1, 19)); ends = sorted(n for n in alln if n not in (nid(1, 1, 0), nid(1, 1, 1))); mid = [nid(1, 1, 0), nid(1, 1, 1)]
        for l in spec_header(a.tag, "solid", "RT-A law identity: 2x2 hex8, homogeneous F."):
            w(l)
        w('\t<Control><analysis>STATIC</analysis><time_steps>5</time_steps><step_size>0.2</step_size>'
          '<solver type="solid"><symmetric_stiffness>symmetric</symmetric_stiffness><max_refs>25</max_refs><dtol>1e-10</dtol><etol>1e-20</etol><rtol>0</rtol><lstol>0.9</lstol>'
          '<min_residual>1e-30</min_residual><qn_method type="full Newton"/></solver></Control>')
        w(f'\t<Material><material id="1" name="Ring" type="Mooney-Rivlin"><density>{RHO_S:g}</density><c1>{MU_S / 2:.15g}</c1><c2>0</c2><k>{K_S:.15g}</k>'
          '<pressure_model>NIKE3D</pressure_model></material></Material>')
        w('\t<Mesh><Nodes name="N">')
        for n, (x, y, z) in enumerate(nodes, 1):
            w(f'\t\t<node id="{n}">{x},{y},{z}</node>')
        w('\t</Nodes><Elements type="hex8" name="Ring">')
        for i, e in enumerate(els, 1):
            w(f'\t\t<elem id="{i}">{lst(e)}</elem>')
        w(f'\t</Elements><NodeSet name="Ends">{lst(ends)}</NodeSet><NodeSet name="All">{lst(alln)}</NodeSet></Mesh>')
        w('\t<MeshDomains><SolidDomain name="Ring" mat="Ring"/></MeshDomains><Boundary>')
        for dof, row in (("x", 0), ("y", 1)):
            e = "+".join(f"({Fm[row, c] - (1 if c == row else 0):.15g})*{'XYZ'[c]}" for c in range(3))
            w(f'\t\t<bc name="u{dof}" node_set="Ends" type="prescribed displacement"><dof>{dof}</dof><value type="math">t*({e})</value><relative>0</relative></bc>')
        w('\t\t<bc name="uz" node_set="All" type="zero displacement"><x_dof>0</x_dof><y_dof>0</y_dof><z_dof>1</z_dof></bc>')
        w('\t</Boundary><Output><plotfile type="febio"><var type="displacement"/><var type="stress"/></plotfile><logfile>'
          '<element_data data="sx;sy;sz;sxy;syz;sxz;J" file="stress.txt" delim=" ">1,2,3,4</element_data>'
          f'<node_data data="x;y;z;ux;uy;uz" file="nodes.txt" delim=" ">{lst(alln)}</node_data></logfile></Output></febio_spec>')
        meta.update(F=Fm.tolist(), note="u = t (F-I) X, t = 0.2..1.0: F_k = I + 0.2k (Fm - I); stress logged at each of the 5 steps")
        return o, None, meta
    if mode == "rotcheck":
        nodes = [(x, y, z) for z in (0.0, 0.1) for (x, y) in ((0.5, 0), (1.0, 0), (1.0, 0.1), (0.5, 0.1))]
        ex, ey = rot_exprs(0.5)
        for l in spec_header(a.tag, "solid", "RT-A rotation-expression check: single hex8 prescribed with the exact rotation expressions."):
            w(l)
        w('\t<Control><analysis>STATIC</analysis><time_steps>40</time_steps><step_size>0.05</step_size>'
          '<solver type="solid"><symmetric_stiffness>symmetric</symmetric_stiffness><max_refs>25</max_refs><dtol>1e-10</dtol><etol>1e-20</etol><rtol>0</rtol><lstol>0.9</lstol>'
          '<min_residual>1e-30</min_residual><qn_method type="full Newton"/></solver></Control>')
        w(f'\t<Material><material id="1" name="Ring" type="Mooney-Rivlin"><density>{RHO_S:g}</density><c1>{MU_S / 2:.15g}</c1><c2>0</c2><k>{K_S:.15g}</k><pressure_model>NIKE3D</pressure_model></material></Material>')
        w('\t<Mesh><Nodes name="N">')
        for n, (x, y, z) in enumerate(nodes, 1):
            w(f'\t\t<node id="{n}">{x},{y},{z}</node>')
        w('\t</Nodes><Elements type="hex8" name="Ring"><elem id="1">1,2,3,4,5,6,7,8</elem></Elements><NodeSet name="All">1,2,3,4,5,6,7,8</NodeSet></Mesh>')
        w('\t<MeshDomains><SolidDomain name="Ring" mat="Ring"/></MeshDomains><Boundary>')
        w(f'\t\t<bc name="ux" node_set="All" type="prescribed displacement"><dof>x</dof><value type="math">{ex}</value><relative>0</relative></bc>')
        w(f'\t\t<bc name="uy" node_set="All" type="prescribed displacement"><dof>y</dof><value type="math">{ey}</value><relative>0</relative></bc>')
        w('\t\t<bc name="uz" node_set="All" type="zero displacement"><x_dof>0</x_dof><y_dof>0</y_dof><z_dof>1</z_dof></bc>')
        w('\t</Boundary><Output><logfile><node_data data="x;y;ux;uy" file="nodes.txt" delim=" ">1,2,3,4</node_data></logfile></Output></febio_spec>')
        return o, None, meta
    # ----------------------------------------------------------------- annulus modes
    g = Grid(nth, nrs, nrf, solid=(mode in ("twist", "fsi") and not (mode == "fsi" and a.no_solid)), fluid=(mode in ("fluid", "fsi")))
    nodes = g.nodes
    allids = list(range(1, len(nodes) + 1))
    dt, tend = a.dt, a.t_end
    nsteps = int(round(tend / dt))
    solid_ids = sorted({n for e in g.solid_el for n in e}); fluid_ids = sorted({n for e in g.fluid_el for n in e})
    inner = g.ring("s", 0) if g.has_solid else []
    w("")
    o.pop()
    if mode == "twist":
        ex, ey = rot_exprs(a.amp)
        for l in spec_header(a.tag, "solid", f"RT-A static torsion of the annulus (amp {a.amp} rad inner, outer clamped)."):
            w(l)
        w(f'\t<Control><analysis>STATIC</analysis><time_steps>{nsteps}</time_steps><step_size>{dt:g}</step_size><plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range>'
          '<plot_level>PLOT_MAJOR_ITRS</plot_level><plot_stride>%d</plot_stride><output_level>OUTPUT_MAJOR_ITRS</output_level>' % a.plot_stride)
        w('\t\t<time_stepper type="default"><max_retries>10</max_retries><opt_iter>10</opt_iter><dtmin>0</dtmin><dtmax>%g</dtmax><aggressiveness>0</aggressiveness><cutback>0.5</cutback><dtforce>0</dtforce></time_stepper>' % dt)
        w('\t\t<solver type="solid"><symmetric_stiffness>symmetric</symmetric_stiffness><max_refs>25</max_refs><dtol>1e-10</dtol><etol>1e-20</etol><rtol>0</rtol><lstol>0.9</lstol>'
          '<min_residual>1e-30</min_residual><qn_method type="full Newton"/></solver></Control>')
        w(f'\t<Material><material id="1" name="Ring" type="Mooney-Rivlin"><density>{RHO_S:g}</density><c1>{MU_S / 2:.15g}</c1><c2>0</c2><k>{K_S:.15g}</k><pressure_model>NIKE3D</pressure_model></material></Material>')
    # mesh block for all annulus modes
    pre = []
    if mode == "fluid":
        for l in spec_header(a.tag, "fluid", "RT-B fixed-geometry fluid annulus with tabulated tangential wall velocity."):
            w(l)
        w(f'\t<Control><analysis>DYNAMIC</analysis><time_steps>{nsteps}</time_steps><step_size>{dt:g}</step_size><plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range>'
          f'<plot_level>PLOT_MAJOR_ITRS</plot_level><output_level>OUTPUT_MAJOR_ITRS</output_level><plot_stride>{a.plot_stride}</plot_stride>')
        w('\t\t<solver type="fluid"><symmetric_stiffness>non-symmetric</symmetric_stiffness><equation_scheme>block</equation_scheme>'
          f'<lstol>0.9</lstol><lsmin>0.01</lsmin><lsiter>5</lsiter><max_refs>{a.max_refs}</max_refs><check_zero_diagonal>0</check_zero_diagonal><reform_each_time_step>1</reform_each_time_step>'
          f'<diverge_reform>1</diverge_reform><min_residual>{a.min_res:g}</min_residual><max_residual>1e+20</max_residual>'
          f'<vtol>{a.vtol:g}</vtol><ftol>{a.ftol:g}</ftol><etol>{a.etol:g}</etol><rtol>{a.rtol:g}</rtol><rhoi>{a.rhoi:g}</rhoi><predictor>0</predictor><min_volume_ratio>0</min_volume_ratio>'
          '<qn_method type="full Newton"/></solver></Control>')
        w(f'\t<Material><material id="1" name="Fluid" type="fluid"><density>{RHO_F:g}</density><k>{a.k:g}</k><viscous type="Newtonian fluid"><kappa>0</kappa><mu>{a.mu_f:g}</mu></viscous></material></Material>')
    if mode == "fsi":
        ex, ey = rot_exprs(0.1)
        for l in spec_header(a.tag, "fluid-FSI", "RT-C coupled torsion member."):
            w(l)
        w(f'\t<Control><analysis>DYNAMIC</analysis><time_steps>{nsteps}</time_steps><step_size>{dt:g}</step_size><plot_zero_state>0</plot_zero_state><plot_range>0,-1</plot_range>'
          f'<plot_level>PLOT_MAJOR_ITRS</plot_level><output_level>OUTPUT_MAJOR_ITRS</output_level><plot_stride>{a.plot_stride}</plot_stride>')
        w('\t\t<solver type="fluid-FSI"><symmetric_stiffness>non-symmetric</symmetric_stiffness><equation_scheme>staggered</equation_scheme><equation_order>default</equation_order><optimize_bw>0</optimize_bw>'
          f'<lstol>0.9</lstol><lsmin>0.01</lsmin><lsiter>5</lsiter><max_refs>{a.max_refs}</max_refs><check_zero_diagonal>0</check_zero_diagonal><reform_each_time_step>1</reform_each_time_step><diverge_reform>1</diverge_reform>'
          f'<min_residual>{a.min_res:g}</min_residual><max_residual>1e+20</max_residual><dtol>{a.dtol:g}</dtol><vtol>{a.vtol:g}</vtol><ftol>{a.ftol:g}</ftol><etol>{a.etol:g}</etol><rtol>{a.rtol:g}</rtol>'
          f'<rhoi>{a.rhoi:g}</rhoi><predictor>0</predictor><min_volume_ratio>0</min_volume_ratio><order>2</order><qn_method type="full Newton"/></solver></Control>')
        w('\t<Material>')
        if a.no_solid:
            pass
        elif a.solid == "svk":
            w(f'\t\t<material id="1" name="Ring" type="isotropic elastic"><density>{RHO_S:g}</density><E>{E_S:g}</E><v>{NU_S:g}</v></material>')
        else:
            w(f'\t\t<material id="1" name="Ring" type="Mooney-Rivlin"><density>{RHO_S:g}</density><c1>{MU_S / 2:.15g}</c1><c2>0</c2><k>{K_S:.15g}</k><pressure_model>NIKE3D</pressure_model></material>')
        w('\t\t<material id="%d" name="Fluid" type="fluid-FSI">' % (1 if a.no_solid else 2))
        w(f'\t\t\t<solid type="neo-Hookean"><density>0</density><E>{a.e_mesh:g}</E><v>0</v></solid>')
        w(f'\t\t\t<fluid type="fluid"><density>{RHO_F:g}</density><k>{a.k:g}</k><viscous type="Newtonian fluid"><kappa>0</kappa><mu>{a.mu_f:g}</mu></viscous></fluid></material>')
        w('\t</Material>')
    w('\t<Mesh>\n\t\t<Nodes name="AllNodes">')
    for n, (x, y, z) in enumerate(nodes, 1):
        w(f'\t\t\t<node id="{n}">{x:.15g},{y:.15g},{z:.15g}</node>')
    w('\t\t</Nodes>')
    eid = 0
    if g.has_solid:
        w('\t\t<Elements type="hex8" name="Ring">')
        for e in g.solid_el:
            eid += 1; w(f'\t\t\t<elem id="{eid}">{lst(e)}</elem>')
        w('\t\t</Elements>')
    if g.has_fluid:
        w('\t\t<Elements type="hex8" name="Fluid">')
        for e in g.fluid_el:
            eid += 1; w(f'\t\t\t<elem id="{eid}">{lst(e)}</elem>')
        w('\t\t</Elements>')
    outer_f = g.ring("f", nrf) if g.has_fluid else []
    iface = g.iface_ring()
    sets = {"AllN": allids}
    if mode == "twist":
        sets["Inner"] = inner; sets["Outer"] = g.ring("s", nrs)
        ray = [g.idx[("s", i, 0, 0)] for i in range(nrs + 1)]; ray2 = [g.idx[("s", i, nth // 8, 0)] for i in range(nrs + 1)]
        sets["Ray0"] = ray; sets["Ray45"] = ray2
    if mode == "fluid":
        sets["Wall"] = iface; sets["OuterF"] = outer_f; sets["FluidAll"] = fluid_ids
        sets["Mon0"] = [n for n in fluid_ids if abs(nodes[n - 1][2]) < 1e-12]
        sets["MonTop"] = [g.idx[("f", i, 0, 1)] if not (i == 0 and g.has_solid) else g.idx[("s", nrs, 0, 1)] for i in range(nrf + 1)]
    if mode == "fsi":
        sets["Iface"] = iface; sets["OuterF"] = outer_f; sets["FluidAll"] = fluid_ids
        if g.has_solid:
            sets["Inner"] = inner; sets["SolidAll"] = solid_ids
            sets["Mon0S"] = [n for n in solid_ids if abs(nodes[n - 1][2]) < 1e-12]
        sets["FluidFree"] = sorted(set(fluid_ids) - set(iface)); sets["FluidInt"] = sorted(set(fluid_ids) - set(iface))
        sets["Mon0F"] = [n for n in fluid_ids if abs(nodes[n - 1][2]) < 1e-12]
        sets["MonTop"] = [g.idx[("f", i, 0, 1)] if (i > 0 or not g.has_solid) else g.idx[("s", nrs, 0, 1)] for i in range(nrf + 1)]
    for nm, ids in sets.items():
        w(f'\t\t<NodeSet name="{nm}">{lst(ids)}</NodeSet>')
    wallpair = {}
    if mode == "fluid":
        # one node set per wall angular position (both z layers): tangential velocity = omega(t) (-Y, X) via per-node amplitude x load curve
        for j in range(nth):
            pair = sorted([g.idx[("f", 0, j, 0)], g.idx[("f", 0, j, 1)]]); wallpair[j] = pair
            w(f'\t\t<NodeSet name="W{j}">{lst(pair)}</NodeSet>')
    if g.has_fluid:
        w('\t\t<Surface name="Interface">')
        for i, q in enumerate(g.iface_faces(), 1):
            w(f'\t\t\t<quad4 id="{i}">{lst(q)}</quad4>')
        w('\t\t</Surface>')
    w('\t</Mesh>\n\t<MeshDomains>')
    if g.has_solid:
        w('\t\t<SolidDomain name="Ring" mat="Ring"/>')
    if g.has_fluid:
        w('\t\t<SolidDomain name="Fluid" mat="%s"/>' % ("Fluid"))
    w('\t</MeshDomains>\n\t<Boundary>')
    if mode == "twist":
        w(f'\t\t<bc name="ux_in" node_set="Inner" type="prescribed displacement"><dof>x</dof><value type="math">{ex}</value><relative>0</relative></bc>')
        w(f'\t\t<bc name="uy_in" node_set="Inner" type="prescribed displacement"><dof>y</dof><value type="math">{ey}</value><relative>0</relative></bc>')
        w('\t\t<bc name="fix_out_x" node_set="Outer" type="prescribed displacement"><dof>x</dof><value>0</value><relative>0</relative></bc>')
        w('\t\t<bc name="fix_out_y" node_set="Outer" type="prescribed displacement"><dof>y</dof><value>0</value><relative>0</relative></bc>')
        w('\t\t<bc name="uz" node_set="AllN" type="zero displacement"><x_dof>0</x_dof><y_dof>0</y_dof><z_dof>1</z_dof></bc>')
    if mode == "fluid":
        w('\t\t<bc name="wz" node_set="AllN" type="zero fluid velocity"><wx_dof>0</wx_dof><wy_dof>0</wy_dof><wz_dof>1</wz_dof></bc>')
        for j in range(nth):
            th = 2 * math.pi * j / nth; X, Y = B_OUT * math.cos(th), B_OUT * math.sin(th)
            w(f'\t\t<bc name="wx{j}" node_set="W{j}" type="prescribed fluid velocity"><dof>wx</dof><value lc="1">{-Y:.16g}</value><relative>0</relative></bc>')
            w(f'\t\t<bc name="wy{j}" node_set="W{j}" type="prescribed fluid velocity"><dof>wy</dof><value lc="1">{X:.16g}</value><relative>0</relative></bc>')
        w('\t\t<bc name="p_out" node_set="OuterF" type="zero fluid dilatation"/>')
    if mode == "fsi":
        w('\t\t<bc name="uz" node_set="AllN" type="zero displacement"><x_dof>0</x_dof><y_dof>0</y_dof><z_dof>1</z_dof></bc>')
        rs = "Iface" if a.no_solid else "Inner"
        w(f'\t\t<bc name="ux_in" node_set="{rs}" type="prescribed displacement"><dof>x</dof><value type="math">{ex}</value><relative>0</relative></bc>')
        w(f'\t\t<bc name="uy_in" node_set="{rs}" type="prescribed displacement"><dof>y</dof><value type="math">{ey}</value><relative>0</relative></bc>')
        w('\t\t<bc name="wz" node_set="FluidAll" type="zero fluid velocity"><wx_dof>0</wx_dof><wy_dof>0</wy_dof><wz_dof>1</wz_dof></bc>')
        w('\t\t<bc name="wxy_if" node_set="Iface" type="zero fluid velocity"><wx_dof>1</wx_dof><wy_dof>1</wy_dof><wz_dof>0</wz_dof></bc>')
        w('\t\t<bc name="p_out" node_set="OuterF" type="zero fluid dilatation"/>')
        if a.clamp_outer_mesh:
            w('\t\t<bc name="mesh_out" node_set="OuterF" type="zero displacement"><x_dof>1</x_dof><y_dof>1</y_dof><z_dof>0</z_dof></bc>')
    w('\t</Boundary>')
    if mode == "fsi" and not a.no_solid:
        w('\t<Loads>\n\t\t<surface_load name="FSItraction" surface="Interface" type="fluid-FSI traction"><shell_bottom>0</shell_bottom></surface_load>\n\t</Loads>')
    if mode == "fluid":
        TT, AL, OM = np.loadtxt(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reference_data", "rtr", "omega_B_table_dt1e-3.csv")).T
        OMs = np.interp(TT + a.bc_time_shift, TT, OM) if a.bc_time_shift else OM      # DIAGNOSTIC: table(tau) = omega(tau + shift) (FEBio evaluates the BC at t_n - dt/2, see RTR_PILOT_RESULTS.md)
        pts = "".join(f"<pt>{t:.10g},{om:.17g}</pt>" for t, om in zip(TT, OMs))
        w('\t<LoadData><load_controller id="1" name="omega_w" type="loadcurve"><interpolate>LINEAR</interpolate><extend>CONSTANT</extend><points>' + pts + '</points></load_controller></LoadData>')
    w('\t<Output>')
    if mode == "twist":
        w('\t\t<plotfile type="febio"><var type="displacement"/><var type="stress"/><var type="reaction forces"/></plotfile>\n\t\t<logfile>')
        w(f'\t\t\t<node_data data="x;y;ux;uy;Rx;Ry" file="inner.txt" delim=" ">{lst(inner)}</node_data>')
        w(f'\t\t\t<node_data data="x;y;ux;uy" file="ray0.txt" delim=" ">{lst(sets["Ray0"])}</node_data>')
        w(f'\t\t\t<node_data data="x;y;ux;uy" file="ray45.txt" delim=" ">{lst(sets["Ray45"])}</node_data>')
        w(f'\t\t\t<node_data data="x;y;ux;uy;Rx;Ry" file="outer.txt" delim=" ">{lst(sets["Outer"])}</node_data>')
        w('\t\t</logfile>')
    if mode == "fluid":
        w('\t\t<plotfile type="febio"><var type="fluid pressure"/><var type="fluid velocity"/></plotfile>\n\t\t<logfile>')
        w(f'\t\t\t<node_data data="wx;wy;ef" file="Mon0.txt" delim=" ">{lst(sets["Mon0"])}</node_data>')
        w(f'\t\t\t<node_data data="wx;wy;ef" file="MonTop.txt" delim=" ">{lst(sets["MonTop"])}</node_data>')
        w('\t\t</logfile>')
    if mode == "fsi":
        w('\t\t<plotfile type="febio"><var type="displacement"/><var type="fluid pressure"/><var type="relative fluid velocity"/><var type="stress"/></plotfile>\n\t\t<logfile>')
        w(f'\t\t\t<node_data data="x;y;vx;vy;wx;wy;ef" file="Mon0F.txt" delim=" ">{lst(sets["Mon0F"])}</node_data>')
        if g.has_solid:
            w(f'\t\t\t<node_data data="x;y;vx;vy" file="Mon0S.txt" delim=" ">{lst(sets["Mon0S"])}</node_data>')
        w(f'\t\t\t<node_data data="x;y;vx;vy;wx;wy;ef" file="MonTop.txt" delim=" ">{lst(sets["MonTop"])}</node_data>')
        w('\t\t</logfile>')
    w('\t</Output>\n</febio_spec>')
    return o, g, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True, choices=["law1", "rotcheck", "twist", "fluid", "fsi"])
    ap.add_argument("--out", required=True); ap.add_argument("--tag", required=True)
    ap.add_argument("--mesh", choices=sorted(MESH), default=None); ap.add_argument("--nth", type=int, default=16); ap.add_argument("--nrs", type=int, default=2); ap.add_argument("--nrf", type=int, default=4)
    ap.add_argument("--dt", type=float, default=0.01); ap.add_argument("--t-end", type=float, default=2.0)
    ap.add_argument("--amp", type=float, default=0.5, help="twist: inner rotation amplitude [rad]")
    ap.add_argument("--k", type=float, default=1e7); ap.add_argument("--e-mesh", type=float, default=1.0); ap.add_argument("--rhoi", type=float, default=1.0)
    ap.add_argument("--vtol", type=float, default=1e-6); ap.add_argument("--ftol", type=float, default=1e-6); ap.add_argument("--dtol", type=float, default=1e-6)
    ap.add_argument("--etol", type=float, default=0.0); ap.add_argument("--rtol", type=float, default=0.0); ap.add_argument("--max-refs", type=int, default=15)
    ap.add_argument("--plot-stride", type=int, default=50); ap.add_argument("--three-field", type=int, default=1, help="1: FEBio default three-field hex8 domain for the uncoupled Mooney-Rivlin; 0: `use_three_field_hex` off (standard displacement-based domain)"); ap.add_argument("--clamp-outer-mesh", type=int, default=0, help="fsi: hold the MESH (skeleton) displacement of the outer fluid circle at zero (fluid velocity/traction BC unchanged); amendment A2"); ap.add_argument("--no-solid", type=int, default=0, help="DIAGNOSTIC: fsi without the ring; the interface circle is prescribed to the rotation"); ap.add_argument("--mu-f", type=float, default=300.0, help="fluid dynamic viscosity (300 = the proposal; other values DIAGNOSTIC ONLY)"); ap.add_argument("--bc-time-shift", type=float, default=0.0, help="fluid mode DIAGNOSTIC: wall-velocity table advanced by this time"); ap.add_argument("--solid", choices=["mr", "svk"], default="mr", help="DIAGNOSTIC ONLY: svk replaces the proposal law")
    ap.add_argument("--min-res", type=float, default=1e-18, help="FEBio min_residual (squared residual norm below which a step is accepted); amendment A1 of the plan")
    a = ap.parse_args()
    lines, g, meta = build(a)
    text = "\n".join(lines) + "\n"
    if a.three_field == 0 and "Mooney-Rivlin" in text:
        text = text.replace('<SolidDomain name="Ring" mat="Ring"/>', '<SolidDomain name="Ring" mat="Ring" three_field="off"/>')
    os.makedirs(a.out, exist_ok=True)
    open(os.path.join(a.out, a.tag + ".feb"), "w", encoding="latin-1").write(text)
    meta.update(feb_sha256=hashlib.sha256(text.encode("latin-1")).hexdigest(), generator="tools/make_rtr_feb.py v1",
                febio_expected_sha256="f5ae1ed6f347f875c392f93d847f2369f437345a9b785296d46b0bd6fbb72f02")
    if g is not None:
        np.savez(os.path.join(a.out, "mesh_rtr.npz"), nodes=np.array(g.nodes), solid_el=np.array(g.solid_el, dtype=int) - 1 if g.solid_el else np.zeros((0, 8), int),
                 fluid_el=np.array(g.fluid_el, dtype=int) - 1 if g.fluid_el else np.zeros((0, 8), int), nth=g.nth, nrs=g.nrs, nrf=g.nrf,
                 iface=np.array(g.iface_ring((0,))) - 1, dz=DZ)
        meta.update(n_nodes=len(g.nodes), n_solid_el=len(g.solid_el), n_fluid_el=len(g.fluid_el),
                    mesh_sha256=hashlib.sha256((repr(g.nodes) + repr(g.solid_el) + repr(g.fluid_el)).encode()).hexdigest())
    json.dump(meta, open(os.path.join(a.out, "model_meta.json"), "w"), indent=1, default=str)
    print(os.path.join(a.out, a.tag + ".feb"), meta["feb_sha256"][:16], {k: meta[k] for k in meta if k in ("n_nodes", "n_solid_el", "n_fluid_el")})


if __name__ == "__main__":
    main()
