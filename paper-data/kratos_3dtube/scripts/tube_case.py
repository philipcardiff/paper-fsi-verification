#!/usr/bin/env python3
"""Write a Kratos input set for the 3dTube benchmark (quarter tube) at any level.

The block layout and divisions follow the solids4foam tutorial
(tutorials/fluidSolidInteraction/3dTube, system/{fluid,solid}/blockMeshDict), scaled by --factor:
- fluid: an O-grid with a square-like core block (10 x 10 cells) with corners (0,0), (a,0), (b,b),
  (0,a), a = 0.003414, b = 0.002828, and two outer blocks between the core and the wall, 5 cells
  radial and 10 around each (20 around the quarter); 80 cells along the 5 cm tube;
- wall: 4 cells through the 1 mm thickness, 20 around, 80 along (hexahedra);
- level factor 1, 2, 4 = solids4foam levels 1, 2, 3; factor 3 is the ratio-1.5 level above level 2.
Every fluid hexahedron is split into two prisms (the same 2D triangulation in every layer) and each
prism into three tetrahedra by the minimum-global-index rule, so the tetrahedral mesh is conforming.
Fluid and wall nodes coincide on the interface.

Usage:
  python tube_case.py <out dir> --factor 1 --dt 1e-4 --scheme euler [--end-time 0.02]
                      [--nl-tol 1e-9] [--vtk-every 1e-3]
--scheme euler: implicit Euler (BDF1) in the fluid (via MainKratos_tube.py), wall and mesh velocity;
--scheme bdf2:  BDF2 in all three.
"""
import argparse
import json
import math
import os

L = 0.05            # tube length (m)
R = 0.005           # inner radius (m)
H = 0.001           # wall thickness (m)
A_CORE = 0.003414   # core block corner on the axes (blockMeshDict vertices 1 and 3)
B_CORE = 0.002828   # core block corner on the diagonal (vertex 2)
RHO_F, MU_F = 1000.0, 1000.0 * 3e-6
RHO_S, E_S, NU_S = 1200.0, 3e5, 0.3
P_IN, T_PULSE, T_RAMP = 1333.2, 3e-3, 1e-4   # 0/fluid/time-series: 1.3332 to 3 ms, 0 at 3.1 ms
STATIONS = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045]
PROBE_OFFSET = 0.0005
FLUID_MDPA, SOLID_MDPA = "tube_fluid", "tube_structure"
TOL = 1e-12


def fmt(v):
    return f"{v:.17g}"


# --------------------------------------------------------------------------------------------- #
# Fluid cross-section (O-grid) and tetrahedra
# --------------------------------------------------------------------------------------------- #

def fluid_section(nc, nr):
    """2D nodes and quads of the fluid cross-section (shared points merged)."""
    pts, index = [], {}

    def node(x, y):
        x = 0.0 if abs(x) < 1e-15 else x
        y = 0.0 if abs(y) < 1e-15 else y
        key = (round(x, 12), round(y, 12))
        if key not in index:
            index[key] = len(pts)
            pts.append((x, y))
        return index[key]

    c0, c1, c2, c3 = (0.0, 0.0), (A_CORE, 0.0), (B_CORE, B_CORE), (0.0, A_CORE)

    def bilinear(u, v):
        return tuple((1 - u) * (1 - v) * c0[d] + u * (1 - v) * c1[d] + u * v * c2[d]
                     + (1 - u) * v * c3[d] for d in range(2))

    quads = []
    core = [[node(*bilinear(i / nc, j / nc)) for j in range(nc + 1)] for i in range(nc + 1)]
    for i in range(nc):
        for j in range(nc):
            quads.append([core[i][j], core[i + 1][j], core[i + 1][j + 1], core[i][j + 1]])

    def outer(p_from, p_to, th0, th1):
        # inner edge p_from -> p_to (straight), outer arc from th0 to th1, linear radial blend
        grid = []
        for t_i in range(nc + 1):
            t = t_i / nc
            inner = (p_from[0] + t * (p_to[0] - p_from[0]), p_from[1] + t * (p_to[1] - p_from[1]))
            th = th0 + t * (th1 - th0)
            wall = (R * math.cos(th), R * math.sin(th))
            grid.append([node(inner[0] + s / nr * (wall[0] - inner[0]),
                              inner[1] + s / nr * (wall[1] - inner[1])) for s in range(nr + 1)])
        for t_i in range(nc):
            for s in range(nr):
                quads.append([grid[t_i][s], grid[t_i][s + 1], grid[t_i + 1][s + 1], grid[t_i + 1][s]])

    outer(c1, c2, 0.0, math.pi / 4)          # block (1 4 5 2)
    outer(c2, c3, math.pi / 4, math.pi / 2)  # block (3 2 5 6)
    # counter-clockwise in the xy plane
    fixed = []
    for q in quads:
        area = 0.5 * sum(pts[q[k]][0] * pts[q[(k + 1) % 4]][1] - pts[q[(k + 1) % 4]][0] * pts[q[k]][1]
                         for k in range(4))
        fixed.append(q if area > 0 else q[::-1])
    return pts, fixed


def tet_volume(p):
    (x0, y0, z0), (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = p
    ax, ay, az = x1 - x0, y1 - y0, z1 - z0
    bx, by, bz = x2 - x0, y2 - y0, z2 - z0
    cx, cy, cz = x3 - x0, y3 - y0, z3 - z0
    return (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)) / 6.0


def prism_to_tets(v):
    """v = [a, b, c, a', b', c'] global ids (bottom triangle, top triangle above it).
    Minimum-index rule (Dompierre et al. 1999): conforming across shared quadrilateral faces."""
    # rotations and the top/bottom swap of the prism
    perms = [(0, 1, 2, 3, 4, 5), (1, 2, 0, 4, 5, 3), (2, 0, 1, 5, 3, 4),
             (3, 5, 4, 0, 2, 1), (4, 3, 5, 1, 0, 2), (5, 4, 3, 2, 1, 0)]
    imin = v.index(min(v))
    perm = next(p for p in perms if p[0] == imin)
    w = [v[k] for k in perm]  # w[0] is the minimum
    if min(w[1], w[5]) < min(w[2], w[4]):
        return [[w[0], w[1], w[2], w[5]], [w[0], w[1], w[5], w[4]], [w[0], w[4], w[5], w[3]]]
    return [[w[0], w[1], w[2], w[4]], [w[0], w[4], w[2], w[5]], [w[0], w[4], w[5], w[3]]]


def write_fluid(path, f):
    nc, nr, nz = 10 * f, 5 * f, 80 * f
    pts2, quads = fluid_section(nc, nr)
    n2 = len(pts2)

    def nid(k, i):
        return k * n2 + i + 1

    coords = {}
    for k in range(nz + 1):
        z = L * k / nz
        for i, (x, y) in enumerate(pts2):
            coords[nid(k, i)] = (x, y, z)
    tris2 = []
    for q in quads:
        tris2 += [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
    tets = []
    for k in range(nz):
        for t in tris2:
            v = [nid(k, t[0]), nid(k, t[1]), nid(k, t[2]), nid(k + 1, t[0]), nid(k + 1, t[1]), nid(k + 1, t[2])]
            for tet in prism_to_tets(v):
                if tet_volume([coords[n] for n in tet]) < 0:
                    tet[1], tet[2] = tet[2], tet[1]
                tets.append(tet)
    # boundary faces: faces of exactly one tetrahedron, oriented outwards
    faces = {}
    for tet in tets:
        for opp in range(4):
            face = [tet[j] for j in range(4) if j != opp]
            key = tuple(sorted(face))
            if key in faces:
                faces[key] = None
            else:
                faces[key] = (face, tet[opp])
    groups = {"inlet": [], "outlet": [], "wall": [], "symx": [], "symy": []}
    for key, val in faces.items():
        if val is None:
            continue
        face, opp = val
        p = [coords[n] for n in face]
        if tet_volume(p + [coords[opp]]) > 0:   # opposite node on the positive side: flip
            face = [face[0], face[2], face[1]]
        zs = [c[2] for c in p]
        rs = [math.hypot(c[0], c[1]) for c in p]
        if all(abs(z) < TOL for z in zs):
            groups["inlet"].append(face)
        elif all(abs(z - L) < TOL for z in zs):
            groups["outlet"].append(face)
        elif all(abs(r - R) < 1e-10 for r in rs):
            groups["wall"].append(face)
        elif all(abs(c[0]) < TOL for c in p):
            groups["symx"].append(face)
        elif all(abs(c[1]) < TOL for c in p):
            groups["symy"].append(face)
        else:
            raise RuntimeError(f"unclassified boundary face {face}")

    with open(path, "w") as fo:
        fo.write("Begin ModelPartData\nEnd ModelPartData\n\nBegin Properties 0\nEnd Properties\n\n")
        fo.write("Begin Properties 1\nEnd Properties\n\nBegin Nodes\n")
        for n in sorted(coords):
            fo.write(f"    {n} {fmt(coords[n][0])} {fmt(coords[n][1])} {fmt(coords[n][2])}\n")
        fo.write("End Nodes\n\nBegin Elements Element3D4N\n")
        for e, tet in enumerate(tets, 1):
            fo.write(f"    {e} 1 {tet[0]} {tet[1]} {tet[2]} {tet[3]}\n")
        fo.write("End Elements\n\n")
        cid, cond_ids = 0, {}
        for g in ("inlet", "outlet", "wall"):
            fo.write("Begin Conditions SurfaceCondition3D3N\n")
            cond_ids[g] = []
            for face in groups[g]:
                cid += 1
                cond_ids[g].append(cid)
                fo.write(f"    {cid} 0 {face[0]} {face[1]} {face[2]}\n")
            fo.write("End Conditions\n\n")

        def nodes_where(pred):
            return [n for n in sorted(coords) if pred(coords[n])]

        smp(fo, "Parts_Fluid", sorted(coords), range(1, len(tets) + 1))
        smp(fo, "Inlet", nodes_where(lambda c: abs(c[2]) < TOL), conditions=cond_ids["inlet"])
        smp(fo, "Outlet", nodes_where(lambda c: abs(c[2] - L) < TOL), conditions=cond_ids["outlet"])
        smp(fo, "FluidInterface", nodes_where(lambda c: abs(math.hypot(c[0], c[1]) - R) < 1e-10),
            conditions=cond_ids["wall"])
        smp(fo, "SymmetryX", nodes_where(lambda c: abs(c[0]) < TOL))
        smp(fo, "SymmetryY", nodes_where(lambda c: abs(c[1]) < TOL))
        smp(fo, "MeshEnds", nodes_where(lambda c: abs(c[2]) < TOL or abs(c[2] - L) < TOL))
    vol = sum(tet_volume([coords[n] for n in tet]) for tet in tets)
    return {"fluid_nodes": len(coords), "fluid_tets": len(tets), "fluid_hex_equivalent": len(quads) * nz,
            "fluid_volume": vol, "boundary_faces": {g: len(v) for g, v in groups.items()},
            "section_nodes": n2, "around_wall": 2 * nc, "radial_outer": nr, "along": nz}


# --------------------------------------------------------------------------------------------- #
# Wall (hexahedra)
# --------------------------------------------------------------------------------------------- #

def hex_volume_sign(p):
    # sign of the Jacobian at the centre for the node order bottom (0-3) and top (4-7) faces
    def sub(a, b):
        return [a[i] - b[i] for i in range(3)]
    dx = sub([(p[1][i] + p[2][i] + p[5][i] + p[6][i]) / 4 for i in range(3)],
             [(p[0][i] + p[3][i] + p[4][i] + p[7][i]) / 4 for i in range(3)])
    dy = sub([(p[2][i] + p[3][i] + p[6][i] + p[7][i]) / 4 for i in range(3)],
             [(p[0][i] + p[1][i] + p[4][i] + p[5][i]) / 4 for i in range(3)])
    dz = sub([(p[4][i] + p[5][i] + p[6][i] + p[7][i]) / 4 for i in range(3)],
             [(p[0][i] + p[1][i] + p[2][i] + p[3][i]) / 4 for i in range(3)])
    return (dx[0] * (dy[1] * dz[2] - dy[2] * dz[1]) - dx[1] * (dy[0] * dz[2] - dy[2] * dz[0])
            + dx[2] * (dy[0] * dz[1] - dy[1] * dz[0]))


def write_solid(path, f, node_offset, elem_offset):
    nt, nw, nz = 20 * f, 4 * f, 80 * f

    def nid(l, m, k):  # layer, radial, around
        return node_offset + 1 + (l * (nw + 1) + m) * (nt + 1) + k

    coords = {}
    for l in range(nz + 1):
        z = L * l / nz
        for m in range(nw + 1):
            r = R + H * m / nw
            for k in range(nt + 1):
                th = 0.5 * math.pi * k / nt
                x, y = r * math.cos(th), r * math.sin(th)
                if k == nt:
                    x = 0.0
                if k == 0:
                    y = 0.0
                coords[nid(l, m, k)] = (x, y, z)
    hexes = []
    for l in range(nz):
        for m in range(nw):
            for k in range(nt):
                h = [nid(l, m, k), nid(l, m + 1, k), nid(l, m + 1, k + 1), nid(l, m, k + 1),
                     nid(l + 1, m, k), nid(l + 1, m + 1, k), nid(l + 1, m + 1, k + 1), nid(l + 1, m, k + 1)]
                if hex_volume_sign([coords[n] for n in h]) < 0:
                    h = [h[0], h[3], h[2], h[1], h[4], h[7], h[6], h[5]]
                hexes.append(h)
    inner = []
    for l in range(nz):
        for k in range(nt):
            # inner face (m = 0) as two triangles: the FSI coupling skin in 3D needs triangles
            q = [nid(l, 0, k), nid(l + 1, 0, k), nid(l + 1, 0, k + 1), nid(l, 0, k + 1)]
            inner += [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
    with open(path, "w") as fo:
        fo.write("Begin ModelPartData\nEnd ModelPartData\n\nBegin Properties 0\nEnd Properties\n\nBegin Nodes\n")
        for n in sorted(coords):
            fo.write(f"    {n} {fmt(coords[n][0])} {fmt(coords[n][1])} {fmt(coords[n][2])}\n")
        fo.write("End Nodes\n\nBegin Elements SmallDisplacementElement3D8N\n")
        eids = []
        for e, h in enumerate(hexes, elem_offset + 1):
            eids.append(e)
            fo.write(f"    {e} 0 " + " ".join(str(n) for n in h) + "\n")
        fo.write("End Elements\n\nBegin Conditions SurfaceCondition3D3N\n")
        for c, q in enumerate(inner, 1):
            fo.write(f"    {c} 0 " + " ".join(str(n) for n in q) + "\n")
        fo.write("End Conditions\n\n")

        def nodes_where(pred):
            return [n for n in sorted(coords) if pred(coords[n])]

        smp(fo, "Parts_Structure", sorted(coords), eids)
        smp(fo, "ClampedEnds", nodes_where(lambda c: abs(c[2]) < TOL or abs(c[2] - L) < TOL))
        smp(fo, "SymmetryX", nodes_where(lambda c: abs(c[0]) < TOL))
        smp(fo, "SymmetryY", nodes_where(lambda c: abs(c[1]) < TOL))
        smp(fo, "StructureInterface", [nid(l, 0, k) for l in range(nz + 1) for k in range(nt + 1)],
            conditions=range(1, len(inner) + 1))
    return {"solid_nodes": len(coords), "solid_hexes": len(hexes), "wall_through": nw, "wall_around": nt,
            "wall_along": nz}


def smp(fo, name, nodes, elements=(), conditions=()):
    fo.write(f"Begin SubModelPart {name}\n    Begin SubModelPartNodes\n")
    for n in nodes:
        fo.write(f"        {n}\n")
    fo.write("    End SubModelPartNodes\n    Begin SubModelPartElements\n")
    for e in elements:
        fo.write(f"        {e}\n")
    fo.write("    End SubModelPartElements\n    Begin SubModelPartConditions\n")
    for c in conditions:
        fo.write(f"        {c}\n")
    fo.write("    End SubModelPartConditions\nEnd SubModelPart\n")


# --------------------------------------------------------------------------------------------- #
# Parameters
# --------------------------------------------------------------------------------------------- #

def vector_process(mp, var, constrained, value):
    return {"python_module": "assign_vector_variable_process", "kratos_module": "KratosMultiphysics",
            "process_name": "AssignVectorVariableProcess",
            "Parameters": {"model_part_name": mp, "variable_name": var, "constrained": constrained,
                           "value": value, "interval": [0.0, "End"]}}


def outlet_process(mp, backflow_velocity):
    """Pressure opening. backflow_velocity > 0 switches on Kratos's outlet inflow contribution
    (backflow stabilisation where u.n < 0) with that user-defined characteristic velocity; the switch
    and the velocity live in the shared ProcessInfo, so inlet and outlet must use the same value."""
    params = {"model_part_name": mp, "variable_name": "PRESSURE", "constrained": True,
              "value": 0.0, "hydrostatic_outlet": False, "h_top": 0.0}
    if backflow_velocity > 0.0:
        params.update({"outlet_inflow_contribution": True,
                       "outlet_inflow_contribution_characteristic_velocity_calculation": "user_defined",
                       "outlet_inflow_contribution_characteristic_velocity_value": backflow_velocity})
    return {"python_module": "apply_outlet_process", "kratos_module": "KratosMultiphysics.FluidDynamicsApplication",
            "process_name": "ApplyOutletProcess", "Parameters": params}


def parameters(args, vtk_steps):
    # The FSI solver requires the mesh-velocity scheme name to equal the fluid's ("bdf2",
    # BDF2TurbulentScheme). For --scheme euler, MainKratos_tube.py overwrites the fluid BDF
    # coefficients with BDF1 ones every step and replaces the mesh-velocity helper by BDF1.
    bdf = "bdf2"
    solid_scheme = "backward_euler" if args.scheme == "euler" else "bdf2"
    fluid = "FluidModelPart"
    p = {
        "analysis_stage": "KratosMultiphysics.FSIApplication.fsi_analysis",
        "problem_data": {"problem_name": "tube3d", "parallel_type": "OpenMP", "echo_level": 0,
                         "start_time": 0.0, "end_time": args.end_time},
        "solver_settings": {
            "solver_type": "partitioned", "coupling_scheme": "dirichlet_neumann", "echo_level": 1,
            "structure_solver_settings": {
                "solver_type": "Dynamic", "model_part_name": "Structure", "domain_size": 3,
                "echo_level": 0, "analysis_type": "linear", "time_integration_method": "implicit",
                "scheme_type": solid_scheme,
                "model_import_settings": {"input_type": "mdpa", "input_filename": SOLID_MDPA},
                "material_import_settings": {"materials_filename": "StructuralMaterials.json"},
                "time_stepping": {"time_step": args.dt},
                "rotation_dofs": False,
                **({"linear_solver_settings": {
                    # Kratos's default for the wall is the single-threaded direct sparse_lu, which at
                    # factor 2 (about 178k DOFs) takes minutes per coupling iteration; the dynamic wall
                    # matrix is symmetric positive definite and mass-dominated, so CG with AMG is cheap
                    "solver_type": "amgcl", "krylov_type": "cg", "smoother_type": "ilu0",
                    "coarsening_type": "aggregation", "use_block_matrices_if_possible": True,
                    "tolerance": args.structure_tol, "max_iteration": 1000, "verbosity": 0}}
                   if args.structure_solver == "amgcl" else {})},
            "fluid_solver_settings": {
                "solver_type": "Monolithic", "model_part_name": fluid, "domain_size": 3,
                "model_import_settings": {"input_type": "mdpa", "input_filename": FLUID_MDPA},
                "material_import_settings": {"materials_filename": "FluidMaterials.json"},
                "echo_level": args.fluid_echo, "compute_reactions": True, "maximum_iterations": 10,
                "relative_velocity_tolerance": args.fluid_rel_tol, "absolute_velocity_tolerance": 1e-10,
                "relative_pressure_tolerance": args.fluid_rel_tol, "absolute_pressure_tolerance": 1e-6,
                "volume_model_part_name": f"{fluid}.Parts_Fluid",
                "skin_parts": [f"{fluid}.Inlet", f"{fluid}.Outlet", f"{fluid}.FluidInterface"],
                "no_skin_parts": [f"{fluid}.SymmetryX", f"{fluid}.SymmetryY", f"{fluid}.MeshEnds"],
                "time_scheme": "bdf2",
                "time_stepping": {"automatic_time_step": False, "time_step": args.dt},
                # QSVMS with the scheme-level BDF2 (BDF2TurbulentScheme); TimeIntegratedQSVMS cannot
                # compute the right-hand side alone, which the FSI reactions need
                "formulation": {"element_type": "qsvms", "use_orthogonal_subscales": False,
                                "dynamic_tau": args.dynamic_tau, "element_manages_time_integration": False},
                "linear_solver_settings": {"solver_type": "amgcl", "tolerance": args.amgcl_tol,
                                           "max_iteration": 1000}},
            "mesh_solver_settings": {
                "echo_level": 0, "domain_size": 3, "model_part_name": fluid,
                "solver_type": "structural_similarity",
                "mesh_velocity_calculation": {"time_scheme": bdf}},
            "coupling_settings": {
                "nl_tol": args.nl_tol, "nl_max_it": 50, "solve_mesh_at_each_iteration": True,
                "coupling_strategy_settings": {"solver_type": args.accelerator},
                "structure_interfaces_list": ["Structure.StructureInterface"],
                "fluid_interfaces_list": [f"{fluid}.FluidInterface"],
                "mapper_settings": [{"mapper_face": "unique",
                                     "fluid_interface_submodelpart_name": f"{fluid}.FluidInterface",
                                     "structure_interface_submodelpart_name": "Structure.StructureInterface"}]}},
        "processes": {
            "structure_constraints_process_list": [
                vector_process("Structure.ClampedEnds", "DISPLACEMENT", [True, True, True], [0.0, 0.0, 0.0]),
                vector_process("Structure.SymmetryX", "DISPLACEMENT", [True, False, False], [0.0, None, None]),
                vector_process("Structure.SymmetryY", "DISPLACEMENT", [False, True, False], [None, 0.0, None])],
            "structure_loads_process_list": [],
            "fluid_initial_conditions_process_list": [],
            "fluid_boundary_conditions_process_list": [
                # The inlet pressure pulse is set every step by MainKratos_tube.py on top of this process
                outlet_process(f"{fluid}.Inlet", args.backflow_velocity),
                outlet_process(f"{fluid}.Outlet", args.backflow_velocity),
                vector_process(f"{fluid}.SymmetryX", "VELOCITY", [True, False, False], [0.0, None, None]),
                vector_process(f"{fluid}.SymmetryY", "VELOCITY", [False, True, False], [None, 0.0, None]),
                vector_process(f"{fluid}.MeshEnds", "MESH_DISPLACEMENT", [True, True, True], [0.0, 0.0, 0.0]),
                vector_process(f"{fluid}.SymmetryX", "MESH_DISPLACEMENT", [True, False, False], [0.0, None, None]),
                vector_process(f"{fluid}.SymmetryY", "MESH_DISPLACEMENT", [False, True, False], [None, 0.0, None])],
            "fluid_gravity": [],
            "fluid_auxiliar_process_list": []},
        "output_processes": {"vtk_output": [
            {"python_module": "vtk_output_process", "kratos_module": "KratosMultiphysics",
             "process_name": "VtkOutputProcess",
             "Parameters": {"model_part_name": mp, "output_control_type": "step", "output_interval": vtk_steps,
                            "file_format": "binary", "output_precision": 8, "output_sub_model_parts": False,
                            "output_path": "vtk_output", "save_output_files_in_folder": True,
                            "nodal_solution_step_data_variables": vars_}}
            for mp, vars_ in (("Structure.Parts_Structure", ["DISPLACEMENT", "VELOCITY"]),
                              (f"{fluid}.Parts_Fluid", ["VELOCITY", "PRESSURE", "MESH_DISPLACEMENT"]))]}}
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--factor", type=int, default=1)
    ap.add_argument("--dt", type=float, required=True)
    ap.add_argument("--scheme", choices=("euler", "bdf2"), required=True)
    ap.add_argument("--end-time", type=float, default=0.02)
    ap.add_argument("--nl-tol", type=float, default=1e-9)
    ap.add_argument("--vtk-every", type=float, default=1e-3)
    ap.add_argument("--fluid-rel-tol", type=float, default=1e-6)
    ap.add_argument("--amgcl-tol", type=float, default=1e-8)
    ap.add_argument("--fluid-echo", type=int, default=0)
    ap.add_argument("--backflow-velocity", type=float, default=0.0,
                    help="characteristic velocity (m/s) of the backflow stabilisation at the openings; 0 = off")
    ap.add_argument("--dynamic-tau", type=float, default=1.0,
                    help="weight of the 1/dt term in the QSVMS stabilisation parameter. 0 (the QSVMS "
                         "default) gives tau of about 1 ms here and damped the Euler L1 peak by about 7%%")
    ap.add_argument("--backflow-after", type=float, default=0.0,
                    help="switch the backflow stabilisation on only for t > this (s). Switching it on at the "
                         "end of the pulse (3.1e-3) made MVQN diverge, so the default is on from the start")
    ap.add_argument("--accelerator", choices=("MVQN", "MVQN_recursive", "MVQN_randomized_SVD"), default="MVQN",
                    help="MVQN keeps the full dense interface Jacobian (about 20k x 20k at factor 2, 18 GB and "
                         "serial updates: the factor-2 run did not finish its first step in 33 min); "
                         "MVQN_recursive rebuilds it from a buffer of the last 10 steps without storing it; "
                         "MVQN_randomized_SVD (low rank) segfaults reproducibly in step 3 in Kratos 10.3.0 "
                         "(null pointer at KratosFSIApplication.so +0x1bd73e), do not use")
    ap.add_argument("--structure-solver", choices=("direct", "amgcl"), default="direct",
                    help="direct = Kratos's default (sparse_lu, as in step 1); amgcl = CG with AMG")
    ap.add_argument("--structure-tol", type=float, default=1e-12,
                    help="relative residual tolerance of the amgcl wall solve")
    args = ap.parse_args()
    if os.path.exists(args.out):
        raise SystemExit(f"{args.out} exists, refusing to overwrite")
    os.makedirs(args.out)
    finfo = write_fluid(os.path.join(args.out, FLUID_MDPA + ".mdpa"), args.factor)
    sinfo = write_solid(os.path.join(args.out, SOLID_MDPA + ".mdpa"), args.factor,
                        node_offset=finfo["fluid_nodes"], elem_offset=finfo["fluid_tets"])
    vtk_steps = max(1, round(args.vtk_every / args.dt))
    with open(os.path.join(args.out, "ProjectParameters.json"), "w") as fo:
        json.dump(parameters(args, vtk_steps), fo, indent=4)
    with open(os.path.join(args.out, "FluidMaterials.json"), "w") as fo:
        json.dump({"properties": [{"model_part_name": "FluidModelPart.Parts_Fluid", "properties_id": 1,
                                   "Material": {"constitutive_law": {"name": "Newtonian3DLaw"},
                                                "Variables": {"DENSITY": RHO_F, "DYNAMIC_VISCOSITY": MU_F},
                                                "Tables": {}}}]}, fo, indent=4)
    with open(os.path.join(args.out, "StructuralMaterials.json"), "w") as fo:
        json.dump({"properties": [{"model_part_name": "Structure.Parts_Structure", "properties_id": 2,
                                   "Material": {"constitutive_law": {"name": "StructuralMechanicsApplication.LinearElastic3DLaw"},
                                                "Variables": {"DENSITY": RHO_S, "YOUNG_MODULUS": E_S,
                                                              "POISSON_RATIO": NU_S},
                                                "Tables": {}}}]}, fo, indent=4)
    case = dict(vars(args), **finfo, **sinfo, generator="tube_case.py",
                inlet={"pressure_Pa": P_IN, "hold_until_s": T_PULSE, "ramp_s": T_RAMP},
                backflow_after_s=args.backflow_after,
                stations_z_m=STATIONS, probe_offset_m=PROBE_OFFSET, inner_radius_m=R)
    with open(os.path.join(args.out, "case.json"), "w") as fo:
        json.dump(case, fo, indent=2)
    print(json.dumps({k: v for k, v in case.items() if k not in ("inlet",)}))


if __name__ == "__main__":
    main()
