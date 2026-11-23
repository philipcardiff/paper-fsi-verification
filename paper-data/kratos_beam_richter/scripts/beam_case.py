#!/usr/bin/env python3
"""Write a Kratos input set for the beamInCrossFlow benchmark (Richter 2012, section 6.2) at any level.

Definition as solids4foam's verification campaign (branch verification/beam-richter,
tutorials/fluidSolidInteraction/beamInCrossFlow/verification/scripts/richter_case.py, 7e8cd7379):
- duct (0, 1.5) x (0, 0.4) x (-0.4, 0) m (half domain, symmetry plane z = 0), plate (0.4, 0.5) x
  (0, 0.2) x (-0.2, 0) clamped at y = 0;
- fluid rho 1000, nu 1e-3; plate St Venant-Kirchhoff, E 1.4 MPa, nu 0.4, rho 1000;
- inlet v_x = U y (0.4 - y) (0.4^2 - z^2) / (0.2^2 0.4^2), ramped by (1 - cos(pi t / 1 s)) / 2 over
  the first second (paraboloidInletVelocity); outlet p = 0 with natural velocity; no-slip walls;
- point A = (0.45, 0.15, -0.15), F = fluid force on the half plate.

Mesh: the graded block family of richter_case.py. Every block cell count is multiplied by --factor and
the counts through the plate thickness (x in 0.4-0.5, solid and the three fluid blocks beside and above
it) by --thickness-multiplier (2 = solids4foam's production family, 8 f solid cells through the
thickness). Per direction the blocks share their divisions and gradings, so the fluid is one graded
tensor-product grid with the plate cells removed: x 12f (ratio 1/8) | 4mf | 24f (ratio 8), y 8f | 6f
(ratio 6), z 6f (ratio 1/6) | 8f. Each fluid hexahedron becomes two prisms (the same xy triangulation
in every z layer) and three tetrahedra per prism by the minimum-index rule (conforming). The plate is
the solid block (4mf, 8f, 8f) of hexahedra; fluid and solid nodes coincide on the interface and A is
a node at every level.

Usage:
  python beam_case.py <out dir> --factor 1 --inlet-peak 0.2 [--thickness-multiplier 2] [--dt 0.05]
                      [--end-time 20] [--scheme euler|bdf2] [--dynamic-tau 0] [--nl-tol 1e-10]
                      [--structure-solver direct|amgcl] [--structure-tol 1e-12]
"""
import argparse
import json
import os

X0, X1, X2, X3 = 0.0, 0.4, 0.5, 1.5
Y0, Y1, Y2 = 0.0, 0.2, 0.4
Z0, Z1, Z2 = -0.4, -0.2, 0.0
RHO_F, MU_F = 1000.0, 1000.0 * 1e-3
RHO_S, E_S, NU_S = 1000.0, 1.4e6, 0.4
RAMP = 1.0
POINT_A = (0.45, 0.15, -0.15)
FLUID_MDPA, SOLID_MDPA = "beam_fluid", "beam_structure"
TOL = 1e-10


def fmt(v):
    return f"{v:.17g}"


def graded(a, b, n, ratio):
    """n + 1 node coordinates from a to b; ratio = last cell / first cell (OpenFOAM simpleGrading)."""
    if abs(ratio - 1.0) < 1e-14:
        sizes = [1.0] * n
    else:
        r = ratio ** (1.0 / (n - 1))
        sizes = [r ** k for k in range(n)]
    total = sum(sizes)
    xs, acc = [a], 0.0
    for s in sizes[:-1]:
        acc += s
        xs.append(a + (b - a) * acc / total)
    xs.append(b)
    return xs


def grid_lines(f, m):
    x = graded(X0, X1, 12 * f, 0.125) + graded(X1, X2, 4 * m * f, 1.0)[1:] + graded(X2, X3, 24 * f, 8.0)[1:]
    y = graded(Y0, Y1, 8 * f, 1.0) + graded(Y1, Y2, 6 * f, 6.0)[1:]
    z = graded(Z0, Z1, 6 * f, 1.0 / 6.0) + graded(Z1, Z2, 8 * f, 1.0)[1:]
    return x, y, z


def in_plate(xc, yc, zc):
    return X1 < xc < X2 and Y0 < yc < Y1 and Z1 < zc < Z2


def tet_volume(p):
    (x0, y0, z0), (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = p
    ax, ay, az = x1 - x0, y1 - y0, z1 - z0
    bx, by, bz = x2 - x0, y2 - y0, z2 - z0
    cx, cy, cz = x3 - x0, y3 - y0, z3 - z0
    return (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)) / 6.0


def prism_to_tets(v):
    """v = [a, b, c, a', b', c'] (bottom triangle, top triangle above it). Minimum-index rule
    (Dompierre et al. 1999), conforming across shared quadrilateral faces (as tube_case.py)."""
    perms = [(0, 1, 2, 3, 4, 5), (1, 2, 0, 4, 5, 3), (2, 0, 1, 5, 3, 4),
             (3, 5, 4, 0, 2, 1), (4, 3, 5, 1, 0, 2), (5, 4, 3, 2, 1, 0)]
    imin = v.index(min(v))
    perm = next(p for p in perms if p[0] == imin)
    w = [v[k] for k in perm]
    if min(w[1], w[5]) < min(w[2], w[4]):
        return [[w[0], w[1], w[2], w[5]], [w[0], w[1], w[5], w[4]], [w[0], w[4], w[5], w[3]]]
    return [[w[0], w[1], w[2], w[4]], [w[0], w[4], w[2], w[5]], [w[0], w[4], w[5], w[3]]]


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


def classify(p):
    """Boundary group of a fluid boundary face from its three node coordinates."""
    def on(axis, value):
        return all(abs(c[axis] - value) < TOL for c in p)
    if on(0, X0):
        return "inlet"
    if on(0, X3):
        return "outlet"
    if on(1, Y0):
        return "floor"
    if on(1, Y2):
        return "top"
    if on(2, Z0):
        return "side"
    if on(2, Z2):
        return "symmetry"
    if on(0, X1) or on(0, X2) or on(1, Y1) or on(2, Z1):
        return "interface"
    raise RuntimeError(f"unclassified boundary face {p}")


def write_fluid(path, faces_path, f, m):
    xs, ys, zs = grid_lines(f, m)
    nx, ny, nz = len(xs), len(ys), len(zs)

    def gid(i, j, k):
        return 1 + i + nx * (j + ny * k)

    tets = []
    hexes = 0
    for k in range(nz - 1):
        for j in range(ny - 1):
            for i in range(nx - 1):
                if in_plate(0.5 * (xs[i] + xs[i + 1]), 0.5 * (ys[j] + ys[j + 1]), 0.5 * (zs[k] + zs[k + 1])):
                    continue
                hexes += 1
                q = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]
                for tri in ((q[0], q[1], q[2]), (q[0], q[2], q[3])):
                    v = [gid(a, b, k) for a, b in tri] + [gid(a, b, k + 1) for a, b in tri]
                    tets += prism_to_tets(v)
    used = sorted({n for t in tets for n in t})
    new = {g: n for n, g in enumerate(used, 1)}

    def xyz(g):
        g -= 1
        return (xs[g % nx], ys[(g // nx) % ny], zs[g // (nx * ny)])

    coords = {new[g]: xyz(g) for g in used}
    tets = [[new[n] for n in t] for t in tets]
    for t in tets:
        if tet_volume([coords[n] for n in t]) < 0:
            t[1], t[2] = t[2], t[1]

    faces = {}
    for e, t in enumerate(tets, 1):
        for opp in range(4):
            face = [t[j] for j in range(4) if j != opp]
            key = tuple(sorted(face))
            faces[key] = None if key in faces else (face, t[opp], e)
    groups = {g: [] for g in ("inlet", "outlet", "floor", "top", "side", "symmetry", "interface")}
    for val in faces.values():
        if val is None:
            continue
        face, opp, e = val
        p = [coords[n] for n in face]
        if tet_volume(p + [coords[opp]]) > 0:   # opposite node on the positive side: flip to outward
            face = [face[0], face[2], face[1]]
        groups[classify(p)].append((face, e))

    nodes_of = {g: sorted({n for face, _ in v for n in face}) for g, v in groups.items()}
    interface = set(nodes_of["interface"])
    walls = sorted(set(nodes_of["floor"]) | set(nodes_of["top"]) | set(nodes_of["side"]))
    edge = sorted(n for n in interface if abs(coords[n][1] - Y0) < TOL)
    edge_set = set(edge)

    with open(path, "w") as fo:
        fo.write("Begin ModelPartData\nEnd ModelPartData\n\nBegin Properties 0\nEnd Properties\n\n")
        fo.write("Begin Properties 1\nEnd Properties\n\nBegin Nodes\n")
        for n in sorted(coords):
            fo.write(f"    {n} {fmt(coords[n][0])} {fmt(coords[n][1])} {fmt(coords[n][2])}\n")
        fo.write("End Nodes\n\nBegin Elements Element3D4N\n")
        for e, t in enumerate(tets, 1):
            fo.write(f"    {e} 1 {t[0]} {t[1]} {t[2]} {t[3]}\n")
        fo.write("End Elements\n\n")
        cid, cond_ids = 0, {}
        for g in ("inlet", "outlet", "interface"):
            fo.write("Begin Conditions SurfaceCondition3D3N\n")
            cond_ids[g] = []
            for face, _ in groups[g]:
                cid += 1
                cond_ids[g].append(cid)
                fo.write(f"    {cid} 0 {face[0]} {face[1]} {face[2]}\n")
            fo.write("End Conditions\n\n")
        smp(fo, "Parts_Fluid", sorted(coords), range(1, len(tets) + 1))
        smp(fo, "Inlet", nodes_of["inlet"], conditions=cond_ids["inlet"])
        smp(fo, "Outlet", nodes_of["outlet"], conditions=cond_ids["outlet"])
        smp(fo, "FluidInterface", nodes_of["interface"], conditions=cond_ids["interface"])
        smp(fo, "NoSlip", [n for n in walls if n not in interface])
        smp(fo, "Symmetry", nodes_of["symmetry"])
        fixed = set(walls) | set(nodes_of["inlet"]) | set(nodes_of["outlet"])
        smp(fo, "MeshFixed", sorted(n for n in fixed if n not in interface))

    # Faces for the force evaluation in MainKratos_beam.py: the interface faces and the floor faces
    # that touch the plate base edge, each with its fluid tetrahedron, oriented out of the fluid
    with open(faces_path, "w") as fo:
        json.dump({"interface": [face + [e] for face, e in groups["interface"]],
                   "floor_at_edge": [face + [e] for face, e in groups["floor"] if edge_set & set(face)],
                   "edge_nodes": edge}, fo)
    vol = sum(tet_volume([coords[n] for n in t]) for t in tets)
    return {"fluid_nodes": len(coords), "fluid_tets": len(tets), "fluid_hex_equivalent": hexes,
            "fluid_volume": vol, "boundary_faces": {g: len(v) for g, v in groups.items()},
            "grid_lines": [nx, ny, nz], "min_dx": min(b - a for a, b in zip(xs, xs[1:])),
            "plate_base_edge_nodes": len(edge)}


def hex_jacobian_sign(p):
    def avg(ids):
        return [sum(p[i][d] for i in ids) / 4 for d in range(3)]
    dx = [a - b for a, b in zip(avg((1, 2, 5, 6)), avg((0, 3, 4, 7)))]
    dy = [a - b for a, b in zip(avg((2, 3, 6, 7)), avg((0, 1, 4, 5)))]
    dz = [a - b for a, b in zip(avg((4, 5, 6, 7)), avg((0, 1, 2, 3)))]
    return (dx[0] * (dy[1] * dz[2] - dy[2] * dz[1]) - dx[1] * (dy[0] * dz[2] - dy[2] * dz[0])
            + dx[2] * (dy[0] * dz[1] - dy[1] * dz[0]))


def solid_grid(f, m):
    xs, ys, zs = grid_lines(f, m)
    return ([x for x in xs if X1 - TOL <= x <= X2 + TOL], [y for y in ys if y <= Y1 + TOL],
            [z for z in zs if z >= Z1 - TOL])


def write_solid(path, f, m, node_offset, elem_offset, element="TotalLagrangianElement3D8N",
                condition="SurfaceCondition3D3N", interface_faces=("x0", "x1", "y1", "z0")):
    xs, ys, zs = solid_grid(f, m)
    nx, ny, nz = len(xs), len(ys), len(zs)

    def nid(i, j, k):
        return node_offset + 1 + i + nx * (j + ny * k)

    coords = {nid(i, j, k): (xs[i], ys[j], zs[k]) for k in range(nz) for j in range(ny) for i in range(nx)}
    hexes = []
    for k in range(nz - 1):
        for j in range(ny - 1):
            for i in range(nx - 1):
                h = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                     nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                if hex_jacobian_sign([coords[n] for n in h]) < 0:
                    h = [h[0], h[3], h[2], h[1], h[4], h[7], h[6], h[5]]
                hexes.append(h)
    # Interface faces as triangles (the 3D FSI coupling skin needs triangles):
    # x0 = upstream face x = 0.4, x1 = downstream x = 0.5, y1 = top y = 0.2, z0 = side z = -0.2
    quads = {"x0": [], "x1": [], "y1": [], "z0": []}
    for k in range(nz - 1):
        for j in range(ny - 1):
            quads["x0"].append([nid(0, j, k), nid(0, j, k + 1), nid(0, j + 1, k + 1), nid(0, j + 1, k)])
            quads["x1"].append([nid(nx - 1, j, k), nid(nx - 1, j + 1, k), nid(nx - 1, j + 1, k + 1),
                                nid(nx - 1, j, k + 1)])
    for k in range(nz - 1):
        for i in range(nx - 1):
            quads["y1"].append([nid(i, ny - 1, k), nid(i, ny - 1, k + 1), nid(i + 1, ny - 1, k + 1),
                                nid(i + 1, ny - 1, k)])
    for j in range(ny - 1):
        for i in range(nx - 1):
            quads["z0"].append([nid(i, j, 0), nid(i, j + 1, 0), nid(i + 1, j + 1, 0), nid(i + 1, j, 0)])
    tris = []
    for name in interface_faces:
        for q in quads[name]:
            tris += [[q[0], q[1], q[2]], [q[0], q[2], q[3]]]
    with open(path, "w") as fo:
        fo.write("Begin ModelPartData\nEnd ModelPartData\n\nBegin Properties 0\nEnd Properties\n\nBegin Nodes\n")
        for n in sorted(coords):
            fo.write(f"    {n} {fmt(coords[n][0])} {fmt(coords[n][1])} {fmt(coords[n][2])}\n")
        fo.write(f"End Nodes\n\nBegin Elements {element}\n")
        eids = []
        for e, h in enumerate(hexes, elem_offset + 1):
            eids.append(e)
            fo.write(f"    {e} 0 " + " ".join(str(n) for n in h) + "\n")
        fo.write(f"End Elements\n\nBegin Conditions {condition}\n")
        for c, t in enumerate(tris, 1):
            fo.write(f"    {c} 0 {t[0]} {t[1]} {t[2]}\n")
        fo.write("End Conditions\n\n")

        def nodes_where(pred):
            return [n for n in sorted(coords) if pred(coords[n])]

        smp(fo, "Parts_Structure", sorted(coords), eids)
        smp(fo, "Clamped", nodes_where(lambda c: abs(c[1] - Y0) < TOL))
        smp(fo, "Symmetry", nodes_where(lambda c: abs(c[2] - Z2) < TOL))
        smp(fo, "StructureInterface", sorted({n for t in tris for n in t}), conditions=range(1, len(tris) + 1))
    return {"solid_nodes": len(coords), "solid_hexes": len(hexes), "solid_cells_xyz": [nx - 1, ny - 1, nz - 1]}


def vector_process(mp, var, constrained, value):
    return {"python_module": "assign_vector_variable_process", "kratos_module": "KratosMultiphysics",
            "process_name": "AssignVectorVariableProcess",
            "Parameters": {"model_part_name": mp, "variable_name": var, "constrained": constrained,
                           "value": value, "interval": [0.0, "End"]}}


def parameters(args, vtk_steps):
    # As tube_case.py: the FSI solver requires the mesh-velocity scheme name to equal the fluid's
    # ("bdf2", BDF2TurbulentScheme); for --scheme euler MainKratos_beam.py sets BDF1 coefficients in
    # the fluid every step and a BDF1 mesh-velocity helper
    solid_scheme = "backward_euler" if args.scheme == "euler" else "bdf2"
    fluid = "FluidModelPart"
    return {
        "analysis_stage": "KratosMultiphysics.FSIApplication.fsi_analysis",
        "problem_data": {"problem_name": "beam", "parallel_type": "OpenMP", "echo_level": 0,
                         "start_time": 0.0, "end_time": args.end_time},
        "solver_settings": {
            "solver_type": "partitioned", "coupling_scheme": "dirichlet_neumann", "echo_level": 1,
            "structure_solver_settings": {
                "solver_type": "Dynamic", "model_part_name": "Structure", "domain_size": 3,
                "echo_level": 0, "analysis_type": "non_linear", "time_integration_method": "implicit",
                "scheme_type": solid_scheme,
                "model_import_settings": {"input_type": "mdpa", "input_filename": SOLID_MDPA},
                "material_import_settings": {"materials_filename": "StructuralMaterials.json"},
                "time_stepping": {"time_step": args.dt},
                "convergence_criterion": "residual_criterion",
                "residual_relative_tolerance": 1e-9, "residual_absolute_tolerance": 1e-11,
                "max_iteration": 20, "rotation_dofs": False,
                **({"linear_solver_settings": {
                    # Kratos's default for the plate is the direct sparse_lu, single-threaded and slow at
                    # factor 4 (about 108k DOFs); the dynamic plate matrix is symmetric positive definite,
                    # so CG with AMG is cheap (as for the 3dTube wall)
                    "solver_type": "amgcl", "krylov_type": "cg", "smoother_type": "ilu0",
                    "coarsening_type": "aggregation", "use_block_matrices_if_possible": True,
                    "tolerance": args.structure_tol, "max_iteration": 1000, "verbosity": 0}}
                   if args.structure_solver == "amgcl" else {})},
            "fluid_solver_settings": {
                "solver_type": "Monolithic", "model_part_name": fluid, "domain_size": 3,
                "model_import_settings": {"input_type": "mdpa", "input_filename": FLUID_MDPA},
                "material_import_settings": {"materials_filename": "FluidMaterials.json"},
                "echo_level": args.fluid_echo, "compute_reactions": True, "maximum_iterations": 20,
                "relative_velocity_tolerance": args.fluid_rel_tol, "absolute_velocity_tolerance": 1e-12,
                "relative_pressure_tolerance": args.fluid_rel_tol, "absolute_pressure_tolerance": 1e-9,
                "volume_model_part_name": f"{fluid}.Parts_Fluid",
                "skin_parts": [f"{fluid}.Inlet", f"{fluid}.Outlet", f"{fluid}.FluidInterface"],
                "no_skin_parts": [f"{fluid}.NoSlip", f"{fluid}.Symmetry", f"{fluid}.MeshFixed"],
                "time_scheme": "bdf2",
                "time_stepping": {"automatic_time_step": False, "time_step": args.dt},
                "formulation": {"element_type": "qsvms", "use_orthogonal_subscales": False,
                                "dynamic_tau": args.dynamic_tau, "element_manages_time_integration": False},
                "linear_solver_settings": {"solver_type": "amgcl", "tolerance": args.amgcl_tol,
                                           "max_iteration": 1000}},
            "mesh_solver_settings": {
                "echo_level": 0, "domain_size": 3, "model_part_name": fluid,
                "solver_type": "structural_similarity",
                "mesh_velocity_calculation": {"time_scheme": "bdf2"}},
            "coupling_settings": {
                "nl_tol": args.nl_tol, "nl_max_it": 50, "solve_mesh_at_each_iteration": True,
                "coupling_strategy_settings": {"solver_type": "MVQN"},
                "structure_interfaces_list": ["Structure.StructureInterface"],
                "fluid_interfaces_list": [f"{fluid}.FluidInterface"],
                "mapper_settings": [{"mapper_face": "unique",
                                     "fluid_interface_submodelpart_name": f"{fluid}.FluidInterface",
                                     "structure_interface_submodelpart_name": "Structure.StructureInterface"}]}},
        "processes": {
            "structure_constraints_process_list": [
                vector_process("Structure.Clamped", "DISPLACEMENT", [True, True, True], [0.0, 0.0, 0.0]),
                vector_process("Structure.Symmetry", "DISPLACEMENT", [False, False, True], [None, None, 0.0])],
            "structure_loads_process_list": [],
            "fluid_initial_conditions_process_list": [],
            "fluid_boundary_conditions_process_list": [
                # The inlet velocity is set every step by MainKratos_beam.py
                {"python_module": "apply_outlet_process", "kratos_module": "KratosMultiphysics.FluidDynamicsApplication",
                 "process_name": "ApplyOutletProcess",
                 "Parameters": {"model_part_name": f"{fluid}.Outlet", "variable_name": "PRESSURE",
                                "constrained": True, "value": 0.0, "hydrostatic_outlet": False, "h_top": 0.0}},
                vector_process(f"{fluid}.NoSlip", "VELOCITY", [True, True, True], [0.0, 0.0, 0.0]),
                vector_process(f"{fluid}.Symmetry", "VELOCITY", [False, False, True], [None, None, 0.0]),
                vector_process(f"{fluid}.MeshFixed", "MESH_DISPLACEMENT", [True, True, True], [0.0, 0.0, 0.0]),
                vector_process(f"{fluid}.Symmetry", "MESH_DISPLACEMENT", [False, False, True], [None, None, 0.0])],
            "fluid_gravity": [],
            "fluid_auxiliar_process_list": []},
        "output_processes": {"vtk_output": [
            {"python_module": "vtk_output_process", "kratos_module": "KratosMultiphysics",
             "process_name": "VtkOutputProcess",
             "Parameters": {"model_part_name": mp, "output_control_type": "step", "output_interval": vtk_steps,
                            "file_format": "binary", "output_precision": 8, "output_sub_model_parts": False,
                            "output_path": "vtk_output", "save_output_files_in_folder": True,
                            "nodal_solution_step_data_variables": vars_}}
            for mp, vars_ in (("Structure.Parts_Structure", ["DISPLACEMENT"]),
                              (f"{fluid}.Parts_Fluid", ["VELOCITY", "PRESSURE", "MESH_DISPLACEMENT", "REACTION"]))]}}


def write_materials(out, solid_law="StructuralMechanicsApplication.LinearElastic3DLaw"):
    with open(os.path.join(out, "FluidMaterials.json"), "w") as fo:
        json.dump({"properties": [{"model_part_name": "FluidModelPart.Parts_Fluid", "properties_id": 1,
                                   "Material": {"constitutive_law": {"name": "Newtonian3DLaw"},
                                                "Variables": {"DENSITY": RHO_F, "DYNAMIC_VISCOSITY": MU_F},
                                                "Tables": {}}}]}, fo, indent=4)
    # TotalLagrangianElement3D8N with LinearElastic3DLaw: S = C : E (Green-Lagrange), St Venant-Kirchhoff
    with open(os.path.join(out, "StructuralMaterials.json"), "w") as fo:
        json.dump({"properties": [{"model_part_name": "Structure.Parts_Structure", "properties_id": 2,
                                   "Material": {"constitutive_law": {"name": solid_law},
                                                "Variables": {"DENSITY": RHO_S, "YOUNG_MODULUS": E_S,
                                                              "POISSON_RATIO": NU_S},
                                                "Tables": {}}}]}, fo, indent=4)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--factor", type=int, default=1, help="1, 2, 4 = solids4foam production L0, L1, L2")
    ap.add_argument("--thickness-multiplier", type=int, default=2)
    ap.add_argument("--inlet-peak", type=float, required=True, help="U (m/s): 0.3 printed, 0.2 Richter's values")
    ap.add_argument("--dt", type=float, default=0.05)
    ap.add_argument("--end-time", type=float, default=20.0)
    ap.add_argument("--scheme", choices=("euler", "bdf2"), default="euler")
    ap.add_argument("--dynamic-tau", type=float, default=0.0,
                    help="0: tau independent of dt, so the steady state does not depend on dt")
    ap.add_argument("--nl-tol", type=float, default=1e-10)
    ap.add_argument("--fluid-rel-tol", type=float, default=1e-8)
    ap.add_argument("--amgcl-tol", type=float, default=1e-10)
    ap.add_argument("--fluid-echo", type=int, default=0)
    ap.add_argument("--vtk-every", type=float, default=5.0)
    ap.add_argument("--steady-window", type=float, default=1.0)
    ap.add_argument("--steady-tol", type=float, default=1e-6,
                    help="stop when u_x(A) and F_x change by less than this (relative) over the window")
    ap.add_argument("--steady-after", type=float, default=4.0)
    ap.add_argument("--structure-solver", choices=("direct", "amgcl"), default="direct",
                    help="direct = Kratos's default (sparse_lu, levels 0 and 1); amgcl = CG with AMG")
    ap.add_argument("--structure-tol", type=float, default=1e-12,
                    help="relative residual tolerance of the amgcl plate solve")
    args = ap.parse_args()
    if os.path.exists(args.out):
        raise SystemExit(f"{args.out} exists, refusing to overwrite")
    os.makedirs(args.out)
    finfo = write_fluid(os.path.join(args.out, FLUID_MDPA + ".mdpa"), os.path.join(args.out, "force_faces.json"),
                        args.factor, args.thickness_multiplier)
    sinfo = write_solid(os.path.join(args.out, SOLID_MDPA + ".mdpa"), args.factor, args.thickness_multiplier,
                        node_offset=finfo["fluid_nodes"], elem_offset=finfo["fluid_tets"])
    vtk_steps = max(1, round(args.vtk_every / args.dt))
    with open(os.path.join(args.out, "ProjectParameters.json"), "w") as fo:
        json.dump(parameters(args, vtk_steps), fo, indent=4)
    write_materials(args.out)
    opts = vars(args).copy()
    if args.structure_solver == "direct":   # keep case.json of the direct runs as it was
        del opts["structure_solver"], opts["structure_tol"]
    case = dict(opts, **finfo, **sinfo, generator="beam_case.py", ramp_s=RAMP, point_a=POINT_A,
                inlet={"y_max": Y2, "z_half": -Z0}, mu=MU_F)
    with open(os.path.join(args.out, "case.json"), "w") as fo:
        json.dump(case, fo, indent=2)
    print(json.dumps(case))


if __name__ == "__main__":
    main()
