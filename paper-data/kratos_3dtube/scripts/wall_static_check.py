#!/usr/bin/env python3
"""Static check of the Kratos wall model: the quarter tube wall of tube_case.py under a uniform internal
pressure, clamped at both ends, with the same symmetry constraints as the FSI case. The mid-length
radial displacement is compared with the Lame thick-cylinder solution with zero axial strain (the
clamped ends keep the mean axial strain at zero) and with free ends (zero axial stress).

Usage (from the tube3d case directory): python tools/wall_static_check.py <work dir> <factor> [...]
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tube_case as tc  # noqa: E402

P = tc.P_IN


def lame(r, plane_strain=True):
    a, b, E, nu = tc.R, tc.R + tc.H, tc.E_S, tc.NU_S
    if plane_strain:
        return (1 + nu) * P * a * a / (E * (b * b - a * a)) * ((1 - 2 * nu) * r + b * b / r)
    return P * a * a / (E * (b * b - a * a)) * ((1 - nu) * r + (1 + nu) * b * b / r)


def write_case(work, f):
    os.makedirs(work)
    path = os.path.join(work, tc.SOLID_MDPA + ".mdpa")
    info = tc.write_solid(path, f, 0, 0)
    text = open(path).read().replace("Begin Conditions SurfaceCondition3D3N",
                                     "Begin Conditions SurfaceLoadCondition3D3N")
    open(path, "w").write(text)
    with open(os.path.join(work, "StructuralMaterials.json"), "w") as fo:
        json.dump({"properties": [{"model_part_name": "Structure.Parts_Structure", "properties_id": 1,
                                   "Material": {"constitutive_law": {"name": "StructuralMechanicsApplication.LinearElastic3DLaw"},
                                                "Variables": {"DENSITY": tc.RHO_S, "YOUNG_MODULUS": tc.E_S,
                                                              "POISSON_RATIO": tc.NU_S}, "Tables": {}}}]}, fo)
    params = {
        "problem_data": {"problem_name": "wall", "parallel_type": "OpenMP", "echo_level": 0,
                         "start_time": 0.0, "end_time": 1.0},
        "solver_settings": {"solver_type": "Static", "model_part_name": "Structure", "domain_size": 3,
                            "echo_level": 0, "analysis_type": "linear",
                            "model_import_settings": {"input_type": "mdpa", "input_filename": tc.SOLID_MDPA},
                            "material_import_settings": {"materials_filename": "StructuralMaterials.json"},
                            "time_stepping": {"time_step": 1.0}, "rotation_dofs": False},
        "processes": {
            "constraints_process_list": [
                tc.vector_process("Structure.ClampedEnds", "DISPLACEMENT", [True, True, True], [0.0, 0.0, 0.0]),
                tc.vector_process("Structure.SymmetryX", "DISPLACEMENT", [True, False, False], [0.0, None, None]),
                tc.vector_process("Structure.SymmetryY", "DISPLACEMENT", [False, True, False], [None, 0.0, None])],
            "loads_process_list": [{
                "python_module": "assign_scalar_variable_to_conditions_process", "kratos_module": "KratosMultiphysics",
                "process_name": "AssignScalarVariableToConditionsProcess",
                "Parameters": {"model_part_name": "Structure.StructureInterface",
                               "variable_name": "POSITIVE_FACE_PRESSURE", "value": P, "interval": [0.0, "End"]}}]},
        "output_processes": {}}
    with open(os.path.join(work, "ProjectParameters.json"), "w") as fo:
        json.dump(params, fo, indent=2)
    return info


def run(work):
    import KratosMultiphysics as KM
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
    cwd = os.getcwd()
    os.chdir(work)
    with open("ProjectParameters.json") as f:
        params = KM.Parameters(f.read())
    model = KM.Model()
    StructuralMechanicsAnalysis(model, params).Run()
    mp = model.GetModelPart("Structure")
    out = {}
    for label, r in (("inner", tc.R), ("outer", tc.R + tc.H)):
        node = min(mp.Nodes, key=lambda n: n.X0**2 + (n.Y0 - r)**2 + (n.Z0 - tc.L / 2)**2)
        out[label] = node.GetSolutionStepValue(KM.DISPLACEMENT_Y)
    os.chdir(cwd)
    return out


def main():
    work_root, factors = sys.argv[1], [int(f) for f in sys.argv[2:]]
    print(f"Lame, zero axial strain: inner {lame(tc.R) * 1e3:.5f} mm, outer {lame(tc.R + tc.H) * 1e3:.5f} mm")
    print(f"Lame, zero axial stress: inner {lame(tc.R, False) * 1e3:.5f} mm, outer {lame(tc.R + tc.H, False) * 1e3:.5f} mm")
    for f in factors:
        work = os.path.join(work_root, f"wall_f{f}")
        info = write_case(work, f)
        u = run(work)
        print(f"factor {f} ({info['wall_through']} x {info['wall_around']} x {info['wall_along']} hexes): "
              f"inner {u['inner'] * 1e3:.5f} mm ({100 * (u['inner'] / lame(tc.R) - 1):+.2f}% vs plane strain), "
              f"outer {u['outer'] * 1e3:.5f} mm ({100 * (u['outer'] / lame(tc.R + tc.H) - 1):+.2f}%)")


if __name__ == "__main__":
    main()
