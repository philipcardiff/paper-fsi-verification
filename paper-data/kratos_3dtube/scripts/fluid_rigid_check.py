#!/usr/bin/env python3
"""Fluid-only check of the pressure openings: the tube_case.py fluid mesh as a rigid quarter tube with
a constant inlet pressure 1333.2 Pa and p = 0 at the outlet, applied exactly as in the FSI case
(ApplyOutletProcess, OUTLET flags, EXTERNAL_PRESSURE). At early times the core accelerates uniformly at
dp/(rho L) = 26.66 m/s^2 (the viscous boundary layer, sqrt(nu t) = 0.05 mm at 1 ms, is thin). The axis
velocity at mid-length is compared with that.

Usage (from the tube3d case directory): python tools/fluid_rigid_check.py <work dir> [--scheme euler|bdf2]
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tube_case as tc  # noqa: E402

DT, STEPS = 1e-4, 10


def main():
    work = sys.argv[1]
    scheme = sys.argv[sys.argv.index("--scheme") + 1] if "--scheme" in sys.argv else "bdf2"
    os.makedirs(work)
    tc.write_fluid(os.path.join(work, tc.FLUID_MDPA + ".mdpa"), 1)
    with open(os.path.join(work, "FluidMaterials.json"), "w") as fo:
        json.dump({"properties": [{"model_part_name": "FluidModelPart.Parts_Fluid", "properties_id": 1,
                                   "Material": {"constitutive_law": {"name": "Newtonian3DLaw"},
                                                "Variables": {"DENSITY": tc.RHO_F, "DYNAMIC_VISCOSITY": tc.MU_F},
                                                "Tables": {}}}]}, fo)
    fl = "FluidModelPart"
    params = {
        "analysis_stage": "KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis",
        "problem_data": {"problem_name": "rigid", "parallel_type": "OpenMP", "echo_level": 0,
                         "start_time": 0.0, "end_time": DT * STEPS},
        "solver_settings": {
            "solver_type": "Monolithic", "model_part_name": fl, "domain_size": 3,
            "model_import_settings": {"input_type": "mdpa", "input_filename": tc.FLUID_MDPA},
            "material_import_settings": {"materials_filename": "FluidMaterials.json"},
            "echo_level": 0, "compute_reactions": True, "maximum_iterations": 10,
            "relative_velocity_tolerance": 1e-8, "absolute_velocity_tolerance": 1e-12,
            "relative_pressure_tolerance": 1e-8, "absolute_pressure_tolerance": 1e-8,
            "volume_model_part_name": f"{fl}.Parts_Fluid",
            "skin_parts": [f"{fl}.Inlet", f"{fl}.Outlet", f"{fl}.FluidInterface"],
            "no_skin_parts": [f"{fl}.SymmetryX", f"{fl}.SymmetryY"],
            "time_scheme": "bdf2", "time_stepping": {"automatic_time_step": False, "time_step": DT},
            "formulation": {"element_type": "qsvms", "use_orthogonal_subscales": False, "dynamic_tau": 0.0,
                            "element_manages_time_integration": False},
            "linear_solver_settings": {"solver_type": "amgcl", "tolerance": 1e-10, "max_iteration": 1000}},
        "processes": {
            "initial_conditions_process_list": [],
            "boundary_conditions_process_list": [
                dict(tc.outlet_process(f"{fl}.Inlet", 0.0)),
                tc.outlet_process(f"{fl}.Outlet", 0.0),
                tc.vector_process(f"{fl}.FluidInterface", "VELOCITY", [True, True, True], [0.0, 0.0, 0.0]),
                tc.vector_process(f"{fl}.SymmetryX", "VELOCITY", [True, False, False], [0.0, None, None]),
                tc.vector_process(f"{fl}.SymmetryY", "VELOCITY", [False, True, False], [None, 0.0, None])],
            "gravity": [], "auxiliar_process_list": []},
        "output_processes": {}}
    params["processes"]["boundary_conditions_process_list"][0]["Parameters"]["value"] = tc.P_IN
    with open(os.path.join(work, "ProjectParameters.json"), "w") as fo:
        json.dump(params, fo, indent=2)

    import KratosMultiphysics as KM
    from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis
    os.chdir(work)
    with open("ProjectParameters.json") as f:
        kp = KM.Parameters(f.read())
    model = KM.Model()

    class Check(FluidDynamicsAnalysis):
        def InitializeSolutionStep(self):
            super().InitializeSolutionStep()
            if scheme == "euler":
                pi = self.model.GetModelPart(fl).ProcessInfo
                dt = pi[KM.DELTA_TIME]
                pi[KM.BDF_COEFFICIENTS] = KM.Vector([1.0 / dt, -1.0 / dt, 0.0])

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            mp = self.model.GetModelPart(fl)
            axis = min(mp.Nodes, key=lambda n: n.X0**2 + n.Y0**2 + (n.Z0 - tc.L / 2)**2)
            inlet = self.model.GetModelPart(f"{fl}.Inlet")
            p_in = [n.GetSolutionStepValue(KM.PRESSURE) for n in inlet.Nodes]
            pe_in = [n.GetSolutionStepValue(KM.EXTERNAL_PRESSURE) for n in inlet.Nodes]
            n_outlet_flag = sum(1 for c in inlet.Conditions if c.Is(KM.OUTLET))
            u = axis.GetSolutionStepValue(KM.VELOCITY_Z)
            exact = tc.P_IN / (tc.RHO_F * tc.L) * self.time
            print(f"CHECK t={self.time:.1e} s  u_axis(L/2)={u:.6e} m/s  uniform-acceleration {exact:.6e} "
                  f"({100 * (u / exact - 1):+.2f}%)  inlet p [{min(p_in):.1f}, {max(p_in):.1f}] "
                  f"ext p [{min(pe_in):.1f}, {max(pe_in):.1f}]  inlet conditions with OUTLET flag "
                  f"{n_outlet_flag}/{inlet.NumberOfConditions()}  BDF {list(self.model.GetModelPart(fl).ProcessInfo[KM.BDF_COEFFICIENTS])}",
                  flush=True)

    Check(model, kp).Run()


if __name__ == "__main__":
    main()
