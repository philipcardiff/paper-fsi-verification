"""Run the Kratos 3dTube case written by tube_case.py and record a per-step history.

Built on the MainKratos.py of the cavity case. Additions for the tube:
  - the inlet pressure pulse of the solids4foam tutorial (0/fluid/time-series): 1333.2 Pa held to
    3 ms, linear to 0 at 3.1 ms, 0 after; PRESSURE is fixed and EXTERNAL_PRESSURE set on the inlet
    nodes every step, after the processes (ApplyOutletProcess sets the OUTLET flags);
  - scheme "euler": the fluid uses QSVMS with BDF2TurbulentScheme, which writes its three BDF
    coefficients to the ProcessInfo in InitializeSolutionStep and reads them from there afterwards;
    they are overwritten with BDF1 ones [1/dt, -1/dt, 0] right after every step's initialisation.
    The mesh solver's velocity helper is replaced by TimeDiscretization.BDF1() (the parameters say
    "bdf2" for both, because the FSI solver requires the two names to match); the wall uses its
    backward-Euler scheme;
  - history.csv: time, coupling solves, converged flag, RMS interface residual at acceptance, wall
    time, inlet pressure, then for each station z the wall displacement at (0, R, z) (u_y radial,
    u_z axial) and the pressure at (offset, offset, z), interpolated in the tetrahedron that contains
    the probe in the initial mesh;
  - run_info.json: Kratos version, threads, host, wall time.

Usage: python MainKratos_tube.py [ProjectParameters.json] [--end-time T]
"""
import importlib
import json
import math
import os
import platform
import socket
import sys
import time

import KratosMultiphysics as KM


def load_parameters(argv):
    path, end_time, args = "ProjectParameters.json", None, list(argv)
    if "--end-time" in args:
        i = args.index("--end-time")
        end_time = float(args[i + 1])
        del args[i:i + 2]
    if args:
        path = args[0]
    with open(path) as f:
        parameters = KM.Parameters(f.read())
    if end_time is not None:
        parameters["problem_data"]["end_time"].SetDouble(end_time)
    return path, parameters


def stage_class(parameters):
    module_name = parameters["analysis_stage"].GetString()
    class_name = "".join(x.title() for x in module_name.split(".")[-1].split("_"))
    return getattr(importlib.import_module(module_name), class_name)


class BDF1Coefficients:
    """BDF1 values for the three coefficients BDF2TurbulentScheme reads from the ProcessInfo."""

    def ComputeAndSaveBDFCoefficients(self, process_info):
        dt = process_info[KM.DELTA_TIME]
        process_info[KM.BDF_COEFFICIENTS] = KM.Vector([1.0 / dt, -1.0 / dt, 0.0])


def inlet_pressure(t, inlet):
    p, hold, ramp = inlet["pressure_Pa"], inlet["hold_until_s"], inlet["ramp_s"]
    if t <= hold + 1e-15:
        return p
    if t < hold + ramp:
        return p * (hold + ramp - t) / ramp
    return 0.0


def barycentric(p, a, b, c, d):
    def det3(u, v, w):
        return (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0]))
    def sub(x, y):
        return [x[i] - y[i] for i in range(3)]
    vol = det3(sub(b, a), sub(c, a), sub(d, a))
    l1 = det3(sub(p, a), sub(c, a), sub(d, a)) / vol
    l2 = det3(sub(b, a), sub(p, a), sub(d, a)) / vol
    l3 = det3(sub(b, a), sub(c, a), sub(p, a)) / vol
    return [1.0 - l1 - l2 - l3, l1, l2, l3]


def create_stage(cls, model, parameters, case):

    class RecordedTubeAnalysis(cls):

        def Initialize(self):
            super().Initialize()
            solver = self._GetSolver()
            if case["scheme"] == "euler":
                scheme_name = type(solver.fluid_solver._GetScheme()).__name__
                if scheme_name != "BDF2TurbulentScheme":
                    raise RuntimeError(f"expected BDF2TurbulentScheme in the fluid, found {scheme_name}")
                mesh = solver.mesh_solver
                if not hasattr(mesh, "time_int_helper"):
                    raise RuntimeError("mesh solver has no time_int_helper to replace")
                mesh.time_int_helper = KM.TimeDiscretization.BDF1()

            structure = self.model.GetModelPart("Structure")
            radius = case["inner_radius_m"]
            self.stations = []
            for z in case["stations_z_m"]:
                node = min(structure.Nodes, key=lambda n: n.X0**2 + (n.Y0 - radius)**2 + (n.Z0 - z)**2)
                if math.sqrt(node.X0**2 + (node.Y0 - radius)**2 + (node.Z0 - z)**2) > 1e-10:
                    raise RuntimeError(f"no wall node at (0, {radius}, {z})")
                self.stations.append(node)

            fluid_mp = self.model.GetModelPart("FluidModelPart.Parts_Fluid")
            off = case["probe_offset_m"]
            self.probes = []
            for z in case["stations_z_m"]:
                p = [off, off, z]
                found = None
                for elem in fluid_mp.Elements:
                    nodes = elem.GetNodes()
                    xyz = [[n.X0, n.Y0, n.Z0] for n in nodes]
                    if (min(c[2] for c in xyz) - 1e-12 > z or max(c[2] for c in xyz) + 1e-12 < z
                            or min(c[0] for c in xyz) - 1e-12 > off or min(c[1] for c in xyz) - 1e-12 > off):
                        continue
                    w = barycentric(p, *xyz)
                    if min(w) >= -1e-10:
                        found = (list(nodes), w)
                        break
                if found is None:
                    raise RuntimeError(f"pressure probe {p} not found in the fluid mesh")
                self.probes.append(found)

            self.inlet = self.model.GetModelPart("FluidModelPart.Inlet")

            self.step_rec = {}
            check = solver._CheckFSIConvergence
            def recorded_check(nl_it, residual_norm):
                ok = check(nl_it, residual_norm)
                self.step_rec["checked_it"] = nl_it
                return ok
            solver._CheckFSIConvergence = recorded_check
            interface = solver._GetFSICouplingInterfaceStructure()
            compute = interface.ComputeResidualVector
            def recorded_compute(*a, **k):
                norm = compute(*a, **k)
                self.step_rec["solves"] = self.step_rec.get("solves", 0) + 1
                self.step_rec["last_norm"] = norm
                return norm
            interface.ComputeResidualVector = recorded_compute
            solve = solver.SolveSolutionStep
            def recorded_solve():
                ok = solve()
                self.step_rec["converged"] = bool(ok)
                return ok
            solver.SolveSolutionStep = recorded_solve
            self.n_dofs = solver._GetPartitionedFSIUtilities().GetInterfaceResidualSize(
                solver._GetStructureInterfaceSubmodelPart())

            self.history = open("history.csv", "w")
            cols = ["time", "coupling_solves", "converged", "rms_residual_at_acceptance", "wall_step_s",
                    "p_inlet"]
            for z in case["stations_z_m"]:
                tag = f"z{round(z * 1e4):03d}"
                cols += [f"ur_{tag}", f"uz_{tag}", f"p_{tag}"]
            self.history.write(",".join(cols) + "\n")
            self.t_wall0 = time.time()
            self.t_step0 = self.t_wall0
            sys.stdout.flush()

        def ApplyBoundaryConditions(self):
            super().ApplyBoundaryConditions()
            self.p_in = inlet_pressure(self.time, case["inlet"])
            for node in self.inlet.Nodes:
                node.SetSolutionStepValue(KM.PRESSURE, self.p_in)
                node.Fix(KM.PRESSURE)
                node.SetSolutionStepValue(KM.EXTERNAL_PRESSURE, self.p_in)

        def InitializeSolutionStep(self):
            self.step_rec = {}
            self.t_step0 = time.time()
            super().InitializeSolutionStep()
            if case["scheme"] == "euler":
                BDF1Coefficients().ComputeAndSaveBDFCoefficients(
                    self.model.GetModelPart("FluidModelPart").ProcessInfo)

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            r = self.step_rec
            rms = r.get("last_norm", float("nan")) / math.sqrt(self.n_dofs)
            vals = [f"{self.time:.10g}", str(r.get("solves", 0)), str(int(r.get("converged", False))),
                    f"{rms:.6e}", f"{time.time() - self.t_step0:.4f}", f"{self.p_in:.6g}"]
            for node, (pnodes, w) in zip(self.stations, self.probes):
                p = sum(wi * n.GetSolutionStepValue(KM.PRESSURE) for wi, n in zip(w, pnodes))
                vals += [f"{node.GetSolutionStepValue(KM.DISPLACEMENT_Y):.12e}",
                         f"{node.GetSolutionStepValue(KM.DISPLACEMENT_Z):.12e}", f"{p:.8e}"]
            self.history.write(",".join(vals) + "\n")
            self.history.flush()
            sys.stdout.flush()

        def Finalize(self):
            super().Finalize()
            self.history.close()
            self.wall_total = time.time() - self.t_wall0

    return RecordedTubeAnalysis(model, parameters)


def kratos_version():
    for getter in (lambda: KM.Kernel.Version(), lambda: KM.KratosGlobals.Kernel.Version()):
        try:
            return getter()
        except Exception:
            pass
    return "unknown"


if __name__ == "__main__":
    params_path, parameters = load_parameters(sys.argv[1:])
    with open(os.path.join(os.path.dirname(os.path.abspath(params_path)), "case.json")) as f:
        case = json.load(f)
    model = KM.Model()
    simulation = create_stage(stage_class(parameters), model, parameters, case)
    t0 = time.time()
    simulation.Run()
    info = {
        "kratos_version": kratos_version(), "python": platform.python_version(),
        "host": socket.gethostname(),
        "cpu": next((l.split(":", 1)[1].strip() for l in open("/proc/cpuinfo")
                     if l.startswith("model name")), "unknown"),
        "omp_num_threads_env": os.environ.get("OMP_NUM_THREADS"),
        "kratos_threads": KM.ParallelUtilities.GetNumThreads(),
        "parameters_file": params_path, "scheme": case["scheme"],
        "end_time": parameters["problem_data"]["end_time"].GetDouble(),
        "wall_total_s": round(time.time() - t0, 3),
        "wall_solution_loop_s": round(getattr(simulation, "wall_total", float("nan")), 3),
    }
    with open("run_info.json", "w") as f:
        json.dump(info, f, indent=2)
    print(json.dumps(info, indent=2))
