"""Run the Kratos beamInCrossFlow case written by beam_case.py and record a per-step history.

Built on MainKratos_tube.py. For the beam:
  - inlet: v_x = U r(t) y (0.4 - y) (0.4^2 - z^2) / (0.2^2 0.4^2), r(t) = (1 - cos(pi t / 1 s)) / 2 for
    t < 1 s, 1 after (solids4foam paraboloidInletVelocity), set and fixed on the inlet nodes every step;
  - scheme "euler" (default): BDF1 coefficients in the fluid every step, BDF1 mesh velocity, backward
    Euler in the plate, as for the tube. The steady state does not depend on the scheme, and with
    dynamic_tau = 0 not on dt either;
  - forces on the half plate, three ways:
      Fr = -sum of the fluid REACTION over the interface nodes (the load the FSI coupling gives the
           plate); at the nodes of the plate base edge (y = 0) it also contains the share of the floor
           traction next to the edge;
      Fc = Fr minus that floor share, int N_i sigma n dA over the floor faces at the edge nodes, from the
           P1 velocity gradient of each face's tetrahedron and the linear pressure (primary value);
      Ff = int sigma n dA over the interface faces directly, from the same P1 fields (cross-check);
  - u_x, u_y, u_z at A = (0.45, 0.15, -0.15), a plate node;
  - stops when u_x(A) and Fc_x change by less than steady_tol (relative) over steady_window seconds,
    after steady_after seconds, or at end_time;
  - history.csv and run_info.json as for the tube.

Usage: python MainKratos_beam.py [ProjectParameters.json] [--end-time T]
"""
import importlib
import json
import math
import os
import platform
import socket
import sys
import time

import numpy as np
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


def set_bdf1(process_info):
    dt = process_info[KM.DELTA_TIME]
    process_info[KM.BDF_COEFFICIENTS] = KM.Vector([1.0 / dt, -1.0 / dt, 0.0])


def ramp(t, t_ramp):
    return 0.5 * (1.0 - math.cos(math.pi * t / t_ramp)) if t < t_ramp else 1.0


class FaceForces:
    """Traction integrals over triangles of the fluid boundary, from the P1 fields of the tetrahedron
    each triangle belongs to. Triangles are oriented out of the fluid; the force on the plate (or
    floor) is int sigma n_wall dA with n_wall = -n_fluid."""

    def __init__(self, fluid_mp, faces, mu, weight_nodes=None):
        self.mu = mu
        self.nodes = {}
        self.tri = np.array([f[:3] for f in faces], dtype=int)
        tets = []
        for f in faces:
            tets.append([n.Id for n in fluid_mp.GetElement(f[3]).GetNodes()])
        self.tet = np.array(tets, dtype=int)
        ids = sorted(set(self.tri.ravel()) | set(self.tet.ravel()))
        self.ids = np.array(ids)
        self.pos = {n: i for i, n in enumerate(ids)}
        self.node_list = [fluid_mp.GetNode(n) for n in ids]
        self.tri_i = np.vectorize(self.pos.get)(self.tri)
        self.tet_i = np.vectorize(self.pos.get)(self.tet)
        # weights of the three vertices: 1 for all (total force) or only for the listed nodes (share)
        if weight_nodes is None:
            self.mask = np.ones(self.tri.shape)
        else:
            w = set(weight_nodes)
            self.mask = np.array([[1.0 if n in w else 0.0 for n in row] for row in self.tri])

    def force(self):
        X = np.array([[n.X, n.Y, n.Z] for n in self.node_list])
        U = np.array([[n.GetSolutionStepValue(KM.VELOCITY_X), n.GetSolutionStepValue(KM.VELOCITY_Y),
                       n.GetSolutionStepValue(KM.VELOCITY_Z)] for n in self.node_list])
        P = np.array([n.GetSolutionStepValue(KM.PRESSURE) for n in self.node_list])
        xt = X[self.tet_i]                                   # (F, 4, 3)
        J = np.transpose(xt[:, 1:, :] - xt[:, :1, :], (0, 2, 1))  # columns x_k - x_0
        Jinv = np.linalg.inv(J)                              # rows are grad N_1..3
        dN = np.concatenate([-Jinv.sum(axis=1, keepdims=True), Jinv], axis=1)  # (F, 4, 3)
        G = np.einsum("fki,fkj->fij", U[self.tet_i], dN)     # G_ij = du_i/dx_j
        S = self.mu * (G + np.transpose(G, (0, 2, 1)))
        xf = X[self.tri_i]
        area_vec = 0.5 * np.cross(xf[:, 1] - xf[:, 0], xf[:, 2] - xf[:, 0])   # A n_fluid
        pf = P[self.tri_i]
        # int N_i p dA = A (2 p_i + p_j + p_k) / 12, int N_i dA = A / 3; weighted by the mask
        m = self.mask
        p_int = (m * (pf + pf.sum(axis=1, keepdims=True))).sum(axis=1) / 12.0
        w_int = m.sum(axis=1) / 3.0
        visc = np.einsum("fij,fj->fi", S, area_vec)
        # force on the wall = int (p n_f - S n_f) dA
        return (p_int[:, None] * area_vec - w_int[:, None] * visc).sum(axis=0)


def create_stage(cls, model, parameters, case):

    class RecordedBeamAnalysis(cls):

        def Initialize(self):
            super().Initialize()
            solver = self._GetSolver()
            if case["scheme"] == "euler":
                scheme_name = type(solver.fluid_solver._GetScheme()).__name__
                if scheme_name != "BDF2TurbulentScheme":
                    raise RuntimeError(f"expected BDF2TurbulentScheme in the fluid, found {scheme_name}")
                if not hasattr(solver.mesh_solver, "time_int_helper"):
                    raise RuntimeError("mesh solver has no time_int_helper to replace")
                solver.mesh_solver.time_int_helper = KM.TimeDiscretization.BDF1()

            structure = self.model.GetModelPart("Structure")
            a = case["point_a"]
            self.node_a = min(structure.Nodes, key=lambda n: (n.X0 - a[0])**2 + (n.Y0 - a[1])**2 + (n.Z0 - a[2])**2)
            if math.dist((self.node_a.X0, self.node_a.Y0, self.node_a.Z0), a) > 1e-10:
                raise RuntimeError(f"no plate node at A = {a}")

            fluid = self.model.GetModelPart("FluidModelPart")
            self.inlet = self.model.GetModelPart("FluidModelPart.Inlet")
            self.interface = self.model.GetModelPart("FluidModelPart.FluidInterface")
            with open("force_faces.json") as f:
                ff = json.load(f)
            self.face_int = FaceForces(fluid, ff["interface"], case["mu"])
            self.floor_share = FaceForces(fluid, ff["floor_at_edge"], case["mu"], weight_nodes=ff["edge_nodes"])

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
            cols = ["time", "coupling_solves", "converged", "rms_residual_at_acceptance", "wall_step_s", "u_peak_inlet",
                    "ux_A", "uy_A", "uz_A"]
            for tag in ("Fr", "Fc", "Ff"):
                cols += [f"{tag}x", f"{tag}y", f"{tag}z"]
            self.history.write(",".join(cols) + "\n")
            self.track = []
            self.stop_reason = "end_time"
            self.t_wall0 = time.time()
            self.t_step0 = self.t_wall0
            sys.stdout.flush()

        def ApplyBoundaryConditions(self):
            super().ApplyBoundaryConditions()
            self.u_now = case["inlet_peak"] * ramp(self.time, case["ramp_s"])
            ymax, zh = case["inlet"]["y_max"], case["inlet"]["z_half"]
            scale = self.u_now / ((0.5 * ymax) ** 2 * zh ** 2)
            for node in self.inlet.Nodes:
                y, z = node.Y0, node.Z0
                node.SetSolutionStepValue(KM.VELOCITY_X, max(0.0, scale * y * (ymax - y) * (zh * zh - z * z)))
                node.SetSolutionStepValue(KM.VELOCITY_Y, 0.0)
                node.SetSolutionStepValue(KM.VELOCITY_Z, 0.0)
                node.Fix(KM.VELOCITY_X)
                node.Fix(KM.VELOCITY_Y)
                node.Fix(KM.VELOCITY_Z)

        def InitializeSolutionStep(self):
            self.step_rec = {}
            self.t_step0 = time.time()
            self.step_index = getattr(self, "step_index", 0) + 1
            mesh = self._GetSolver().mesh_solver
            first_bdf1 = case["scheme"] == "bdf2" and self.step_index == 1
            if case["scheme"] == "bdf2":
                mesh.time_int_helper = KM.TimeDiscretization.BDF1() if first_bdf1 else KM.TimeDiscretization.BDF2()
            super().InitializeSolutionStep()
            if case["scheme"] == "euler" or first_bdf1:
                set_bdf1(self.model.GetModelPart("FluidModelPart").ProcessInfo)

        def FinalizeSolutionStep(self):
            super().FinalizeSolutionStep()
            r = self.step_rec
            rms = r.get("last_norm", float("nan")) / math.sqrt(self.n_dofs)
            fr = np.zeros(3)
            for node in self.interface.Nodes:
                fr -= np.array([node.GetSolutionStepValue(KM.REACTION_X), node.GetSolutionStepValue(KM.REACTION_Y),
                                node.GetSolutionStepValue(KM.REACTION_Z)])
            fc = fr - self.floor_share.force()
            ff = self.face_int.force()
            u = [self.node_a.GetSolutionStepValue(v) for v in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y, KM.DISPLACEMENT_Z)]
            vals = [f"{self.time:.10g}", str(r.get("solves", 0)), str(int(r.get("converged", False))),
                    f"{rms:.6e}", f"{time.time() - self.t_step0:.4f}", f"{self.u_now:.8g}"]
            vals += [f"{v:.12e}" for v in u]
            vals += [f"{v:.12e}" for v in list(fr) + list(fc) + list(ff)]
            self.history.write(",".join(vals) + "\n")
            self.history.flush()
            self.track.append((self.time, u[0], fc[0]))
            sys.stdout.flush()

        def KeepAdvancingSolutionLoop(self):
            if not super().KeepAdvancingSolutionLoop():
                return False
            t, ux, fx = self.track[-1] if self.track else (0.0, 0.0, 0.0)
            window = case["steady_window"]
            if t >= case["steady_after"] and t >= window:
                old = min(self.track, key=lambda rec: abs(rec[0] - (t - window)))
                du = abs(ux - old[1]) / max(abs(ux), 1e-30)
                df = abs(fx - old[2]) / max(abs(fx), 1e-30)
                if du < case["steady_tol"] and df < case["steady_tol"]:
                    self.stop_reason = f"steady at t = {t:.6g} s (change over {window} s: u_x {du:.2e}, F_x {df:.2e})"
                    return False
            return True

        def Finalize(self):
            super().Finalize()
            self.history.close()
            self.wall_total = time.time() - self.t_wall0

    return RecordedBeamAnalysis(model, parameters)


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
        "stop_reason": getattr(simulation, "stop_reason", "unknown"),
        "wall_total_s": round(time.time() - t0, 3),
        "wall_solution_loop_s": round(getattr(simulation, "wall_total", float("nan")), 3),
    }
    with open("run_info.json", "w") as f:
        json.dump(info, f, indent=2)
    print(json.dumps(info, indent=2))
