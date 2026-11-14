"""GATE B1 material-law verification (serial). Writes JSON + prints PASS/FAIL against the preregistered thresholds."""
import sys, json
import numpy as np
sys.path.insert(0, __import__("os").environ.get("HESS_ROOT", "."))
import ufl, basix.ufl
from dolfinx import mesh, fem
from mpi4py import MPI
from src.hessenthaler import nh3d as nh

rng = np.random.default_rng(1234)
def rand_F(scale):
    while True:
        F = np.eye(3) + scale * rng.standard_normal((3, 3))
        if np.linalg.det(F) > 0.3: return F
def rand_Q():
    q, r = np.linalg.qr(rng.standard_normal((3, 3))); q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0: q[:, 0] *= -1
    return q
res = {}
# UFL P_ufl vs P_analytic vs numpy, evaluated at quadrature of a single tet with constant F
msh = mesh.create_unit_cube(MPI.COMM_SELF, 1, 1, 1, mesh.CellType.tetrahedron)
V = fem.functionspace(msh, ("DG", 0, (3, 3)))
errs, errs_np = [], []
Fs = [rand_F(s) for s in (0.05, 0.2, 0.4, 0.6) for _ in range(5)]
for F in Fs:
    Fc = fem.Constant(msh, F)
    out = []
    for expr in (nh.P_ufl(Fc), nh.P_analytic(Fc)):
        f = fem.Function(V); f.interpolate(fem.Expression(expr, V.element.interpolation_points)); out.append(f.x.array[:9].reshape(3, 3).copy())
    ref = nh.np_P(F)
    errs.append(np.abs(out[0] - out[1]).max() / np.abs(ref).max())
    errs_np.append(max(np.abs(out[1] - ref).max(), np.abs(out[0] - ref).max()) / np.abs(ref).max())
res["Pufl_vs_Panalytic_rel"] = float(max(errs)); res["P_vs_numpy_rel"] = float(max(errs_np))
# P = dW/dF by central FD of the numpy energy
fd = []
for F in Fs:
    h = 1e-6; G = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            E = np.zeros((3, 3)); E[i, j] = h; G[i, j] = (nh.np_W(F + E) - nh.np_W(F - E)) / (2 * h)
    fd.append(np.abs(G - nh.np_P(F)).max() / np.abs(nh.np_P(F)).max())
res["P_vs_dW_fd_rel"] = float(max(fd))
# consistent tangent: UFL derivative of P_analytic vs central FD of P (symmetry of the tangent A_ijkl = A_klij)
tt, sym = [], []
for F in Fs:
    Fv = ufl.variable(ufl.as_matrix(F.tolist())) if False else None
    A_fd = nh.np_tangent_fd(F, h=1e-6)
    # analytic tangent by differentiating the ufl P wrt F with a constant-matrix variable, via directional derivative
    Fc = fem.Constant(msh, F)
    cols = []
    for k in range(9):
        E = np.zeros(9); E[k] = 1.0; Ec = fem.Constant(msh, E.reshape(3, 3))
        # dP[E] = d/de P(F + e E) at e = 0 using a UFL variable
        e = ufl.variable(fem.Constant(msh, 0.0))
        Pe = nh.P_analytic(Fc + e * Ec); dP = ufl.diff(Pe, e)
        f = fem.Function(V); f.interpolate(fem.Expression(dP, V.element.interpolation_points)); cols.append(f.x.array[:9].copy())
    A = np.array(cols).T
    tt.append(np.abs(A - A_fd).max() / np.abs(A).max()); sym.append(np.abs(A - A.T).max() / np.abs(A).max())
res["tangent_vs_FD_rel"] = float(max(tt)); res["tangent_major_symmetry_rel"] = float(max(sym))
# frame indifference and rotation
fi, rot, si = [], [], []
for F in Fs:
    Q = rand_Q()
    fi.append(abs(nh.np_W(Q @ F) - nh.np_W(F)) / abs(nh.np_W(F)))   # W(QF)=W(F)
    s1, s2 = nh.np_cauchy(Q @ F), Q @ nh.np_cauchy(F) @ Q.T
    si.append(np.abs(s1 - s2).max() / np.abs(s2).max())
    rot.append(np.abs(nh.np_cauchy(Q)).max() / nh.MU)
    # sigma = P F^T / J consistency
    J = np.linalg.det(F); si.append(np.abs(nh.np_P(F) @ F.T / J - nh.np_cauchy(F)).max() / np.abs(nh.np_cauchy(F)).max())
res["W_frame_indiff_rel"] = float(max(fi)); res["sigma_covariance_and_PFt_rel"] = float(max(si)); res["rotation_stress_over_mu"] = float(max(rot))
# stress-free reference, small-strain limit: P ~ 2 mu eps + lambda tr(eps) I with lambda = K - 2mu/3
eps = 1e-6; Fs_ = np.eye(3) + eps * np.array([[0.3, 0.1, 0.0], [0.2, -0.2, 0.05], [0.0, 0.4, 0.5]])
e_ = 0.5 * (Fs_ + Fs_.T) - np.eye(3); lam = nh.K_BULK - 2 * nh.MU / 3
res["small_strain_limit_rel"] = float(np.abs(nh.np_P(Fs_) - (2 * nh.MU * e_ + lam * np.trace(e_) * np.eye(3))).max() / (2 * nh.MU * eps))
res["P_at_I"] = float(np.abs(nh.np_P(np.eye(3))).max())
# volumetric: pure dilation F = s I  -> compare with analytic p
s = 1.1; res["bulk_pressure_check"] = float(abs(np.trace(nh.np_cauchy(s * np.eye(3))) / 3 - 0.5 * nh.K_BULK * (s ** 6 - 1) / s ** 3))
res["K_bulk"] = nh.K_BULK
thr = {"Pufl_vs_Panalytic_rel": 1e-10, "P_vs_numpy_rel": 1e-10, "tangent_vs_FD_rel": 1e-6, "W_frame_indiff_rel": 1e-12,
       "sigma_covariance_and_PFt_rel": 1e-10, "rotation_stress_over_mu": 1e-9}
ok = True
for k, t in thr.items():
    st = res[k] <= t; ok &= st; print(f"{k:34s} {res[k]:.3e} <= {t:.0e}  {'PASS' if st else 'FAIL'}")
for k in ("P_vs_dW_fd_rel", "tangent_major_symmetry_rel", "small_strain_limit_rel", "P_at_I", "bulk_pressure_check", "K_bulk"):
    print(f"{k:34s} {res[k]:.3e} (reported)")
res["GATE_B1"] = "PASS" if ok else "FAIL"; print("GATE B1:", res["GATE_B1"])
json.dump(res, open(sys.argv[1] if len(sys.argv) > 1 else "test_material.json", "w"), indent=1)
