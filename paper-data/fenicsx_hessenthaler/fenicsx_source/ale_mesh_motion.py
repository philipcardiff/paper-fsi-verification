"""Standalone mesh-motion utilities: harmonic extension with stiffness weighting,
Jacobian positivity check."""
import numpy as np
import ufl
from dolfinx import fem
from mpi4py import MPI
from .forms import kin, harmonic_residual


def stiffness_weight(msh, mode="invvol"):
    """DG0 stiffness weight. 'const' -> 1 ; 'invvol' -> (mean cell volume)/(cell volume)
    (small cells stiffer, protects boundary-layer cells)."""
    Q0 = fem.functionspace(msh, ("DG", 0))
    k = fem.Function(Q0)
    if mode == "const":
        k.x.array[:] = 1.0
        return k
    vol = fem.Function(Q0)
    v = fem.form(ufl.TestFunction(Q0) * ufl.dx)
    vv = fem.assemble_vector(v)
    vol.x.array[:] = vv.array
    n = msh.topology.index_map(msh.topology.dim).size_local
    mean = msh.comm.allreduce(vol.x.array[:n].mean(), op=MPI.SUM) / msh.comm.size
    k.x.array[:] = mean / vol.x.array
    return k


def harmonic_extension(V, bcs, k=None, dx=None):
    """Solve  int k grad d : grad psi = 0  with Dirichlet bcs on d (returns Function)."""
    from dolfinx.fem.petsc import LinearProblem
    import inspect
    d = ufl.TrialFunction(V)
    psi = ufl.TestFunction(V)
    dxx = dx if dx is not None else ufl.dx
    kk = k if k is not None else 1.0
    a = kk * ufl.inner(ufl.grad(d), ufl.grad(psi)) * dxx
    L = ufl.inner(fem.Constant(V.mesh, (0.0,) * V.mesh.geometry.dim), psi) * dxx
    out = fem.Function(V)
    kw = dict(bcs=bcs, u=out, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    if "petsc_options_prefix" in inspect.signature(LinearProblem.__init__).parameters:
        kw["petsc_options_prefix"] = "ale_he_"
    LinearProblem(a, L, **kw).solve()
    return out


def jacobian_stats(d, dx=ufl.dx):
    """min/max of J=det(I+grad d) sampled at cell DG0 (midpoint) and P1-interpolated nodes."""
    msh = d.function_space.mesh
    _, J = kin(d)
    Q = fem.functionspace(msh, ("DG", 1))
    Jf = fem.Function(Q)
    Jf.interpolate(fem.Expression(J, Q.element.interpolation_points))
    lo = msh.comm.allreduce(Jf.x.array.min(), op=MPI.MIN)
    hi = msh.comm.allreduce(Jf.x.array.max(), op=MPI.MAX)
    return lo, hi
