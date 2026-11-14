"""Generic PETSc SNES driver for a single (possibly mixed) UFL residual.

Jacobian is derived automatically with ufl.derivative (consistent tangent).
Linear solves: direct (MUMPS) by default, configurable through PETSc options.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import ufl
from petsc4py import PETSc
from dolfinx import fem
from dolfinx.fem import petsc as fpetsc


@dataclass
class SNESResult:
    converged: bool
    reason: int
    iterations: int
    fnorm: float


class SNESSolver:
    """Solve F(w; test)=0 for w (fem.Function) with Dirichlet bcs."""

    def __init__(self, F, w: fem.Function, bcs, prefix="fsi_", rtol=1e-9, atol=1e-10,
                 max_it=25, linesearch="bt", options=None, J=None,
                 F_masked=None, masked_dofs=None):
        self.w, self.bcs = w, list(bcs)
        V = w.function_space
        if J is None:
            J = ufl.derivative(F, w, ufl.TrialFunction(V))
        self.F_form = fem.form(F)
        self.J_form = fem.form(J)
        # Optional extra residual part whose ROWS at `masked_dofs` are dropped (used for the
        # harmonic mesh equation: interface mesh-displacement rows are given by the kinematic eq.)
        self.Fm_form = self.Jm_form = None
        self.masked_dofs = None
        if F_masked is not None:
            self.Fm_form = fem.form(F_masked)
            self.Jm_form = fem.form(ufl.derivative(F_masked, w, ufl.TrialFunction(V)))
            md = np.asarray(masked_dofs, dtype=np.int64)
            im, bs = V.dofmap.index_map, V.dofmap.index_map_bs
            nown = im.size_local * bs
            self.masked_dofs = md[md < nown].astype(np.int32)             # owned local (vector entries)
            gl = im.local_to_global((md // bs).astype(np.int32)).astype(np.int64) * bs + md % bs
            self.masked_rows_global = np.unique(gl).astype(np.int32)      # global rows (matrix)
            self.Am = fpetsc.create_matrix(self.Jm_form)
            self.bm = fpetsc.create_vector(V) if not _has_create_vector_form() else \
                fpetsc.create_vector(self.Fm_form)
        self.A = fpetsc.create_matrix(self.J_form)
        self.b = fpetsc.create_vector(self.F_form) if _has_create_vector_form() else \
            fpetsc.create_vector(V)
        snes = PETSc.SNES().create(w.function_space.mesh.comm)
        snes.setOptionsPrefix(prefix)
        opts = PETSc.Options()
        defaults = {"snes_type": "newtonls", "snes_linesearch_type": linesearch,
                    "snes_rtol": rtol, "snes_atol": atol, "snes_max_it": max_it,
                    "ksp_type": "preonly", "pc_type": "lu",
                    "pc_factor_mat_solver_type": "mumps"}
        defaults.update(options or {})
        for k, v in defaults.items():
            opts[f"{prefix}{k}"] = v
        snes.setFunction(self._residual, self.b)
        snes.setJacobian(self._jacobian, self.A)
        snes.setFromOptions()
        self.snes = snes

    # --- callbacks -------------------------------------------------------
    def _residual(self, snes, x, b):
        x.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        x.copy(self.w.x.petsc_vec)
        self.w.x.scatter_forward()
        with b.localForm() as bl:
            bl.set(0.0)
        fpetsc.assemble_vector(b, self.F_form)
        if self.Fm_form is not None:
            bm = self.bm
            with bm.localForm() as bl:
                bl.set(0.0)
            fpetsc.assemble_vector(bm, self.Fm_form)
            bm.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
            bm.array[self.masked_dofs] = 0.0
            b.axpy(1.0, bm)
        fpetsc.apply_lifting(b, [self.J_form], bcs=[self.bcs], x0=[x], alpha=-1.0)
        b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        fpetsc.set_bc(b, self.bcs, x, -1.0)

    def _jacobian(self, snes, x, A, P):
        x.ghostUpdate(addv=PETSc.InsertMode.INSERT, mode=PETSc.ScatterMode.FORWARD)
        x.copy(self.w.x.petsc_vec)
        self.w.x.scatter_forward()
        A.zeroEntries()
        fpetsc.assemble_matrix(A, self.J_form, bcs=self.bcs)
        A.assemble()
        if self.Jm_form is not None:
            Am = self.Am
            Am.zeroEntries()
            fpetsc.assemble_matrix(Am, self.Jm_form, bcs=self.bcs, diag=0.0)
            Am.assemble()
            Am.zeroRows(self.masked_rows_global, diag=0.0)
            A.axpy(1.0, Am, structure=PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN)

    def solve(self) -> SNESResult:
        x = self.w.x.petsc_vec.copy()
        # make initial guess satisfy bcs
        fpetsc.set_bc(x, self.bcs)
        self.snes.solve(None, x)
        x.copy(self.w.x.petsc_vec)
        self.w.x.scatter_forward()
        r = self.snes.getConvergedReason()
        out = SNESResult(r > 0, r, self.snes.getIterationNumber(),
                         self.snes.getFunctionNorm())
        x.destroy()
        return out


def _has_create_vector_form():
    import inspect
    try:
        sig = inspect.signature(fpetsc.create_vector)
        p = list(sig.parameters)
        return "L" in p or "form" in p
    except Exception:
        return False
