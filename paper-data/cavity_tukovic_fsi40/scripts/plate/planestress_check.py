"""Plate-only check: StVK and linear with the plane-stress lambda of the
foam-extend fsi constitutiveModel (lambda* = nu E/((1+nu)(1-nu)))."""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "cavity_plate_material_law"))
import plate_laws as pl
q = 1.315
for tag, lam in (("planeStrain", pl.LAM), ("planeStress", pl.E*pl.NU/((1+pl.NU)*(1-pl.NU)))):
    pl.LAM = lam
    for law in ("stvk", "linear"):
        r = pl.run(20, 2, law, q, nsteps=20 if law != "linear" else 1)
        print(tag, law, lam, r, flush=True)
