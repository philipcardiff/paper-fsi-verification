"""Extract the small surface data for the 3-D tube perspective figure from a
reconstructed solids4foam 3dTube case (run with ParaView's pvbatch).

Usage: pvbatch extract.py CASE.foam INIT.foam OUTDIR [--time 0.005] [--rho 1000]
  CASE.foam: empty stub file in the case directory, which holds the
             reconstructed time directory (with the moved fluid mesh)
  INIT.foam: stub file in a directory holding only (links to) the case's
             constant/ and system/, so that the reader returns the initial mesh
             (the OpenFOAM reader does not list a 0 directory without fields
             for every region)
  OUTDIR:    receives, all on the UNDEFORMED geometry (coordinates in mm):
     fluid_symy.vtp, fluid_symx.vtp, fluid_inlet.vtp, fluid_outlet.vtp
         with p (Pa, = rho * kinematic p) and disp (mesh motion, mm)
     solid.vtp  (outer surface of the wall) with disp (pointD, mm) and
         ur (radial displacement, mm)
  so that render.py can warp both by the same exaggeration factor.
"""
import argparse

import numpy as np
from paraview import servermanager
from paraview.simple import OpenFOAMReader, MergeBlocks, ExtractSurface, GetTimeKeeper
from vtkmodules.util.numpy_support import numpy_to_vtk, vtk_to_numpy
from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter

ap = argparse.ArgumentParser()
ap.add_argument("case")
ap.add_argument("init")
ap.add_argument("out")
ap.add_argument("--time", type=float, default=0.005)
ap.add_argument("--rho", type=float, default=1000.0)
a = ap.parse_args()


def fetch(regions, t, arrays=("p", "D"), case=None):
    r = OpenFOAMReader(FileName=case or a.case)
    r.MeshRegions = regions
    r.CellArrays = list(arrays)
    r.SkipZeroTime = 0
    r.Createcelltopointfiltereddata = 1
    s = ExtractSurface(Input=MergeBlocks(Input=r))
    GetTimeKeeper().Time = t          # Fetch() updates at the time-keeper time
    s.UpdatePipeline(t)
    return servermanager.Fetch(s)


def write(pd, name, arrays):
    from vtkmodules.vtkCommonDataModel import vtkPolyData
    out = vtkPolyData()
    out.DeepCopy(pd)
    pts = vtk_to_numpy(out.GetPoints().GetData()).copy()
    out.GetPoints().SetData(numpy_to_vtk(pts*1e3, deep=True))
    for k in [out.GetPointData().GetArrayName(i) for i in range(out.GetPointData().GetNumberOfArrays())]:
        out.GetPointData().RemoveArray(k)
    for k in [out.GetCellData().GetArrayName(i) for i in range(out.GetCellData().GetNumberOfArrays())]:
        out.GetCellData().RemoveArray(k)
    for k, v in arrays.items():
        arr = numpy_to_vtk(np.ascontiguousarray(v, dtype=np.float32), deep=True)
        arr.SetName(k)
        out.GetPointData().AddArray(arr)
    w = vtkXMLPolyDataWriter()
    w.SetFileName(f"{a.out}/{name}.vtp")
    w.SetInputData(out)
    w.SetDataModeToBinary()
    w.SetCompressorTypeToZLib()
    w.Write()


for patch, name in (("symmetry-y", "fluid_symy"), ("symmetry-x", "fluid_symx"),
                    ("inlet", "fluid_inlet"), ("outlet", "fluid_outlet")):
    reg = [f"/fluid/patch/{patch}"]
    p0, p1 = fetch(reg, 0.0, (), a.init), fetch(reg, a.time, ("p",))
    X0 = vtk_to_numpy(p0.GetPoints().GetData())
    X1 = vtk_to_numpy(p1.GetPoints().GetData())
    assert X0.shape == X1.shape
    p = vtk_to_numpy(p1.GetPointData().GetArray("p"))*a.rho
    print(name, "max |mesh motion| (mm):", abs(X1 - X0).max()*1e3)
    write(p0, name, {"p": p, "disp": (X1 - X0)*1e3})

for patch, name in (("outer-wall", "solid"), ("symmetry-y", "solid_symy"),
                    ("symmetry-x", "solid_symx"), ("inlet", "solid_inlet"),
                    ("outlet", "solid_outlet")):
    s = fetch([f"/solid/patch/{patch}"], a.time, ("D",))
    X = vtk_to_numpy(s.GetPoints().GetData())
    Dp = vtk_to_numpy(s.GetPointData().GetArray("D"))
    r = np.hypot(X[:, 0], X[:, 1])
    ur = (X[:, 0]*Dp[:, 0] + X[:, 1]*Dp[:, 1])/np.maximum(r, 1e-12)
    write(s, name, {"disp": Dp*1e3, "ur": ur*1e3})
print("t =", a.time)
