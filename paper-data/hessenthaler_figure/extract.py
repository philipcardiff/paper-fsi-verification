"""Extract the small surface data for the Hessenthaler 3-D figure from a
reconstructed solids4foam case (run with ParaView's pvbatch).

Usage: pvbatch extract.py CASE.foam OUTDIR [--time 45]
  CASE.foam: empty stub file in the case directory (reconstructed time dir)
  OUTDIR:    receives wall.vtp, inlets.vtp, outlet.vtp, flap0.vtp,
             flap.vtp (deformed, with |D|) and slice_x0.vtp (with |U|)
All coordinates are written in mm.
"""
import argparse
from paraview.simple import (OpenFOAMReader, ExtractBlock, MergeBlocks, ExtractSurface,
                             Transform, Calculator, WarpByVector, Slice, PassArrays,
                             SaveData, UpdatePipeline)

ap = argparse.ArgumentParser()
ap.add_argument("case")
ap.add_argument("out")
ap.add_argument("--time", type=float, default=45.0)
a = ap.parse_args()
t = a.time


def reader(regions):
    r = OpenFOAMReader(FileName=a.case)
    r.MeshRegions = regions
    r.CellArrays = ["U", "D", "p"]
    r.Createcelltopointfiltereddata = 1
    r.UpdatePipeline(t)
    return r


def mm(src):
    tr = Transform(Input=src)
    tr.Transform.Scale = [1e3, 1e3, 1e3]
    return tr


def save(src, name, arrays=()):
    s = PassArrays(Input=src)
    s.PointDataArrays = list(arrays)
    s.CellDataArrays = []
    s.FieldDataArrays = []
    s.UpdatePipeline(t)
    SaveData(f"{a.out}/{name}.vtp", proxy=s, CompressorType="ZLib", DataMode="Binary")


def surf(regions, name, arrays=()):
    m = MergeBlocks(Input=reader(regions))
    save(mm(ExtractSurface(Input=m)), name, arrays)


surf(["/fluid/patch/wall"], "wall")
surf(["/fluid/patch/upperInlet", "/fluid/patch/lowerInlet"], "inlets")
surf(["/fluid/patch/outlet"], "outlet")

# solid: undeformed and deformed (warped by pointD), with |D|
sol = MergeBlocks(Input=reader(["/solid/internalMesh"]))
ss = ExtractSurface(Input=sol)
save(mm(ss), "flap0")
w = WarpByVector(Input=ss)
w.Vectors = ["POINTS", "pointD"]
c = Calculator(Input=w)
c.ResultArrayName = "magD"
c.Function = "mag(pointD)*1e3"
save(mm(c), "flap", ["magD"])

# fluid mid-plane slice x = 0 (the fluid mesh at time t has moved with the flap)
fl = MergeBlocks(Input=reader(["/fluid/internalMesh"]))
sl = Slice(Input=fl)
sl.SliceType = "Plane"
sl.SliceType.Origin = [0, 0, 0]
sl.SliceType.Normal = [1, 0, 0]
sl.Triangulatetheslice = 0
c2 = Calculator(Input=sl)
c2.ResultArrayName = "magU"
c2.Function = "mag(U)"
save(mm(c2), "slice_x0", ["magU"])
