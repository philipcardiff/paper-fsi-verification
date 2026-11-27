"""Render the 3-D perspective view of the quarter tube (run with pvbatch).

Usage: pvbatch render.py DATADIR OUT.png [--width W] [--height H] [--scale S]
                         [--pmin P0] [--pmax P1] [--urmax U] [--camera CAM.json]
  DATADIR: the .vtp files written by extract.py (undeformed geometry, mm)
  --scale: displacement exaggeration applied to wall and fluid mesh (default 10)
The computed quarter (x, y >= 0) is drawn mirrored about y = 0 and seen
from +y, so that the flow runs left to right. The fluid symmetry plane y = 0 (and the inlet/outlet faces) are coloured by p
(viridis, P0..P1 Pa); the wall surfaces by the radial displacement u_r
(Oranges, 0..U mm). Also writes OUT.json with image coordinates of a few
points, for the labels drawn by fig.py.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from paraview import servermanager
from paraview.simple import (XMLPolyDataReader, WarpByVector, GetActiveViewOrCreate, Show,
                             GetColorTransferFunction, SaveScreenshot, Render, FeatureEdges,
                             Transform)
from vtkmodules.vtkRenderingCore import vtkCoordinate

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("data")
ap.add_argument("out")
ap.add_argument("--width", type=int, default=3200)
ap.add_argument("--height", type=int, default=1200)
ap.add_argument("--scale", type=float, default=10.0)
ap.add_argument("--pmin", type=float, default=0.0)
ap.add_argument("--pmax", type=float, default=1500.0)
ap.add_argument("--urmax", type=float, default=0.16)
ap.add_argument("--camera", default=str(HERE / "camera.json"))
a = ap.parse_args()
cfg = json.load(open(a.camera))

v = GetActiveViewOrCreate("RenderView")
v.ViewSize = [a.width, a.height]
v.Background = [1, 1, 1]
v.UseColorPaletteForBackground = 0
v.OrientationAxesVisibility = 0
v.UseFXAA = 1
v.KeyLightIntensity = 0.7


def rd(n):
    """Reader, mirrored y -> -y: the computed quarter (x, y >= 0) is shown as
    its mirror image about the symmetry plane y = 0 (the quarter x >= 0, y <= 0,
    same solution by symmetry), so that, seen from +y, the flow runs left to right."""
    t = Transform(Input=XMLPolyDataReader(FileName=f"{a.data}/{n}.vtp"))
    t.Transform.Scale = [1, -1, 1]
    return t


def warped(n):
    w = WarpByVector(Input=rd(n))
    w.Vectors = ["POINTS", "disp"]
    w.ScaleFactor = a.scale
    return w


def lut(name, preset, lo, hi):
    L = GetColorTransferFunction(name)
    L.ApplyPreset(preset, True)
    L.AutomaticRescaleRangeMode = "Never"
    L.RescaleTransferFunction(lo, hi)
    return L


Lp = lut("p", "Viridis", a.pmin, a.pmax)
Lu = lut("ur", "Oranges", 0.0, a.urmax)
Lu.InvertTransferFunction()        # light (u_r = 0) to dark orange (peak), as matplotlib "Oranges"
for n in ("fluid_symy", "fluid_symx", "fluid_inlet", "fluid_outlet"):
    d = Show(warped(n), v)
    d.ColorArrayName = ["POINTS", "p"]
    d.LookupTable = Lp
    d.Ambient, d.Diffuse, d.Specular = 1.0, 0.0, 0.0     # unlit: colours match the colour bar
for n in ("solid", "solid_symy", "solid_symx", "solid_inlet", "solid_outlet"):
    w = warped(n)
    d = Show(w, v)
    d.ColorArrayName = ["POINTS", "ur"]
    d.LookupTable = Lu
    d.Ambient, d.Diffuse, d.Specular = 0.3, 0.75, 0.15
    e = FeatureEdges(Input=w)
    e.BoundaryEdges, e.FeatureEdges, e.NonManifoldEdges, e.ManifoldEdges = 1, 0, 0, 0
    d = Show(e, v)
    d.ColorArrayName = ["POINTS", ""]
    d.AmbientColor = d.DiffuseColor = [0.25, 0.25, 0.25]
    d.LineWidth = 1.5
# undeformed outer wall and inner wall outline on the plane y = 0, for scale
for n in ("solid_symy",):
    e = FeatureEdges(Input=rd(n))
    e.BoundaryEdges, e.FeatureEdges, e.NonManifoldEdges, e.ManifoldEdges = 1, 0, 0, 0
    d = Show(e, v)
    d.ColorArrayName = ["POINTS", ""]
    d.AmbientColor = d.DiffuseColor = [0.0, 0.0, 0.0]
    d.LineWidth = 1.2
    d.Opacity = 0.6

cam = v.GetActiveCamera()
cam.SetFocalPoint(*cfg["focal"])
cam.SetPosition(*cfg["pos"])
cam.SetViewUp(*cfg["up"])
cam.SetViewAngle(cfg.get("angle", 30))
Render(v)
SaveScreenshot(a.out, v, ImageResolution=[a.width, a.height], TransparentBackground=0)

# image coordinates (pixels from the top left) of a few world points, for fig.py
from vtkmodules.util.numpy_support import vtk_to_numpy
s = servermanager.Fetch(XMLPolyDataReader(FileName=f"{a.data}/solid.vtp"))
X = vtk_to_numpy(s.GetPoints().GetData())
ur = vtk_to_numpy(s.GetPointData().GetArray("ur"))
D = vtk_to_numpy(s.GetPointData().GetArray("disp"))
ipk = np.argmax(ur)
world = {
    "urPeak": ((X[ipk] + a.scale*D[ipk])*[1, -1, 1]).tolist(),
    "plane": [2.5, 0.0, 12.0], "inlet": [2.5, -2.5, 0.0], "outlet": [3.0, 0.0, 50.0],
    "wallMid": [6.0*np.cos(0.3), -6.0*np.sin(0.3), 25.0],
    "axis0": [0.0, 0.0, 0.0], "axis1": [0.0, 0.0, 50.0],
}
o = np.array([-3.0, 0.0, -6.0])
for k, e in (("o", [0, 0, 0]), ("ex", [3, 0, 0]), ("ey", [0, -3, 0]), ("ez", [0, 0, 3])):
    world["triad_" + k] = (o + e).tolist()
ren = v.GetRenderer()
c = vtkCoordinate()
c.SetCoordinateSystemToWorld()
px = {}
for k, w in world.items():
    c.SetValue(*w)
    x, y = c.GetComputedDoubleDisplayValue(ren)
    px[k] = [x, a.height - y]
json.dump({"size": [a.width, a.height], "px": px, "world": world, "scale": a.scale,
           "pmin": a.pmin, "pmax": a.pmax, "urmax": a.urmax,
           "urPeak_mm": float(ur[ipk]), "urPeak_z_mm": float(X[ipk, 2])},
          open(str(Path(a.out).with_suffix(".json")), "w"), indent=1)
