"""Render the 3-D perspective view of the Hessenthaler flap (run with pvbatch).

Usage: pvbatch render.py DATADIR OUT.png [--width W] [--height H] [--umax U]
                         [--zmax Z] [--near-opacity A] [--camera CAM.json]
  DATADIR: the .vtp files written by extract.py (coordinates in mm)
  --zmax:  cut the flow domain at z = Z mm (default: no cut)
  CAM.json: {"focal": [...], "pos": [...], "up": [...], "angle": deg,
             "parallel_scale": s (optional)}; default: camera.json next to
             this script
The half of the flow domain away from the camera is drawn as a light shaded
surface, the near half as a faint translucent one; the mid-plane x = 0 is
coloured by |U| (viridis, 0..U m/s); the deformed flap is opaque orange, the
undeformed flap a grey ghost. No text: labels and colour bar come from fig.py.
"""
import argparse
import json
from pathlib import Path

from paraview.simple import (XMLPolyDataReader, Clip, Slice, FeatureEdges, GetActiveViewOrCreate,
                             Show, GetColorTransferFunction, SaveScreenshot, Render)

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("data")
ap.add_argument("out")
ap.add_argument("--width", type=int, default=3200)
ap.add_argument("--height", type=int, default=1300)
ap.add_argument("--umax", type=float, default=0.65)
ap.add_argument("--zmax", type=float, default=None)
ap.add_argument("--near-opacity", type=float, default=0.2)
ap.add_argument("--ring", type=float, default=20.0, help="ring spacing (mm)")
ap.add_argument("--camera", default=str(HERE / "camera.json"))
a = ap.parse_args()
D = a.data
cam_cfg = json.load(open(a.camera))
sgn = 1.0 if cam_cfg["pos"][0] > 0 else -1.0     # camera side of the mid-plane

v = GetActiveViewOrCreate("RenderView")
v.ViewSize = [a.width, a.height]
v.Background = [1, 1, 1]
v.UseColorPaletteForBackground = 0
v.OrientationAxesVisibility = 0
v.UseFXAA = 1
v.KeyLightIntensity = 0.7


def rd(n):
    return XMLPolyDataReader(FileName=f"{D}/{n}.vtp")


def plain(d, rgb, **kw):
    d.ColorArrayName = ["POINTS", ""]
    d.AmbientColor = d.DiffuseColor = rgb
    for k, val in kw.items():
        setattr(d, k, val)
    return d


def zclip(src):
    if a.zmax is None:
        return src
    c = Clip(Input=src)
    c.ClipType = "Plane"
    c.ClipType.Origin = [0, 0, a.zmax]
    c.ClipType.Normal = [0, 0, 1]
    c.Invert = 1
    return c


def half(src, near):
    """Half of src on the camera side (near) or the far side of x = 0."""
    c = Clip(Input=src)
    c.ClipType = "Plane"
    c.ClipType.Origin = [0, 0, 0]
    c.ClipType.Normal = [sgn, 0, 0]
    c.Invert = 0 if near else 1
    return c


def edges(src, angle=None):
    e = FeatureEdges(Input=src)
    e.BoundaryEdges = 1
    e.FeatureEdges = 1 if angle else 0
    if angle:
        e.FeatureAngle = angle
    e.NonManifoldEdges = 0
    e.ManifoldEdges = 0
    return e


wall = zclip(rd("wall"))
plain(Show(half(wall, False), v), [0.88, 0.88, 0.88], Specular=0.1, Ambient=0.25, Diffuse=0.8)
plain(Show(half(wall, True), v), [0.70, 0.72, 0.78], Opacity=a.near_opacity, Specular=0.2)
plain(Show(edges(wall), v), [0.3, 0.3, 0.3], LineWidth=2.0)       # cut rims
# circumferential rings every a.ring mm convey the pipe shape
rings = Slice(Input=wall)
rings.SliceType = "Plane"
rings.SliceType.Origin = [0, 0, 0]
rings.SliceType.Normal = [0, 0, 1]
rings.SliceOffsetValues = [z for z in range(-20, 400, int(a.ring))
                           if a.zmax is None or z < a.zmax]
plain(Show(rings, v), [0.45, 0.47, 0.52], LineWidth=1.2, Opacity=0.8)

INLET, OUTLET = [0.0, 0.447, 0.698], [0.835, 0.369, 0.0]          # Okabe-Ito
plain(Show(rd("inlets"), v), INLET, Ambient=0.3)
plain(Show(edges(rd("inlets")), v), INLET, LineWidth=3.0)
if a.zmax is None:
    plain(Show(rd("outlet"), v), OUTLET, Ambient=0.3)
    plain(Show(edges(rd("outlet")), v), OUTLET, LineWidth=3.0)

d = Show(zclip(rd("slice_x0")), v)
d.ColorArrayName = ["POINTS", "magU"]
lut = GetColorTransferFunction("magU")
lut.ApplyPreset("Viridis", True)
lut.AutomaticRescaleRangeMode = "Never"
lut.RescaleTransferFunction(0.0, a.umax)
d.LookupTable = lut
d.Ambient, d.Diffuse, d.Specular = 1.0, 0.0, 0.0     # unlit: colours match the colour bar
d.Translation = [-sgn*0.05, 0, 0]        # a hair behind the flap mid-plane

plain(Show(rd("flap0"), v), [0.3, 0.3, 0.3], Opacity=0.55)
plain(Show(rd("flap"), v), [0.878, 0.541, 0.173], Specular=0.3, Ambient=0.25)
plain(Show(edges(rd("flap"), 40), v), [0.2, 0.12, 0.05], LineWidth=1.5)

cam = v.GetActiveCamera()
cam.SetFocalPoint(*cam_cfg["focal"])
cam.SetPosition(*cam_cfg["pos"])
cam.SetViewUp(*cam_cfg["up"])
cam.SetViewAngle(cam_cfg.get("angle", 30))
Render(v)
SaveScreenshot(a.out, v, ImageResolution=[a.width, a.height], TransparentBackground=0)

# image coordinates (pixels from the top left) of a few world points, for the
# labels drawn by fig.py; written next to the PNG as OUT.json
import numpy as np
from paraview import servermanager
from vtkmodules.vtkRenderingCore import vtkCoordinate
from vtkmodules.util.numpy_support import vtk_to_numpy
fl = servermanager.Fetch(rd("flap"))
P = vtk_to_numpy(fl.GetPoints().GetData())
tip = P[P[:, 2] > P[:, 2].max() - 0.5].mean(axis=0)
zc = a.zmax if a.zmax is not None else 339.5
world = {
    "upperInlet": [0, 27.15, -29.5], "lowerInlet": [0, -27.15, -29.5],
    "flapTip": tip.tolist(), "flapMid": P[np.argmin(abs(P[:, 2] - 40))].tolist(),
    "flap0Tip": [sgn*5.5, -1.0, 62.0], "flap0Mid": [sgn*5.5, -1.0, 45.0],
    "cutTop": [0, 38.1, zc], "cutNear": [sgn*38.1, 0, zc], "cutBottom": [0, -38.1, zc],
    "wallNear": [sgn*38.1, 0, 120.0], "sliceMid": [0, 0, 130.0],
}
o = np.array([0.0, -60.0, 200.0])
for k, e in (("o", [0, 0, 0]), ("ex", [15, 0, 0]), ("ey", [0, 15, 0]), ("ez", [0, 0, 15])):
    world["triad_" + k] = (o + e).tolist()
ren = v.GetRenderer()
c = vtkCoordinate()
c.SetCoordinateSystemToWorld()
px = {}
for k, w in world.items():
    c.SetValue(*w)
    x, y = c.GetComputedDoubleDisplayValue(ren)
    px[k] = [x, a.height - y]
json.dump({"size": [a.width, a.height], "px": px, "world": world, "zcut": zc,
           "umax": a.umax, "tip": tip.tolist()},
          open(str(Path(a.out).with_suffix(".json")), "w"), indent=1)
