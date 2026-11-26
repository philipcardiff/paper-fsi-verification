"""Small readers for OpenFOAM/solids4foam cases, shared by the panel scripts.

They use the VTK Python bindings (vtkOpenFOAMReader), so no ParaView is
needed. A region ("fluid" or "solid") of a reconstructed multi-region case is
cut by a plane, or one of its patches is extracted, and returned as a
matplotlib triangulation with point fields (cell fields are interpolated to
the points by the reader)."""
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from matplotlib.tri import Triangulation


def _reader(case, time):
    r = vtk.vtkOpenFOAMReader()
    r.SetFileName(f"{case}/system/controlDict")
    r.CreateCellToPointOn()
    r.ReadZonesOff()
    r.UpdateInformation()
    times = vtk_to_numpy(r.GetTimeValues())
    if time is None:
        time = times.max()
    k = int(np.argmin(np.abs(times - time)))
    if abs(times[k] - time) > 1e-6 * max(1.0, abs(time)):
        raise RuntimeError(f"time {time} not found in {case}: {times}")
    r.UpdateTimeStep(times[k])
    return r, times[k]


def _select(r, wanted):
    r.DisableAllPatchArrays()
    for i in range(r.GetNumberOfPatchArrays()):
        name = r.GetPatchArrayName(i)
        if name == wanted:
            r.SetPatchArrayStatus(name, 1)
            return
    names = [r.GetPatchArrayName(i) for i in range(r.GetNumberOfPatchArrays())]
    raise RuntimeError(f"{wanted} not in {names}")


def _block(r):
    r.Update()
    it = r.GetOutput().NewIterator()
    it.InitTraversal()
    while not it.IsDoneWithTraversal():
        b = it.GetCurrentDataObject()
        if b is not None and b.GetNumberOfCells() > 0:
            return b
        it.GoToNextItem()
    raise RuntimeError("empty selection")


def _to_tri(poly, axes, fields):
    tf = vtk.vtkTriangleFilter()
    tf.SetInputData(poly)
    tf.Update()
    out = tf.GetOutput()
    pts = vtk_to_numpy(out.GetPoints().GetData()).astype(float)
    conn = vtk_to_numpy(out.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    data = {}
    pd = out.GetPointData()
    for f in fields:
        a = pd.GetArray(f)
        if a is None:
            avail = [pd.GetArrayName(i) for i in range(pd.GetNumberOfArrays())]
            raise RuntimeError(f"field {f} missing; have {avail}")
        data[f] = vtk_to_numpy(a).astype(float)
    return pts, Triangulation(pts[:, axes[0]], pts[:, axes[1]], conn), data


def plane_cut(case, region, time=None, fields=(), origin=None,
              normal=(0, 0, 1), axes=(0, 1)):
    """Cut the internal mesh of `region` with a plane (default: the mid-plane
    z = const of a 2-D case). Returns (points (N,3), triangulation in the
    `axes` coordinates, {field: point values}, time)."""
    r, t = _reader(case, time)
    _select(r, f"/{region}/internalMesh")
    grid = _block(r)
    if origin is None:
        b = grid.GetBounds()
        origin = (0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3]), 0.5 * (b[4] + b[5]))
    pl = vtk.vtkPlane()
    pl.SetOrigin(*origin)
    pl.SetNormal(*normal)
    cut = vtk.vtkCutter()
    cut.SetCutFunction(pl)
    cut.SetInputData(grid)
    cut.GenerateTrianglesOff()
    cut.Update()
    pts, tri, data = _to_tri(cut.GetOutput(), axes, fields)
    return pts, tri, data, t


def patch(case, region, name, time=None, fields=(), axes=(0, 2)):
    """Faces of patch `name` of `region` (e.g. a symmetry plane of a 3-D
    case), with point fields."""
    r, t = _reader(case, time)
    _select(r, f"/{region}/patch/{name}")
    blk = _block(r)
    geo = vtk.vtkGeometryFilter()
    geo.SetInputData(blk)
    geo.Update()
    pts, tri, data = _to_tri(geo.GetOutput(), axes, fields)
    return pts, tri, data, t


def boundary_loops(tri):
    """Closed boundary loops of a triangulation (lists of point indices),
    longest first; used to outline the domains."""
    edges = {}
    for t in tri.triangles:
        for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            key = (min(a, b), max(a, b))
            edges[key] = edges.get(key, 0) + 1
    nbr = {}
    for (a, b), n in edges.items():
        if n == 1:
            nbr.setdefault(a, []).append(b)
            nbr.setdefault(b, []).append(a)
    loops, seen = [], set()
    for start in nbr:
        if start in seen:
            continue
        loop, prev, cur = [start], None, start
        seen.add(start)
        while True:
            nxt = [n for n in nbr[cur] if n != prev]
            if not nxt:
                break
            prev, cur = cur, nxt[0]
            if cur == start or cur in seen:
                break
            loop.append(cur)
            seen.add(cur)
        loops.append(loop + [start])
    loops.sort(key=len, reverse=True)
    return loops


def grid_interp(tri, values, x, y):
    """Linear interpolation of point values on a regular grid (masked
    outside the triangulation), for streamplot."""
    from matplotlib.tri import LinearTriInterpolator
    X, Y = np.meshgrid(x, y)
    return LinearTriInterpolator(tri, values)(X, Y)
