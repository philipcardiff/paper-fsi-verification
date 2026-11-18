"""Feature-preserving decimation of the published fluid-domain surface (wall patch only; inlets, outlet, interface kept exact).
The original STL is a high-resolution anisotropic triangulation (123156 tris, 14272 edges < 10 um, areas 2e-9..357 mm^2) on which
gmsh reparametrisation fails. We decimate the wall patch (boundary vertices fixed: inlet/outlet circles and flap hole untouched) and report
the one-sided Hausdorff deviation of the decimated wall from the original (geometry error of the mesh generator input)."""
import sys, os, gzip, re, json
import numpy as np, vtk
from vtk.util import numpy_support as ns
STL = __import__("os").environ.get("HESS_ROOT", ".") + "/docs/s4f_ref_sources/tut/constant/triSurface/fluidDomain.stl.gz"
def read():
    P = {}
    for s in re.split(r"\nendsolid[^\n]*\n?", gzip.open(STL, "rt").read()):
        m = re.match(r"\s*solid (\S+)", s)
        if m: P[m.group(1)] = np.array(re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", s), float).reshape(-1, 3, 3)
    return P
def to_poly(T):
    pts = T.reshape(-1, 3); pd = vtk.vtkPolyData(); p = vtk.vtkPoints(); p.SetData(ns.numpy_to_vtk(pts, deep=True)); pd.SetPoints(p)
    cells = vtk.vtkCellArray()
    for i in range(len(T)):
        cells.InsertNextCell(3); [cells.InsertCellPoint(3 * i + k) for k in range(3)]
    pd.SetPolys(cells); return pd
def main(reduction, out):
    W = read()["wall"]
    pd = to_poly(W)
    cl = vtk.vtkCleanPolyData(); cl.SetInputData(pd); cl.SetTolerance(0.0); cl.ConvertPolysToLinesOff(); cl.Update()   # merge coincident points only
    c0 = cl.GetOutput()
    tri = vtk.vtkTriangleFilter(); tri.SetInputData(c0); tri.PassLinesOff(); tri.Update()
    # remove degenerate by cleaning again
    dec = vtk.vtkQuadricDecimation(); dec.SetInputData(tri.GetOutput()); dec.SetTargetReduction(reduction); dec.SetVolumePreservation(1)
    dec.Update(); d = dec.GetOutput()
    print("original points", c0.GetNumberOfPoints(), "tris", c0.GetNumberOfCells(), "-> decimated points", d.GetNumberOfPoints(), "tris", d.GetNumberOfCells())
    # Hausdorff-type deviation: decimated vertices to original surface
    loc = vtk.vtkCellLocator(); loc.SetDataSet(c0); loc.BuildLocator()
    pts = ns.vtk_to_numpy(d.GetPoints().GetData()); dev = []
    cp = [0.0] * 3; cid = vtk.mutable(0); sid = vtk.mutable(0); d2 = vtk.mutable(0.0)
    for p in pts:
        loc.FindClosestPoint(p, cp, cid, sid, d2); dev.append(np.sqrt(d2.get()))
    dev = np.array(dev)
    # two-sided deviation: original vertices -> decimated surface; decimated centroids -> original surface
    loc2 = vtk.vtkCellLocator(); loc2.SetDataSet(d); loc2.BuildLocator()
    ov = ns.vtk_to_numpy(c0.GetPoints().GetData()); d_o2d = []
    for p in ov[::7]:
        loc2.FindClosestPoint(p, cp, cid, sid, d2); d_o2d.append(np.sqrt(d2.get()))
    d_o2d = np.array(d_o2d)
    tris = ns.vtk_to_numpy(d.GetPolys().GetData()).reshape(-1, 4)[:, 1:]
    T = pts[tris]
    a = 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1)
    cen = T.mean(1); dc = []
    for p in cen:
        loc.FindClosestPoint(p, cp, cid, sid, d2); dc.append(np.sqrt(d2.get()))
    dc = np.array(dc)
    rep = dict(orig_to_dec_max_m=float(d_o2d.max()), orig_to_dec_p99_m=float(np.percentile(d_o2d, 99)), dec_centroid_to_orig_max_m=float(dc.max()), dec_centroid_to_orig_p99_m=float(np.percentile(dc, 99)),
               reduction=reduction, ntri=int(len(T)), 
               area_decimated_m2=float(a.sum()), area_original_m2=float(0.5 * np.linalg.norm(np.cross(W[:, 1] - W[:, 0], W[:, 2] - W[:, 0]), axis=1).sum()),
               min_tri_area_m2=float(a.min()))
    print(rep); np.save(out + ".npy", T); json.dump(rep, open(out + ".json", "w"), indent=1)
if __name__ == "__main__":
    main(float(sys.argv[1]), sys.argv[2])


def assemble(npy, out_stl):
    """Wall (decimated) + snapped boundary loops + planar discs for inlets/outlet + exact flap interface + clamp; returns report."""
    import collections
    from scipy.spatial import Delaunay
    T = np.load(npy)
    key = lambda v: tuple(np.round(v * 1e9).astype(np.int64))
    # snap planar boundaries
    for z0, lo, hi in ((-0.0295, -0.0296, -0.0294), (0.3395, 0.3394, 0.3396)):
        m = (T[:, :, 2] > lo) & (T[:, :, 2] < hi); T[m, 2] = z0
    # boundary loops
    E = collections.defaultdict(list)
    for i, t in enumerate(T):
        for a, b in ((0, 1), (1, 2), (2, 0)): E[frozenset((key(t[a]), key(t[b])))].append((i, a, b))
    bd = [(e, v[0]) for e, v in E.items() if len(v) == 1]
    adj = collections.defaultdict(list)
    for e, _ in bd:
        a, b = tuple(e); adj[a].append(b); adj[b].append(a)
    seen, loops = set(), []
    for s0 in adj:
        if s0 in seen: continue
        st, comp = [s0], []
        while st:
            u = st.pop()
            if u in seen: continue
            seen.add(u); comp.append(u); st += adj[u]
        loops.append(np.array(comp) / 1e9)
    discs, extra, info = {"upper": None, "lower": None, "outlet": None}, [], []
    for L in loops:
        zc = L[:, 2].mean()
        if abs(zc + 0.0295) < 1e-6:
            tri = Delaunay(L[:, :2]).simplices; P = L.copy()
            D = P[tri]; key_ = "upper" if L[:, 1].mean() > 0 else "lower"
            discs[key_] = D; info.append((key_, len(L)))
        elif abs(zc - 0.3395) < 1e-6:
            tri = Delaunay(L[:, :2]).simplices; discs["outlet"] = L[tri]; info.append(("outlet", len(L)))
        elif len(L) == 4 and abs(zc) < 1e-9:
            info.append(("flapHole", 4))
        elif len(L) == 3:
            extra.append(L); info.append(("smallHole", 3))
        else: info.append(("UNKNOWN", len(L)))
    # tiny features of the CAD triangulation (seam hole of 3 vertices, micro-edges): weld wall vertices closer than WELD (0.2 mm),
    # excluding vertices on the flap-hole loop, then drop collapsed triangles.
    from scipy.spatial import cKDTree
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse import coo_matrix
    WELD = float(os.environ.get("WELD", "2e-4"))
    pts = T.reshape(-1, 3)
    uniq, inv = np.unique(np.round(pts * 1e9).astype(np.int64), axis=0, return_inverse=True)
    U = uniq / 1e9
    protect = np.zeros(len(U), bool)
    for L in loops:
        if len(L) == 4: 
            for q in L: protect |= (np.abs(U - q).max(1) < 2e-9)
    pr = cKDTree(U).query_pairs(WELD, output_type="ndarray")
    pr = pr[~(protect[pr[:, 0]] | protect[pr[:, 1]])]
    g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(U), len(U)))
    nc, lab = connected_components(g, directed=False)
    R = np.zeros((nc, 3)); np.add.at(R, lab, U); R /= np.bincount(lab, minlength=nc)[:, None]
    for i in np.where(protect)[0]: R[lab[i]] = U[i]
    newid = lab[inv.reshape(-1)].reshape(-1, 3)
    keep = (newid[:, 0] != newid[:, 1]) & (newid[:, 1] != newid[:, 2]) & (newid[:, 0] != newid[:, 2])
    T = R[newid[keep]]
    ar = 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1)
    T = T[ar > 1e-12]
    # refill 3-vertex holes that remain after the weld
    E2 = collections.defaultdict(list)
    for t_ in T:
        for a_, b_ in ((0, 1), (1, 2), (2, 0)): E2[frozenset((key(t_[a_]), key(t_[b_])))].append(1)
    free = [tuple(e) for e, v in E2.items() if len(v) == 1]
    adj2 = collections.defaultdict(list)
    for a_, b_ in free: adj2[a_].append(b_); adj2[b_].append(a_)
    seen2, fl = set(), []
    for s_ in adj2:
        if s_ in seen2: continue
        st, comp = [s_], []
        while st:
            u = st.pop()
            if u in seen2: continue
            seen2.add(u); comp.append(u); st += adj2[u]
        if len(comp) == 3 and abs(np.array(comp)[:, 2].mean() / 1e9) > 1e-3 and abs(np.array(comp)[:, 2].mean() / 1e9 - 0.3395) > 1e-3: fl.append(np.array(comp) / 1e9)
    fills = np.array(fl).reshape(-1, 3, 3) if fl else np.zeros((0, 3, 3)); extra = []
    P = read()
    iface = P["interface"]
    a, b = 5.5e-3, 1.0e-3
    clamp = np.array([[[-a, -b, 0], [a, -b, 0], [a, b, 0]], [[-a, -b, 0], [a, b, 0], [-a, b, 0]]], float)
    patches = {"wall": np.concatenate([T, fills]), "upperInlet": discs["upper"], "lowerInlet": discs["lower"], "outlet": discs["outlet"], "interface": iface, "clamp": clamp}
    # global weld (1e-7 m) + exact planes, applied identically to all patches
    from scipy.spatial import cKDTree
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse import coo_matrix
    allp = np.concatenate([v.reshape(-1, 3) for v in patches.values()])
    for z0 in (-0.0295, 0.3395): allp[np.abs(allp[:, 2] - z0) < 5e-5, 2] = z0
    pr = cKDTree(allp).query_pairs(2e-7, output_type="ndarray")
    n = len(allp); g = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(n, n))
    nc, lab = connected_components(g, directed=False)
    rep_pts = np.zeros((nc, 3)); np.add.at(rep_pts, lab, allp); rep_pts /= np.bincount(lab, minlength=nc)[:, None]
    # keep exact flap-box corners: representatives within 1e-7 of a corner/plane take the exact coordinate
    ex = lambda x, vals: vals[np.abs(vals - x[:, None]).argmin(1)] if False else x
    allp = rep_pts[lab]
    for k, (name, v) in enumerate(list(patches.items())):
        m = len(v); patches[name] = allp[:3 * m].reshape(-1, 3, 3); allp = allp[3 * m:]
    allt = np.concatenate(list(patches.values()))
    area = lambda X: 0.5 * np.linalg.norm(np.cross(X[:, 1] - X[:, 0], X[:, 2] - X[:, 0]), axis=1)
    rep = dict(loops=info, areas_m2={k: float(area(v).sum()) for k, v in patches.items()}, ntri={k: int(len(v)) for k, v in patches.items()})
    # watertightness
    cnt = collections.Counter()
    for t in allt:
        for i, j in ((0, 1), (1, 2), (2, 0)):
            ka, kb = key(t[i]), key(t[j]); cnt[(min(ka, kb), max(ka, kb))] += 1
    rep["edge_usage"] = dict(collections.Counter(cnt.values()))
    return patches, rep
