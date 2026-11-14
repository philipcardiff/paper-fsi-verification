"""Build a conforming single-mesh (fluid + solid tets) for the Hessenthaler Phase I benchmark with gmsh.
Input: s4f tutorial fluidDomain.stl.gz (published surface, metres; patches upperInlet, lowerInlet, outlet, wall, interface).
Usage: build_mesh.py LEVEL OUTDIR   with LEVEL in {Gc,G0,G1,G2}  (all target sizes scale by r^k)
Tags (cells): 1 fluid, 2 solid.  Facets: 1 upperInlet, 2 lowerInlet, 3 outlet, 4 wall, 5 interface, 6 clamp (flap base).
No chamfer: the interface is the exact 11x2x65 mm box (5 faces); the flap base rectangle at z=0 is the clamp surface.
"""
import sys, gzip, re, os, json, hashlib
import numpy as np
import gmsh

STL = __import__("os").environ.get("HESS_ROOT", ".") + "/docs/s4f_ref_sources/tut/constant/triSurface/fluidDomain.stl.gz"
LEVELS = {"Gc": 1.5, "Gm": 1.22, "G0": 1.0, "G1": 1.0 / 1.5, "G2": 1.0 / 2.25}
BASE = dict(flap=1.0e-3, sweep=2.0e-3, wall=4.0e-3, far=10.0e-3, inlet=1.5e-3, solid=1.0e-3)

def read_patches():
    txt = gzip.open(STL, "rt").read()
    P = {}
    for s in re.split(r"\nendsolid[^\n]*\n?", txt):
        m = re.match(r"\s*solid (\S+)", s)
        if m:
            P[m.group(1)] = np.array(re.findall(r"vertex\s+(\S+)\s+(\S+)\s+(\S+)", s), float).reshape(-1, 3, 3)
    return P

def write_stl(path, tris):
    with open(path, "w") as f:
        f.write("solid all\n")
        for t in tris:
            n = np.cross(t[1] - t[0], t[2] - t[0]); nn = np.linalg.norm(n)
            n = n / nn if nn > 0 else n
            f.write("facet normal %.17g %.17g %.17g\nouter loop\n" % tuple(n))
            for v in t: f.write("vertex %.17g %.17g %.17g\n" % tuple(v))
            f.write("endloop\nendfacet\n")
        f.write("endsolid all\n")

def main(level, out):
    os.makedirs(out, exist_ok=True)
    s = LEVELS[level]; h = {k: v * s for k, v in BASE.items()}
    if os.environ.get('INLET_FIXED'): h['inlet'] = min(h['inlet'], BASE['inlet'])   # keep the inlet circumference resolution of level G0 (Gate A inlet-area criterion)
    # inlet circumference is resolved no coarser than 1.5 mm *s (see GATE A inlet-area check)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import decimate_surface as D
    D.main(0.985, os.path.join(out, "wall_dec"))
    patches, rep = D.assemble(os.path.join(out, "wall_dec.npy"), None)
    json.dump(rep, open(os.path.join(out, "surface_assembly.json"), "w"), indent=1)
    allt = np.concatenate(list(patches.values()))
    stl = os.path.join(out, "surface_all.stl"); write_stl(stl, allt)
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 1); gmsh.option.setNumber("General.NumThreads", 1)
    gmsh.model.add("hess")
    gmsh.merge(stl)
    ang = float(os.environ.get('ANGLE', '40')) * np.pi / 180
    gmsh.model.mesh.classifySurfaces(ang, True, True, np.pi)
    gmsh.model.mesh.createGeometry()
    surfs = gmsh.model.getEntities(2)
    print("n surfaces", len(surfs))
    # classify each surface by its centroid / bounding box
    cls = {}
    for (d, t) in surfs:
        bb = np.array(gmsh.model.getBoundingBox(d, t)).reshape(2, 3)
        lo, hi = bb[0], bb[1]; c = 0.5 * (lo + hi)
        tol = 1e-6
        if abs(hi[2] - lo[2]) < tol and abs(c[2] + 0.0295) < tol: tag = 1 if c[1] > 0 else 2
        elif abs(hi[2] - lo[2]) < tol and abs(c[2] - 0.3395) < tol: tag = 3
        elif abs(hi[2] - lo[2]) < tol and abs(c[2]) < tol and hi[0] <= 5.5e-3 + tol and lo[0] >= -5.5e-3 - tol and hi[1] <= 1e-3 + tol and lo[1] >= -1e-3 - tol: tag = 6
        elif (lo[0] >= -5.5e-3 - tol and hi[0] <= 5.5e-3 + tol and lo[1] >= -1e-3 - tol and hi[1] <= 1e-3 + tol and lo[2] >= -tol and hi[2] <= 0.065 + tol): tag = 5
        else: tag = 4
        cls[t] = tag
    for t, g in sorted(cls.items()): print("surface", t, "->", g)
    byTag = {g: [t for t, gg in cls.items() if gg == g] for g in range(1, 7)}
    assert len(byTag[5]) >= 5 and len(byTag[6]) >= 1 and byTag[1] and byTag[2] and byTag[3], byTag
    # volumes: solid = interface + clamp ; fluid = wall + inlets + outlet + interface
    sl_s = gmsh.model.geo.addSurfaceLoop(byTag[5] + byTag[6]) if False else None
    sl_s = gmsh.model.geo.addSurfaceLoop(byTag[5] + byTag[6])
    vs = gmsh.model.geo.addVolume([sl_s])
    sl_f = gmsh.model.geo.addSurfaceLoop(byTag[1] + byTag[2] + byTag[3] + byTag[4] + byTag[5])
    vf = gmsh.model.geo.addVolume([sl_f])
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(3, [vf], 1, "fluid"); gmsh.model.addPhysicalGroup(3, [vs], 2, "solid")
    names = {1: "upperInlet", 2: "lowerInlet", 3: "outlet", 4: "wall", 5: "interface", 6: "clamp"}
    for g in range(1, 7): gmsh.model.addPhysicalGroup(2, byTag[g], g, names[g])
    # --- size fields
    F = gmsh.model.mesh.field
    iface_s = byTag[5]
    F.add("Distance", 1); F.setNumbers(1, "SurfacesList", iface_s); F.setNumber(1, "Sampling", 200)
    F.add("Threshold", 2); F.setNumber(2, "InField", 1); F.setNumber(2, "SizeMin", h["flap"]); F.setNumber(2, "SizeMax", h["far"]); F.setNumber(2, "DistMin", 0.0); F.setNumber(2, "DistMax", 0.06)
    F.add("Box", 3); F.setNumber(3, "VIn", h["sweep"]); F.setNumber(3, "VOut", h["far"])
    F.setNumber(3, "XMin", -0.012); F.setNumber(3, "XMax", 0.012); F.setNumber(3, "YMin", -0.006); F.setNumber(3, "YMax", 0.036)
    F.setNumber(3, "ZMin", -0.004); F.setNumber(3, "ZMax", 0.075); F.setNumber(3, "Thickness", 0.02)
    F.add("Box", 4); F.setNumber(4, "VIn", h["wall"]); F.setNumber(4, "VOut", h["far"])
    F.setNumber(4, "XMin", -0.04); F.setNumber(4, "XMax", 0.04); F.setNumber(4, "YMin", -0.04); F.setNumber(4, "YMax", 0.04)
    F.setNumber(4, "ZMin", -0.0295); F.setNumber(4, "ZMax", 0.11); F.setNumber(4, "Thickness", 0.06)
    # inlet pipes: distance from inlet surfaces
    F.add("Distance", 5); F.setNumbers(5, "SurfacesList", byTag[1] + byTag[2]); F.setNumber(5, "Sampling", 100)
    F.add("Threshold", 6); F.setNumber(6, "InField", 5); F.setNumber(6, "SizeMin", h["inlet"]); F.setNumber(6, "SizeMax", h["wall"]); F.setNumber(6, "DistMin", 0.004); F.setNumber(6, "DistMax", 0.03)
    F.add("Min", 7); F.setNumbers(7, "FieldsList", [2, 3, 4, 6]); F.setAsBackgroundMesh(7)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0); gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0); gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    # solid is meshed with the same flap size (interface is shared); interior of solid: limit to h['solid']
    F.add("Constant", 8); F.setNumber(8, "VIn", h["solid"]); F.setNumbers(8, "VolumesList", [vs]) ; F.setNumber(8, "VOut", 1.0)
    F.add("Min", 9); F.setNumbers(9, "FieldsList", [2, 3, 4, 6, 8]); F.setAsBackgroundMesh(9)
    gmsh.option.setNumber("Mesh.Algorithm", int(os.environ.get("ALG2D", "6"))); gmsh.option.setNumber("Mesh.Algorithm3D", int(os.environ.get("ALG3D", "1")))
    gmsh.option.setNumber("Mesh.Optimize", 1); gmsh.option.setNumber("Mesh.OptimizeNetgen", 0)
    gmsh.model.mesh.generate(3)
    def _minq():
        et, tg, _ = gmsh.model.mesh.getElements(3)
        return np.array(gmsh.model.mesh.getElementQualities(np.concatenate(tg), "minSICN"))
    for it in range(10):
        for meth in ("Netgen", "Relocate3D", "", "Netgen"):
            gmsh.model.mesh.optimize(meth, True)
        q_ = _minq(); print("opt iter", it, "min SICN", q_.min(), "n<0.05", int((q_ < 0.05).sum()))
        if q_.min() >= 0.05: break
    et, tg, _ = gmsh.model.mesh.getElements(3)
    allt_ = np.concatenate(tg)
    sicn = np.array(gmsh.model.mesh.getElementQualities(allt_, "minSICN"))
    json.dump(dict(n_tets=int(len(sicn)), min_SICN=float(sicn.min()), p1_SICN=float(np.percentile(sicn, 1)), mean_SICN=float(sicn.mean()), n_SICN_lt_0p05=int((sicn < 0.05).sum())),
              open(os.path.join(out, "gmsh_quality.json"), "w"), indent=1)
    print("minSICN", sicn.min(), "p1", np.percentile(sicn, 1))
    gmsh.model.mesh.removeDuplicateNodes()
    gmsh.write(os.path.join(out, f"hess_{level}.msh"))
    types, tags, _ = gmsh.model.mesh.getElements(3)
    print("tets", sum(len(t) for t in tags))
    gmsh.finalize()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
