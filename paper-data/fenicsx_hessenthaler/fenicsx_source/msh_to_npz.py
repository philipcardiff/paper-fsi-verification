"""gmsh .msh -> raw npz (x, tets, tet_phys, tris, tri_phys) via the gmsh API (physical groups preserved). usage: msh_to_npz.py in.msh out.npz"""
import sys, numpy as np, gmsh
gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.open(sys.argv[1])
tags, xyz, _ = gmsh.model.mesh.getNodes(); xyz = xyz.reshape(-1, 3)
idx = {int(t): i for i, t in enumerate(tags)}
lut = np.full(int(tags.max()) + 1, -1, np.int64); lut[tags.astype(np.int64)] = np.arange(len(tags))
T, TP, F, FP = [], [], [], []
for dim, phys in gmsh.model.getPhysicalGroups():
    for ent in gmsh.model.getEntitiesForPhysicalGroup(dim, phys):
        et, el, en = gmsh.model.mesh.getElements(dim, ent)
        for t, e, n in zip(et, el, en):
            if dim == 3 and t == 4: a = lut[n.astype(np.int64)].reshape(-1, 4); T.append(a); TP.append(np.full(len(a), phys))
            if dim == 2 and t == 2: a = lut[n.astype(np.int64)].reshape(-1, 3); F.append(a); FP.append(np.full(len(a), phys))
gmsh.finalize()
np.savez_compressed(sys.argv[2], x=xyz, tets=np.concatenate(T), tet_phys=np.concatenate(TP), tris=np.concatenate(F), tri_phys=np.concatenate(FP))
print("tets", sum(len(a) for a in T), {int(p): int((np.concatenate(FP) == p).sum()) for p in np.unique(np.concatenate(FP))})
