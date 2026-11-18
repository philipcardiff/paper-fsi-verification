"""GATE A geometry/mesh audit of a Hessenthaler mesh (serial). usage: gate_a.py MESH.msh OUT.json"""
import sys, json
import numpy as np
from mpi4py import MPI
from dolfinx.io import gmsh as dgmsh
from dolfinx import mesh as dmesh
sys.path.insert(0, __import__("os").environ.get("HESS_ROOT", "."))
from src.hessenthaler import mesh_io
m = mesh_io.load(sys.argv[1], MPI.COMM_SELF, 0)
msh, ct, ft = m.mesh, m.cell_tags, m.facet_tags
X = msh.geometry.x; tdim = msh.topology.dim
msh.topology.create_connectivity(tdim - 1, 0); msh.topology.create_connectivity(tdim - 1, tdim); msh.topology.create_connectivity(tdim, 0)
c2v = msh.topology.connectivity(tdim, 0).array.reshape(-1, 4); f2v = msh.topology.connectivity(tdim - 1, 0)
f2c = msh.topology.connectivity(tdim - 1, tdim)
P = X[c2v]
vol6 = np.einsum("ij,ij->i", np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), P[:, 3] - P[:, 0]); vol = vol6 / 6
el = sum(np.sum((P[:, i] - P[:, j]) ** 2, 1) for i in range(4) for j in range(i + 1, 4))
q = 12 * (3 * np.abs(vol)) ** (2 / 3) / el           # mean-ratio quality, 1 = regular
res = dict(ncell=int(len(c2v)), npoints=int(len(X)), n_fluid=int((ct.values == 1).sum()), n_solid=int((ct.values == 2).sum()))
res["min_cell_volume_m3"] = float(np.abs(vol).min()); res["n_inverted_or_zero"] = int((np.abs(vol) < 1e-18).sum())
res["mean_ratio_quality_min"] = float(q.min()); res["mean_ratio_quality_p1"] = float(np.percentile(q, 1)); res["mean_ratio_quality_mean"] = float(q.mean())
res["vol_fluid_m3"] = float(np.abs(vol)[ct.indices[ct.values == 1]].sum()); res["vol_solid_m3"] = float(np.abs(vol)[ct.indices[ct.values == 2]].sum())
names = {1: "upperInlet", 2: "lowerInlet", 3: "outlet", 4: "wall", 5: "interface", 6: "clamp"}
fa = {}
for tag, nm in names.items():
    idx = ft.find(tag); v = np.array([f2v.links(f) for f in idx]); T = X[v]
    area = 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1)
    fa[nm] = dict(n=int(len(idx)), area_m2=float(area.sum()), bbox_min=T.reshape(-1, 3).min(0).tolist(), bbox_max=T.reshape(-1, 3).max(0).tolist())
    if nm in ("interface", "clamp"):
        # each facet touches exactly one fluid (interface) and one solid cell; clamp only solid
        nf = ns = 0
        for f in idx:
            cs = f2c.links(f); tg = set(ct.values[np.searchsorted(ct.indices, cs)]) if False else None
        cmap = np.zeros(msh.topology.index_map(tdim).size_local + msh.topology.index_map(tdim).num_ghosts, int); cmap[ct.indices] = ct.values
        cnt = {"fluid_only": 0, "solid_only": 0, "both": 0}
        for f in idx:
            tg = set(cmap[f2c.links(f)])
            cnt["both" if tg == {1, 2} else ("fluid_only" if tg == {1} else "solid_only")] += 1
        fa[nm]["cell_adjacency"] = cnt
        # coplanarity with flap faces
        Tn = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); Tn /= np.linalg.norm(Tn, axis=1)[:, None]
        fa[nm]["normal_axes_abs_max_offaxis"] = float(np.min(np.abs(Tn).max(1)) )
res["patches"] = fa
iv = np.unique(f2v.array.reshape(-1, 3)[ft.find(5)].ravel()); Pi = X[iv]
res["flap_interface_extent_m"] = dict(x=[float(Pi[:, 0].min()), float(Pi[:, 0].max())], y=[float(Pi[:, 1].min()), float(Pi[:, 1].max())], z=[float(Pi[:, 2].min()), float(Pi[:, 2].max())])
# interface vertices lie on the 5 box planes
dist = np.minimum.reduce([np.abs(np.abs(Pi[:, 0]) - 5.5e-3), np.abs(np.abs(Pi[:, 1]) - 1.0e-3), np.abs(Pi[:, 2] - 65e-3)])
res["max_interface_vertex_distance_to_nearest_box_plane_m"] = float(dist.max())
inl = {}
for nm, R in (("upperInlet", 10.95e-3), ("lowerInlet", 10.95e-3), ("outlet", 38.1e-3)):
    inl[nm] = dict(area_over_pi_R2=fa[nm]["area_m2"] / (np.pi * R ** 2))
res["circle_area_ratios"] = inl
res["flap_area_error_rel"] = fa["interface"]["area_m2"] / 1712e-6 - 1.0
res["solid_volume_error_rel"] = res["vol_solid_m3"] / (11e-3 * 2e-3 * 65e-3) - 1.0
import gmsh
z = np.load(sys.argv[1]); gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("q")
e_ = gmsh.model.addDiscreteEntity(3); gmsh.model.mesh.addNodes(3, e_, np.arange(1, len(z["x"]) + 1), z["x"].ravel()); gmsh.model.mesh.addElementsByType(e_, 4, [], (z["tets"] + 1).ravel())
sic = np.array(gmsh.model.mesh.getElementQualities(gmsh.model.mesh.getElementsByType(4)[0], "minSICN")); gmsh.finalize()
res["gmsh_minSICN"] = float(sic.min()); res["gmsh_n_SICN_lt_0p05"] = int((sic < 0.05).sum())
res["GATE_A"] = "PASS" if (abs(res["flap_area_error_rel"]) < 1e-3 and res["n_inverted_or_zero"] == 0 and res["gmsh_minSICN"] > 0.05 and
  fa["interface"]["cell_adjacency"]["both"] == fa["interface"]["n"] and abs(res["solid_volume_error_rel"]) < 1e-3 and
  all(abs(v["area_over_pi_R2"] - 1) < 0.005 for v in inl.values())) else "FAIL (see fields)"
json.dump(res, open(sys.argv[2], "w"), indent=1); print(json.dumps(res, indent=1))
