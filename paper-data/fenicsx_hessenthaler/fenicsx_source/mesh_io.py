"""Load the cleaned conforming Hessenthaler mesh (npz from geometry/scripts/mesh_clean.py) into dolfinx with cell/facet tags.
Cells: 1 fluid, 2 solid. Facets: 1 upperInlet 2 lowerInlet 3 outlet 4 wall 5 interface 6 clamp. Returns an object with .mesh, .cell_tags, .facet_tags."""
from types import SimpleNamespace
import numpy as np
import ufl, basix.ufl
from dolfinx import mesh as dmesh
from dolfinx.io.utils import distribute_entity_data
from dolfinx.graph import adjacencylist

def load(path, comm, rank=0):
    if comm.rank == rank:
        z = np.load(path); x = z["x"].astype(np.float64); tets = z["tets"].astype(np.int64); tp = z["tet_phys"].astype(np.int32)
        tris = z["tris"].astype(np.int64); fp = z["tri_phys"].astype(np.int32)
    else:
        x = np.empty((0, 3)); tets = np.empty((0, 4), np.int64); tp = np.empty(0, np.int32); tris = np.empty((0, 3), np.int64); fp = np.empty(0, np.int32)
    dom = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    msh = dmesh.create_mesh(comm, tets, dom, x, dmesh.create_cell_partitioner(dmesh.GhostMode.shared_facet))
    msh.topology.create_entities(2); msh.topology.create_connectivity(2, 3)
    out = {}
    for dim, ent, val in ((3, tets, tp), (2, tris, fp)):
        le, lv = distribute_entity_data(msh, dim, ent, val)
        msh.topology.create_connectivity(dim, 3)
        out[dim] = dmesh.meshtags_from_entities(msh, dim, adjacencylist(le), lv.astype(np.int32))
    return SimpleNamespace(mesh=msh, cell_tags=out[3], facet_tags=out[2])
