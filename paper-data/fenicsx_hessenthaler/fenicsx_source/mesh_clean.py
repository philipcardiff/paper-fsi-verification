"""Collapse micro-edges (< THR) of a gmsh tet mesh by node merging (relabel), drop degenerate tets/facets, report quality before/after.
Nodes on the flap interface/clamp, inlets and outlet are preferred as representatives (they are never moved; merging two such nodes is refused).
usage: mesh_clean.py in.msh out.msh THR_m"""
import sys, json
import numpy as np
src, dst, thr = sys.argv[1], sys.argv[2], float(sys.argv[3])
z = np.load(src)
X = z["x"].copy(); tets = z["tets"].copy(); tris = z["tris"].copy(); ptet = z["tet_phys"].copy(); ptri = z["tri_phys"].copy()
def quality(X, T):
    P = X[T]; v = np.einsum("ij,ij->i", np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), P[:, 3] - P[:, 0]) / 6
    el = sum(np.sum((P[:, i] - P[:, j]) ** 2, 1) for i in range(4) for j in range(i + 1, 4))
    return 12 * (3 * np.abs(v)) ** (2 / 3) / el, v
q0, v0 = quality(X, tets)
prot = np.zeros(len(X), bool)
for t_, p_ in zip(tris, ptri):
    if p_ in (1, 2, 3, 5, 6): prot[t_] = True
rep = np.arange(len(X))
def find(a):
    while rep[a] != a: rep[a] = rep[rep[a]]; a = rep[a]
    return a
E = np.unique(np.sort(np.concatenate([tets[:, [i, j]] for i in range(4) for j in range(i + 1, 4)]), axis=1), axis=0)
L = np.linalg.norm(X[E[:, 0]] - X[E[:, 1]], axis=1)
short = E[L < thr]; nm = 0; refused = 0
for a, b in short:
    ra, rb = find(a), find(b)
    if ra == rb: continue
    if prot[ra] and prot[rb]: refused += 1; continue
    keep, drop = (ra, rb) if (prot[ra] or (not prot[rb] and ra < rb)) else (rb, ra)
    rep[drop] = keep; nm += 1
for k in range(len(X)): rep[k] = find(k)
tets2 = rep[tets]; tris2 = rep[tris]
ok = np.array([len(set(t)) == 4 for t in tets2]); tets2, ptet2 = tets2[ok], ptet[ok]
# second pass: remaining poor tets (mean-ratio < QTHR): collapse their shortest edge (interior/boundary rule as above), up to 30 sweeps
QTHR = 0.08; extra_merged = 0
for sweep in range(30):
    qq, vv = quality(X, tets2); bad = np.where(qq < QTHR)[0]
    if len(bad) == 0: break
    rep2 = np.arange(len(X)); fnd = lambda a: a
    changed = False
    for k in bad:
        t = tets2[k]; best = None
        for i in range(4):
            for j in range(i + 1, 4):
                a_, b_ = t[i], t[j]
                if prot[a_] and prot[b_]: continue
                L_ = np.linalg.norm(X[a_] - X[b_])
                if best is None or L_ < best[0]: best = (L_, a_, b_)
        if best is None: continue
        _, a_, b_ = best; keep, drop = (a_, b_) if (prot[a_] or (not prot[b_] and a_ < b_)) else (b_, a_)
        tets2 = np.where(tets2 == drop, keep, tets2); tris = np.where(tris == drop, keep, tris); tris2 = np.where(tris2 == drop, keep, tris2); extra_merged += 1
        ok2 = np.array([len(set(t_)) == 4 for t_ in tets2]); tets2, ptet2 = tets2[ok2], ptet2[ok2]; changed = True
        break      # one collapse per sweep per worst tet (indices change)
    if not changed: break
ok = np.ones(len(tets2), bool)
okf = np.array([len(set(t)) == 3 for t in tris2]); tris2, ptri2 = tris2[okf], ptri[okf]
q1, v1 = quality(X, tets2)
rep_info = dict(second_pass_collapses=int(extra_merged), n_short_edges=int(len(short)), merged=int(nm), refused=int(refused), dropped_tets=int((~ok).sum()), dropped_tris=int((~okf).sum()),
                min_q_before=float(q0.min()), min_q_after=float(q1.min()), n_q_lt_0p05_before=int((q0 < 0.05).sum()), n_q_lt_0p05_after=int((q1 < 0.05).sum()),
                n_inverted_after=int((np.sign(v1) != np.sign(np.median(v1))).sum()))
print(json.dumps(rep_info)); json.dump(rep_info, open(dst + ".clean.json", "w"), indent=1)
used = np.unique(np.concatenate([tets2.ravel(), tris2.ravel()])); idx = -np.ones(len(X), int); idx[used] = np.arange(len(used))
np.savez_compressed(dst, x=X[used], tets=idx[tets2], tet_phys=ptet2, tris=idx[tris2], tri_phys=ptri2)
