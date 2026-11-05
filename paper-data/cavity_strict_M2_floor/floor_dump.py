import re, sys, numpy as np
case = sys.argv[1]; T = "200"
def foamlist(path, kind):
    s = open(path).read()
    s = s[s.index("FoamFile"):]; s = s[s.index("}")+1:]
    m = re.search(r"\n(\d+)\s*\n\(", s); n = int(m.group(1)); body = s[m.end():]
    if kind == "pts":
        a = np.array(re.findall(r"\(([^()]*)\)", body)[:n], dtype=object)
        return np.array([list(map(float, x.split())) for x in a])
    out = []
    for mm in re.finditer(r"(\d+)\(([^()]*)\)", body):
        out.append(list(map(int, mm.group(2).split())))
        if len(out) == n: break
    return out
pts = foamlist(f"{case}/{T}/fluid/polyMesh/points", "pts")
faces = foamlist(f"{case}/constant/fluid/polyMesh/faces", "faces")
b = open(f"{case}/constant/fluid/polyMesh/boundary").read()
m = re.search(r"interface\s*\{[^}]*nFaces\s+(\d+);[^}]*startFace\s+(\d+);", b)
nF, sF = int(m.group(1)), int(m.group(2))
Sf = []; Cf = []
for f in faces[sF:sF+nF]:
    P = pts[f]; c = P.mean(0); a = np.zeros(3)
    for i in range(len(P)): a += 0.5*np.cross(P[i]-c, P[(i+1)%len(P)]-c)
    Sf.append(a); Cf.append(c)
Sf = np.array(Sf); Cf = np.array(Cf)
def raw(name):
    return np.loadtxt(f"{case}/postProcessing/floorSamples/{T}/{name}_interface.raw")
p = raw("p"); wss = raw("wallShearStress"); U = raw("U")
assert np.allclose(p[:, :2], Cf[:, :2], atol=1e-9), "face order mismatch"
Fp = (p[:, 3:4]*Sf).sum(0)            # rho = 1 (kinematic p, rhoInf 1)
Fv = (-wss[:, 3:6]*np.linalg.norm(Sf, axis=1)[:, None]).sum(0)
print("pressure force from dump", Fp, "\nviscous force from dump", Fv)
mag = np.linalg.norm(Sf, axis=1); n = Sf/mag[:, None]
order = np.argsort(Cf[:, 0])
with open(sys.argv[2], "w") as o:
    o.write("# cavityFlexibleBottom strict partitioned level 2 (Aitken, requireAllResidualMeasures yes, dt 0.25), t = 200 s\n")
    o.write("# fluid interface (plate top) faces, deformed configuration; Sf points out of the fluid (into the plate); depth 1 m; rho = 1\n")
    o.write("# x y Sf_x Sf_y |Sf| p wallShearStress_x wallShearStress_y U_x U_y\n")
    for i in order:
        o.write(" ".join(f"{v:.9e}" for v in [Cf[i,0], Cf[i,1], Sf[i,0], Sf[i,1], mag[i], p[i,3], wss[i,3], wss[i,4], U[i,3], U[i,4]]) + "\n")
