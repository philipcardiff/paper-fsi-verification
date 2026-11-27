"""Build a foam-extend-4.0 fsiFoam (Fsi_40) case of the flexible-bottom cavity.
Usage: python3 -I mkcase.py OUTDIR R DT ENDTIME
R: refinement factor relative to the solids4foam tutorial mesh (L1: R=1, dx=0.1).
Fluid mesh is the solids4foam tutorial blockMeshDict with all cell counts x R;
solid is (40R x 4R). StVK plate via unsTotalLagrangianSolid nonLinear totalLagrangian."""
import os, re, sys
here = os.path.dirname(os.path.abspath(__file__))
out, R, DT, END = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
K = sys.argv[5] if len(sys.argv) > 5 else "0.1"
HDR = "FoamFile { version 2.0; format ascii; class %s; object %s; }\n"

def w(path, cls, obj, body):
    p = os.path.join(out, path)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    open(p, "w").write(HDR % (cls, obj) + body)

# ---------- fluid mesh
src = open(os.path.join(here, "s4f_fluid_blockMeshDict")).read()
src = re.sub(r"/\*.*?\*/", "", src, flags=re.S)
src = re.sub(r"//[^\n]*", "", src)
verts = re.search(r"vertices\s*\((.*?)\n\);", src, re.S).group(1)
blocks = re.search(r"blocks\s*\((.*?)\n\);", src, re.S).group(1)
blocks = re.sub(r"(hex \([^)]*\)) \((\d+) (\d+) 1\)", lambda m: "%s (%d %d 1)" % (m.group(1), int(m.group(2))*R, int(m.group(3))*R), blocks)
edges = re.search(r"edges\s*\((.*?)\n\);", src, re.S).group(1)
bnd = re.search(r"boundary\s*\((.*)\n\);\s*mergePatchPairs", src, re.S).group(1)
patches = []
for name, typ, faces in re.findall(r"(\w+)\s*\{\s*type\s+(\w+);\s*faces\s*\((.*?)\);\s*\}", bnd, re.S):
    patches.append("    %s %s\n    (%s)\n" % (typ, name, faces))
w("fluid/constant/polyMesh/blockMeshDict", "dictionary", "blockMeshDict",
  "convertToMeters 1;\nvertices\n(%s\n);\nblocks\n(%s\n);\nedges\n(%s\n);\npatches\n(\n%s);\nmergePatchPairs ();\n"
  % (verts, blocks, edges, "".join(patches)))

w("fluid/0/U", "volVectorField", "U", """dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField {
  inlet { type transitionalParabolicVelocity; maxValue 1.5; n (1 0 0); y (0 1 0);
          transitionPeriod 10; boundBoxMin (0 1 0); boundBoxMax (0 2 1); value uniform (0 0 0); }
  outlet { type zeroGradient; }
  interface { type movingWallVelocity; value uniform (0 0 0); }
  walls { type fixedValue; value uniform (0 0 0); }
  frontAndBack { type empty; }
}
""")
w("fluid/0/p", "volScalarField", "p", """dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField {
  inlet { type extrapolatedPressure; value uniform 0; }
  outlet { type fixedValue; value uniform 0; }
  interface { type extrapolatedPressure; value uniform 0; }
  walls { type extrapolatedPressure; value uniform 0; }
  frontAndBack { type empty; }
}
""")
w("fluid/0/pointMotionU", "pointVectorField", "pointMotionU", """dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField {
  inlet { type fixedValue; value uniform (0 0 0); }
  outlet { type fixedValue; value uniform (0 0 0); }
  interface { type fixedValue; value uniform (0 0 0); }
  walls { type fixedValue; value uniform (0 0 0); }
  frontAndBack { type empty; }
}
""")
w("fluid/constant/fsiProperties", "dictionary", "fsiProperties", """solidPatch interface;
solidZone interfaceZone;
fluidPatch interface;
fluidZone interfaceZone;
relaxationFactor 0.1;
interfaceDeformationLimit 0;
outerCorrTolerance 1e-5;
nOuterCorr 30;
interpolatorUpdateFrequency 0;
couplingScheme Aitken;
couplingReuse 0;
coupled yes;
predictor yes;
""")
w("fluid/constant/fluidProperties", "dictionary", "fluidProperties", """fluidSolver consistentIcoFluid;
icoFluidCoeffs { nCorrectors 3; nNonOrthogonalCorrectors 1; nOuterCorrectors 1; }
consistentIcoFluidCoeffs { nCorrectors 3; nNonOrthogonalCorrectors 1; }
""")
w("fluid/constant/dynamicMeshDict", "dictionary", "dynamicMeshDict", """dynamicFvMesh dynamicMotionSolverFvMesh;
solver refVelocityLaplacian;
diffusivity quadratic inverseDistance 1(interface);
nNonOrthogonalCorrectors 2;
leastSquaresVolPoint yes;
""")
w("fluid/constant/transportProperties", "dictionary", "transportProperties",
  "nu nu [0 2 -1 0 0 0 0] 0.01;\nrho rho [1 -3 0 0 0 0 0] 1;\n")
w("fluid/system/controlDict", "dictionary", "controlDict", """application fsiFoam;
startFrom startTime; startTime 0; stopAt endTime; endTime %s; deltaT %s;
writeControl runTime; writeInterval 25; purgeWrite 0;
writeFormat ascii; writePrecision 10; writeCompression uncompressed;
timeFormat general; timePrecision 6; runTimeModifiable yes; adjustTimeStep no; maxCo 10;
""" % (END, DT))
w("fluid/system/fvSchemes", "dictionary", "fvSchemes", """ddtSchemes { default Euler; }
gradSchemes { default leastSquares; }
divSchemes { default none; div(phi,U) Gauss linearUpwind cellLimited leastSquares 1; }
laplacianSchemes { default none; laplacian(nu,U) Gauss linear skewCorrected 1;
  laplacian((1|A(U)),p) Gauss linear skewCorrected 1; laplacian(diffusivity,cellMotionU) Gauss linear skewCorrected 1; }
interpolationSchemes { default linear; interpolate(U) skewCorrected linear; }
snGradSchemes { default skewCorrected 1; }
fluxRequired { default no; p; }
""")
gamg = """solver GAMG; tolerance 1e-08; relTol 0; minIter 1; maxIter 1000; smoother GaussSeidel;
 nPreSweeps 0; nPostSweeps 2; nFinestSweeps 2; scaleCorrection true; directSolveCoarsest false;
 cacheAgglomeration true; nCellsInCoarsestLevel 20; agglomerator faceAreaPair; mergeLevels 1;"""
w("fluid/system/fvSolution", "dictionary", "fvSolution", """solvers {
 p { %s }
 cellMotionU { %s }
 U { solver PBiCG; preconditioner DILU; tolerance 1e-08; relTol 0; minIter 1; }
}
""" % (gamg, gamg))
w("fluid/setBatch", "dictionary", "x", "")
open(os.path.join(out, "fluid/setBatch"), "w").write("faceSet interfaceZone new patchToFace interface\nquit\n")

# ---------- solid
nx, ny = 40*R, 4*R
w("solid/constant/polyMesh/blockMeshDict", "dictionary", "blockMeshDict", """convertToMeters 1;
vertices ( (2 -0.4 0) (6 -0.4 0) (6 0 0) (2 0 0) (2 -0.4 1) (6 -0.4 1) (6 0 1) (2 0 1) );
blocks ( hex (0 1 2 3 4 5 6 7) (%d %d 1) simpleGrading (1 1 1) );
edges ();
patches (
  patch interface ( (3 7 6 2) )
  patch fixed ( (0 4 7 3) (2 6 5 1) )
  patch bottom ( (4 0 1 5) )
  empty frontAndBack ( (4 5 6 7) (3 2 1 0) )
);
mergePatchPairs ();
""" % (nx, ny))
w("solid/constant/rheologyProperties", "dictionary", "rheologyProperties", """planeStress no;
rheology { type linearElastic; rho rho [1 -3 0 0 0 0 0] 1000; E E [1 -1 -2 0 0 0 0] 500; nu nu [0 0 0 0 0 0 0] 0.4; }
""")
w("solid/constant/solidProperties", "dictionary", "solidProperties", """solidSolver unsTotalLagrangianSolid;
unsTotalLagrangianSolidCoeffs
{ nCorrectors 1000; convergenceTolerance 1e-7; relConvergenceTolerance 1e-3; nonLinear totalLagrangian; debug no; K K [0 0 -1 0 0 0 0] %s; }
""" % K)
w("solid/0/D", "volVectorField", "D", """dimensions [0 1 0 0 0 0 0];
internalField uniform (0 0 0);
boundaryField {
  interface { type tractionDisplacement; traction uniform (0 0 0); pressure uniform 0; value uniform (0 0 0); }
  fixed { type fixedDisplacement; value uniform (0 0 0); }
  bottom { type tractionDisplacement; traction uniform (0 0 0); pressure uniform 0; value uniform (0 0 0); }
  frontAndBack { type empty; }
}
""")
w("solid/0/pointD", "pointVectorField", "pointD", """dimensions [0 1 0 0 0 0 0];
internalField uniform (0 0 0);
boundaryField {
  interface { type calculated; value uniform (0 0 0); }
  fixed { type fixedValue; value uniform (0 0 0); }
  bottom { type calculated; value uniform (0 0 0); }
  frontAndBack { type empty; }
}
""")
w("solid/system/controlDict", "dictionary", "controlDict", "application fsiFoam; startFrom startTime; startTime 0; stopAt endTime; endTime %s; deltaT %s; writeControl runTime; writeInterval 25;\n" % (END, DT))
w("solid/system/fvSchemes", "dictionary", "fvSchemes", """d2dt2Schemes { default Euler; }
ddtSchemes { default Euler; }
gradSchemes { default none; }
divSchemes { default none; }
laplacianSchemes { default none; laplacian(DD,D) Gauss linear skewCorrected 1; }
snGradSchemes { default none; snGrad(D) skewCorrected 1; }
interpolationSchemes { default none; interpolate(mu) linear; interpolate(lambda) linear; }
""")
w("solid/system/fvSolution", "dictionary", "fvSolution", """solvers { D { solver PCG; preconditioner DIC; tolerance 1e-09; relTol 0.1; } }
relaxationFactors { fields { D 1; } equations { } }
""")
open(os.path.join(out, "solid/setBatch"), "w").write("faceSet interfaceZone new patchToFace interface\nquit\n")
open(os.path.join(out, "Allrun"), "w").write("""#!/bin/bash
source ~/tukovic-fsi40/env-fe40.sh
cd solid && blockMesh > log.blockMesh 2>&1 && setSet -batch setBatch > log.setSet 2>&1 && setsToZones -noFlipMap > log.setsToZones 2>&1 && cd ..
cd fluid && blockMesh > log.blockMesh 2>&1 && setSet -batch setBatch > log.setSet 2>&1 && setsToZones -noFlipMap > log.setsToZones 2>&1
(cd constant && ln -sf ../../solid/constant solid) ; (cd system && ln -sf ../../solid/system solid) ; (cd 0 && ln -sf ../../solid/0 solid)
fsiFoam > log.fsiFoam 2>&1
echo fin > ../done
""")
os.chmod(os.path.join(out, "Allrun"), 0o755)

# ---------- optional: incremental total Lagrangian solid (6th arg "incr")
if len(sys.argv) > 6 and sys.argv[6] == "incr":
    w("solid/constant/solidProperties", "dictionary", "solidProperties", """solidSolver unsIncrTotalLagrangianSolid;
unsIncrTotalLagrangianSolidCoeffs
{ nCorrectors 1000; convergenceTolerance 1e-7; relConvergenceTolerance 1e-3; nonLinear totalLagrangian; debug no; K K [0 0 -1 0 0 0 0] %s; }
""" % K)
    w("solid/0/DD", "volVectorField", "DD", """dimensions [0 1 0 0 0 0 0];
internalField uniform (0 0 0);
boundaryField {
  interface { type tractionDisplacementIncrement; traction uniform (0 0 0); pressure uniform 0; value uniform (0 0 0); }
  fixed { type fixedDisplacement; value uniform (0 0 0); }
  bottom { type tractionDisplacementIncrement; traction uniform (0 0 0); pressure uniform 0; value uniform (0 0 0); }
  frontAndBack { type empty; }
}
""")
    s = open(os.path.join(out, "solid/0/pointD")).read().replace("object pointD", "object pointDD")
    open(os.path.join(out, "solid/0/pointDD"), "w").write(s)
    os.remove(os.path.join(out, "solid/0/D"))
