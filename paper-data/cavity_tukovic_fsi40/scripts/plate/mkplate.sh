#!/bin/bash
# Usage: mkplate.sh DIR NX NY NONLINEAR PLANESTRESS
# Plate-only check for the fsi_40 unsTotalLagrangianSolid: 4 m x 0.4 m,
# clamped ends, uniform follower pressure 1.315 Pa on the top face.
set -e
D=$1; NX=$2; NY=$3; NL=$4; PS=$5
H='FoamFile { version 2.0; format ascii; class dictionary; object'
mkdir -p $D/0 $D/constant/polyMesh $D/system
cat > $D/constant/polyMesh/blockMeshDict <<B
$H blockMeshDict; }
convertToMeters 1;
vertices ( (2 -0.4 0) (6 -0.4 0) (6 0 0) (2 0 0) (2 -0.4 0.05) (6 -0.4 0.05) (6 0 0.05) (2 0 0.05) );
blocks ( hex (0 1 2 3 4 5 6 7) ($NX $NY 1) simpleGrading (1 1 1) );
edges ();
patches (
  patch fixed ( (0 4 7 3) (1 2 6 5) )
  patch top ( (3 7 6 2) )
  patch bottom ( (0 1 5 4) )
  empty frontAndBack ( (0 3 2 1) (4 5 6 7) )
);
mergePatchPairs ();
B
cat > $D/constant/rheologyProperties <<B
$H rheologyProperties; }
planeStress $PS;
rheology { type linearElastic; rho rho [1 -3 0 0 0 0 0] 1000; E E [1 -1 -2 0 0 0 0] 500; nu nu [0 0 0 0 0 0 0] 0.4; }
B
cat > $D/constant/solidProperties <<B
$H solidProperties; }
solidSolver unsTotalLagrangianSolid;
unsTotalLagrangianSolidCoeffs
{ nCorrectors 100000; convergenceTolerance 1e-10; relConvergenceTolerance 0; nonLinear $NL; debug no; }
B
cat > $D/0/D <<B
$H D; }
dimensions [0 1 0 0 0 0 0];
internalField uniform (0 0 0);
boundaryField {
  fixed { type fixedDisplacement; value uniform (0 0 0); }
  top { type tractionDisplacement; traction uniform (0 0 0); pressure uniform 1.315; value uniform (0 0 0); }
  bottom { type tractionDisplacement; traction uniform (0 0 0); pressure uniform 0; value uniform (0 0 0); }
  frontAndBack { type empty; }
}
B
sed -i 's/class dictionary; object D;/class volVectorField; object D;/' $D/0/D
cat > $D/0/pointD <<B
FoamFile { version 2.0; format ascii; class pointVectorField; object pointD; }
dimensions [0 1 0 0 0 0 0];
internalField uniform (0 0 0);
boundaryField {
  fixed { type fixedValue; value uniform (0 0 0); }
  top { type calculated; value uniform (0 0 0); }
  bottom { type calculated; value uniform (0 0 0); }
  frontAndBack { type empty; }
}
B
cat > $D/system/controlDict <<B
$H controlDict; }
application solidFoam; startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 1;
writeControl timeStep; writeInterval 1; writeFormat ascii; writePrecision 10; writeCompression uncompressed;
timeFormat general; timePrecision 6; runTimeModifiable no;
B
cat > $D/system/fvSchemes <<B
$H fvSchemes; }
d2dt2Schemes { default steadyState; }
ddtSchemes { default steadyState; }
gradSchemes { default none; }
divSchemes { default none; }
laplacianSchemes { default none; laplacian(DD,D) Gauss linear skewCorrected 1; }
snGradSchemes { default none; snGrad(D) skewCorrected 1; }
interpolationSchemes { default none; interpolate(mu) linear; interpolate(lambda) linear; }
B
cat > $D/system/fvSolution <<B
$H fvSolution; }
solvers { D { solver PCG; preconditioner DIC; tolerance 1e-12; relTol 0.01; } }
relaxationFactors { fields { D 0.7; } equations { } }
B
