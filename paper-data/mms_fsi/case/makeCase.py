#!/usr/bin/env python3
"""
Write a solids4foam case for the manufactured FSI solution (../derive).

Geometry (2-D, one cell thick):
  fluid  [0, 2] x [0, 1]      patches inlet, outlet, walls (y = 0), interface (y = 1)
  solid  [0, 2] x [1, 1.5]    patches interface (y = 1), top, left, right (all clamped)

Usage:
  makeCase.py <dir> --coupling qm|iqnils|aitken|robin|fluid|solid --nx 32 --dt 0.01 --end 0.4

  qm      newtonQuasiMonolithic: newtonIcoFluid + velocity TL solid (needs solids4foam
          with #751, feature-petsc-snes-quasi-monolithic db3039b76 or later)
  aitken, robin, iqnils
          partitioned: pimpleFluid + displacement TL solid
  fluid   fluid alone on the fluid mesh, with the exact interface motion and
          velocity prescribed (--fluidModel newtonIcoFluid | pimpleFluid)
  solid   solid alone, with the exact interface traction applied

The solids4foam build is chosen with --env (an environment file that loads
OpenFOAM and sets WM_PROJECT_USER_DIR to the solids4foam build); the default
is $MMS_FSI_ENV.

nx is the number of cells along the interface; the fluid has nx x nx/2 cells
and the solid nx x nx/4 (square cells, conformal interface).
"""
import argparse
import os
import textwrap

ap = argparse.ArgumentParser()
ap.add_argument("dir")
ap.add_argument("--coupling", default="qm", choices=["qm", "iqnils", "aitken", "robin", "fluid", "solid"])
ap.add_argument("--nx", type=int, default=32)
ap.add_argument("--dt", type=float, default=0.01)
ap.add_argument("--end", type=float, default=0.4)
ap.add_argument("--motion", default="displacementLaplacian",
                choices=["displacementLaplacian", "displacementSBRStress", "velocityLaplacian"])
ap.add_argument("--diffusivity", default="quadratic inverseDistance (interface)")
ap.add_argument("--snesRtol", default="1e-10")
ap.add_argument("--pStab", default="1.0", help="QM Rhie-Chow scale factor")
ap.add_argument("--kspRtol", default="1e-6")
ap.add_argument("--relax", default="0.05")
ap.add_argument("--outerTol", default="1e-7", help="partitioned FSI outer tolerance")
ap.add_argument("--wallP", default="fixedGradient",
                help="QM bottom-wall p BC: fixedGradient (momentum balance) or zeroGradient")
ap.add_argument("--solidStab", default="0.1")
ap.add_argument("--variant", default="full", choices=["full", "flat", "mu10"],
                help="manufactured solution variant (flat: flat interface, orthogonal ALE mesh)")
ap.add_argument("--uCorrected", action="store_true",
                help="QM: fixedValueCorrected U (boundary non-orthogonal correction of snGrad) on interface and wall")
ap.add_argument("--pOracle", default=None, choices=["exact", "correction"],
                help="QM diagnostic: interface p value replaced by the exact value, or given the exact k.grad(p)")
ap.add_argument("--env", default=os.environ.get("MMS_FSI_ENV"),
                help="solids4foam build environment sourced by Allrun (default $MMS_FSI_ENV)")
ap.add_argument("--rcHistory", action="store_true",
                help="QM: RhieChowTimeHistory yes (time-history-consistent Rhie-Chow)")
ap.add_argument("--pCorrected", action="store_true",
                help="QM: fixedGradientCorrected (boundary non-orthogonal correction) on interface and wall p")
ap.add_argument("--interfaceP", default="fixedFluxPressure",
                help="pimpleFluid interface p BC: fixedFluxPressure, zeroGradient or movingWallPressure")
ap.add_argument("--writeInterval", type=float, default=None)
ap.add_argument("--grad", default="leastSquares", help="fluid gradSchemes default (QM and fluid-only)")
ap.add_argument("--interp", default="linear", help="QM fluid interpolationSchemes default")
ap.add_argument("--solidSurfGrad", action="store_true",
                help="velocity solid: tractionPatchSurfaceGradient yes (solids4foam #756)")
ap.add_argument("--np", type=int, default=1, help="number of MPI ranks (simple decomposition in x)")
ap.add_argument("--fluidModel", default="newtonIcoFluid", choices=["newtonIcoFluid", "pimpleFluid"],
                help="fluid-only mode (--coupling fluid): the fluid model")
ap.add_argument("--outlet", default="neumannU", choices=["neumannU", "dirichletU"],
                help="fluid outlet: exact Neumann U + exact Dirichlet p (default), or exact Dirichlet U and p "
                     "(no Neumann-U/Dirichlet-p outlet next to the wall)")
ap.add_argument("--ddtPhiCoeff", default=None, help="pimpleFluid: ddt scheme 'backward <ddtPhiCoeff>'")
ap.add_argument("--pimpleGrad", default="Gauss linear", help="pimpleFluid: gradSchemes default")
ap.add_argument("--ddtCorr", default=None, choices=["yes", "no"], help="pimpleFluid: PIMPLE ddtCorr")
ap.add_argument("--fsiExtra", default="", help="partitioned: extra entries for the coupling coeffs")
ap.add_argument("--solidModel", default="velocity", choices=["velocity", "totalDisp"],
                help="solid-only mode: velocity (QM solid) or totalDisp (displacement-form TL solid)")
ap.add_argument("--nominalTraction", action="store_true",
                help="solid-only mode: exact nominal traction with useUndeformedArea yes")
ap.add_argument("--solidExtra", default="",
                help="extra entries for the solid model coeffs dictionary")
ap.add_argument("--tractionExtra", default="",
                help="solid-only mode: extra entries for the solid interface traction condition")
a = ap.parse_args()
if not a.env:
    ap.error("give --env or set MMS_FSI_ENV")
a.env = os.path.abspath(a.env)
SSG = "\n    tractionPatchSurfaceGradient yes;" if a.solidSurfGrad else ""
RC = "\n    RhieChowTimeHistory yes;" if a.rcHistory else ""


# pimpleFluid dictionaries (partitioned couplings and the pimpleFluid fluid-only mode)
DDTCORR = ("\n    ddtCorr " + a.ddtCorr + ";") if a.ddtCorr else ""
PIMPLE_SCHEMES = """
ddtSchemes { default backward""" + (" " + a.ddtPhiCoeff if a.ddtPhiCoeff else "") + """; }
gradSchemes { default """ + a.pimpleGrad + """; }
divSchemes
{
    default            none;
    div(phi,U)         Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
    div((nuEff*dev(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
fluxRequired { default no; p; }
"""
PIMPLE_SOLUTION = f"""
solvers
{{
    "p|pFinal|pcorr|pcorrFinal"
    {{
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-12;
        relTol          0;
    }}
    "(cellMotionU|cellDisplacement)(Final)?"
    {{
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-12;
        relTol          0;
    }}
    "U|UFinal"
    {{
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-12;
        relTol          0;
    }}
}}
PIMPLE
{{
    momentumPredictor   on;
    nOuterCorrectors    1;
    nCorrectors         4;
    nNonOrthogonalCorrectors 1;{DDTCORR}
}}
"""

# Fluid outlet velocity: exact Neumann (default) or exact Dirichlet
OUTLET_U = ("type mmsFsiVelocity;\nvalue uniform (0 0 0);" if a.outlet == "dirichletU"
            else "type mmsFsiVelocityGradient;\nvalue uniform (0 0 0);")


def endMotion(variant):
    """Inlet and outlet point motion of the fluid mesh. Flat variant: the exact
    linear stretching u_y = eta(t) y0/H (a slip condition lets the corner point
    lag the interface and folds the corner cell); y0 = y - (last displacement)."""
    if variant != "flat":
        return "type fixedValue;\nvalue uniform (0 0 0);"
    return """type codedFixedValue;
value uniform (0 0 0);
name flatEndMotion;
codeInclude #{ #include "mmsFsiExact.H" #};
codeOptions #{ -I%s #};
code
#{
    const scalar t = this->db().time().value();
    double h[1];
    mmsFsi::channelHeight(0.0, t, h);
    const scalar eta = h[0] - mmsFsi::param_H;
    const pointField& pts = this->patch().localPoints();
    vectorField disp(pts.size());
    forAll(pts, i)
    {
        const scalar y0 = pts[i].y() - (*this)[i].y();
        disp[i] = vector(0, eta*y0/mmsFsi::param_H, 0);
    }
    this->operator==(disp);
#};""" % os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "derive", "generated_flat"))

QM = a.coupling == "qm"
nx, nyf, nys = a.nx, a.nx // 2, a.nx // 4
assert nx % 4 == 0

HDR = """/*--------------------------------*- C++ -*----------------------------------*\\
  Manufactured-solution FSI case (generated by makeCase.py)
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def w(rel, body, cls="dictionary", obj=None):
    path = os.path.join(a.dir, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(HDR.format(cls=cls, obj=obj or os.path.basename(rel)))
        fh.write(textwrap.dedent(body))


# ---------------------------------------------------------------------------
# meshes
# ---------------------------------------------------------------------------
def block(x0, x1, y0, y1, n1, n2, patches):
    verts = f"""
    ({x0} {y0} 0) ({x1} {y0} 0) ({x1} {y1} 0) ({x0} {y1} 0)
    ({x0} {y0} 0.1) ({x1} {y0} 0.1) ({x1} {y1} 0.1) ({x0} {y1} 0.1)
    """
    faces = {"bottom": "(0 1 5 4)", "top": "(3 7 6 2)", "left": "(0 4 7 3)", "right": "(1 2 6 5)"}
    bnd = ""
    for name, (ptype, side) in patches.items():
        bnd += f"    {name} {{ type {ptype}; faces ({faces[side]}); }}\n"
    bnd += "    frontAndBack { type empty; faces ((0 3 2 1) (4 5 6 7)); }\n"
    return f"""
scale 1;
vertices
(
{verts}
);
blocks ( hex (0 1 2 3 4 5 6 7) ({n1} {n2} 1) simpleGrading (1 1 1) );
edges ();
boundary
(
{bnd});
mergePatchPairs ();
"""


w("system/fluid/blockMeshDict", block(0, 2, 0, 1, nx, nyf, {
    "inlet": ("patch", "left"), "outlet": ("patch", "right"),
    "walls": ("wall", "bottom"), "interface": ("wall", "top")}))
w("system/solid/blockMeshDict", block(0, 2, 1, 1.5, nx, nys, {
    "interface": ("patch", "bottom"), "top": ("patch", "top"),
    "left": ("patch", "left"), "right": ("patch", "right")}))

if a.coupling == "solid":
    # Solid only: the velocity solid of the QM coupling with the exact interface
    # traction (mmsFsiSolidTraction)
    import shutil
    LIBS = {"full": "mmsFsi", "flat": "mmsFsiFlat", "mu10": "mmsFsiMu10"}[a.variant]
    os.rename(os.path.join(a.dir, "system/solid/blockMeshDict"), os.path.join(a.dir, "system/blockMeshDict"))
    shutil.rmtree(os.path.join(a.dir, "system/fluid")); shutil.rmtree(os.path.join(a.dir, "system/solid"))
    w("constant/physicsProperties", "type solid;\n")
    w("constant/g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n", "uniformDimensionedVectorField")
    w("constant/mechanicalProperties", """
planeStress     no;
mechanical ( layer { type neoHookeanElastic; rho rho [1 -3 0 0 0 0 0] 1; E E [1 -1 -2 0 0 0 0] 10; nu nu [0 0 0 0 0 0 0] 0.3; } );
""")
    w("constant/dynamicMeshDict", "dynamicFvMesh staticFvMesh;\n")
    w("constant/fvOptions", "mmsSource { type mmsFsiSource; region solid; fields (%s); }\n"
      % ("U" if a.solidModel == "velocity" else "D"))
    SM = {"velocity": "nonLinearGeometryTotalLagrangianVelocity",
          "totalDisp": "nonLinearGeometryTotalLagrangianTotalDisplacement"}[a.solidModel]
    w("constant/solidProperties", f"""
solidModel      {SM};
{SM}Coeffs
{{{SSG}
    solutionAlgorithm PETScSNES;
    stabilisation {{ momentum {{ type diffStencilLaplacian; scaleFactor {a.solidStab}; }} }}
    {a.solidExtra}
}}
""")
    wi = a.writeInterval if a.writeInterval else a.end
    w("system/controlDict", f"""
application     solids4Foam;
libs            ({LIBS});
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {a.end};
deltaT          {a.dt};
writeControl    adjustableRunTime;
writeInterval   {wi};
writeFormat     ascii;
writePrecision  12;
timeFormat      general;
timePrecision   10;
runTimeModifiable no;
functions
{{
    mmsErrors {{ type mmsFsiErrors; libs ({LIBS}); fluidRegion none; solidRegion region0; }}
}}
""")
    w("system/fvSchemes", """
d2dt2Schemes { default backward; }
ddtSchemes { default backward; }
gradSchemes { default leastSquaresS4f; }
divSchemes { default Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
snGradSchemes { default corrected; }
interpolationSchemes { default linear; }
""")
    w("system/fvSolution", """
solvers
{
    SOLVEDFIELD
    {
        solver    petsc;
        options
        {
            snes_type newtonls;
            snes_linesearch_type basic;
            snes_rtol "1e-10";
            snes_atol "1e-12";
            snes_stol "0";
            snes_max_it "60";
            snes_converged_reason;
            snes_mf_operator;
            ksp_type gmres;
            ksp_gmres_restart "200";
            ksp_rtol "1e-8";
            pc_type lu;
        }
    }
}
""".replace("SOLVEDFIELD", "U" if a.solidModel == "velocity" else "D")
       .replace("pc_type lu;", "pc_type lu;" if a.np == 1 else "pc_type redundant;\n            redundant_pc_type lu;"))
    def fieldS(rel, cls, dims, internal, bcs):
        body = f"dimensions {dims};\ninternalField {internal};\nboundaryField\n{{\n"
        for name, bc in bcs.items():
            body += f"    {name}\n    {{\n" + "".join(f"        {l}\n" for l in bc.strip().split("\n")) + "    }\n"
        body += "    frontAndBack { type empty; }\n}\n"
        w(rel, body, cls, os.path.basename(rel))
    clampS = "type fixedDisplacement;\nvalue uniform (0 0 0);"
    sideS = clampS if a.variant != "flat" else "type mmsFsiSolidDisplacement;\nvalue uniform (0 0 0);"
    fieldS("0/D", "volVectorField", "[0 1 0 0 0 0 0]", "uniform (0 0 0)", {
        "interface": "type mmsFsiSolidTraction;\ntraction uniform (0 0 0);\npressure uniform 0;\nvalue uniform (0 0 0);"
                     + ("\nnominalTraction true;\nuseUndeformedArea yes;" if a.nominalTraction else "")
                     + ("\n" + a.tractionExtra if a.tractionExtra else ""),
        "top": clampS, "left": sideS, "right": sideS})
    fieldS("0/pointD", "pointVectorField", "[0 1 0 0 0 0 0]", "uniform (0 0 0)", {
        p: "type calculated;\nvalue uniform (0 0 0);" for p in ("interface", "top", "left", "right")})
    fv0 = "type fixedValue;\nvalue uniform (0 0 0);"
    sideU = fv0 if a.variant != "flat" else "type mmsFsiSolidVelocity;\nvalue uniform (0 0 0);"
    if a.solidModel == "velocity":
        fieldS("0/U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", {
            "interface": "type zeroGradient;", "top": fv0, "left": sideU, "right": sideU})
    RUNS = "solids4Foam > log.solids4Foam 2>&1\n"
    if a.np > 1:
        w("system/decomposeParDict", f"""
numberOfSubdomains {a.np};
method          simple;
coeffs {{ n ({a.np} 1 1); }}
simpleCoeffs {{ n ({a.np} 1 1); delta 0.001; }}
""")
        RUNS = ("decomposePar > log.decomposePar 2>&1\n"
                f"mpirun -np {a.np} solids4Foam -parallel > log.solids4Foam 2>&1\n")
    with open(os.path.join(a.dir, "Allrun"), "w") as fh:
        fh.write("#!/bin/bash\ncd \"${0%/*}\" || exit\nsource ENVF > /dev/null\nblockMesh > log.blockMesh 2>&1\n"
                 .replace("ENVF", a.env) + RUNS)
    os.chmod(os.path.join(a.dir, "Allrun"), 0o755)
    print("written", a.dir)
    raise SystemExit(0)

if a.coupling == "fluid":
    # Fluid only: newtonIcoFluid on the fluid mesh, with the exact interface
    # motion and velocity prescribed (no solid)
    import shutil
    LIBF = {"full": "mmsFsi", "flat": "mmsFsiFlat", "mu10": "mmsFsiMu10"}[a.variant]
    GEN = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "derive",
                          {"full": "generated", "flat": "generated_flat", "mu10": "generated_mu10"}[a.variant]))
    os.rename(os.path.join(a.dir, "system/fluid/blockMeshDict"), os.path.join(a.dir, "system/blockMeshDict"))
    shutil.rmtree(os.path.join(a.dir, "system/fluid")); shutil.rmtree(os.path.join(a.dir, "system/solid"))
    w("constant/physicsProperties", "type fluid;\n")
    w("constant/g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n", "uniformDimensionedVectorField")
    w("constant/transportProperties", """
transportModel  Newtonian;
nu              nu [0 2 -1 0 0 0 0] %s;
rho             rho [1 -3 0 0 0 0 0] 1;
""" % ("0.1" if a.variant == "mu10" else "0.01"))
    w("constant/turbulenceProperties", "simulationType laminar;\n")
    w("constant/fvOptions", "mmsSource { type mmsFsiSource; region fluid; fields (U); }\n")
    w("constant/dynamicMeshDict", f"""
dynamicFvMesh   dynamicMotionSolverFvMesh;
motionSolver    {a.motion};
diffusivity     {a.diffusivity};
""")
    PIMPLE = a.fluidModel == "pimpleFluid"
    if PIMPLE:
        w("constant/fluidProperties", "fluidModel    pimpleFluid;\npimpleFluidCoeffs {}\n")
    else:
      w("constant/fluidProperties", f"""
fluidModel    newtonIcoFluid;
newtonIcoFluidCoeffs
{{
    forceImplicitFlux yes;{RC}
    fluidFluxExtrapolationAlgorithm1 yes;
    stabilisation {{ type RhieChow; scaleFactor {a.pStab}; }}
}}
""")
    wi = a.writeInterval if a.writeInterval else a.end
    w("system/controlDict", f"""
application     solids4Foam;
libs            ({LIBF});
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {a.end};
deltaT          {a.dt};
writeControl    adjustableRunTime;
writeInterval   {wi};
writeFormat     ascii;
writePrecision  12;
timeFormat      general;
timePrecision   10;
runTimeModifiable no;
functions
{{
    mmsErrors {{ type mmsFsiErrors; libs ({LIBF}); fluidRegion region0; solidRegion none; writeFields true; }}
}}
""")
    if PIMPLE:
        w("system/fvSchemes", PIMPLE_SCHEMES)
        w("system/fvSolution", PIMPLE_SOLUTION)
    else:
        w("system/fvSchemes", """
ddtSchemes { default none; ddt(U) backward; ddt(dU) backward; ddt(p) backward; }
gradSchemes { default GRADSCHEME; }
divSchemes
{
    default            none;
    jacobian-div(phi,U) Gauss upwind;
    div(phi,U)         Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
    div(U)             Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default INTERP; }
snGradSchemes { default corrected; }
fluxRequired { default no; p; }
""".replace("INTERP", a.interp).replace("GRADSCHEME", "leastSquares" if a.grad == "leastSquares"
             else "leastSquares;\n    grad(U) %s;\n    grad(p) %s" % (a.grad, a.grad)))
        w("system/fvSolution", f"""
solvers
{{
    Up
    {{
        solver    petsc;
        options
        {{
            snes_type newtonls;
            snes_linesearch_type basic;
            snes_rtol "{a.snesRtol}";
            snes_atol "1e-12";
            snes_stol "0";
            snes_max_it "60";
            snes_monitor;
            snes_converged_reason;
            snes_max_linear_solve_fail "100";
            snes_mf_operator;
            ksp_type lgmres;
            ksp_gmres_restart "200";
            ksp_rtol "{a.kspRtol}";
            ksp_max_it "1000";
            pc_type fieldsplit;
            pc_fieldsplit_type schur;
            pc_fieldsplit_schur_fact_type upper;
            pc_fieldsplit_schur_precondition selfp;
            pc_fieldsplit_block_size "3";
            pc_fieldsplit_0_fields "0,1";
            pc_fieldsplit_1_fields "2";
            fieldsplit_0_ksp_type preonly;
            fieldsplit_0_pc_type bjacobi;
            fieldsplit_0_sub_pc_type ilu;
            fieldsplit_1_ksp_type preonly;
            fieldsplit_1_pc_type hypre;
            fieldsplit_1_pc_hypre_type boomeramg;
        }}
    }}
    "cellDisplacement|dU|dp"
    {{
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-12;
        relTol          0;
    }}
}}
""")
    def fieldF(rel, cls, dims, internal, bcs):
        body = f"dimensions {dims};\ninternalField {internal};\nboundaryField\n{{\n"
        for name, bc in bcs.items():
            body += f"    {name}\n    {{\n" + "".join(f"        {l}\n" for l in bc.strip().split("\n")) + "    }\n"
        body += "    frontAndBack { type empty; }\n}\n"
        w(rel, body, cls, os.path.basename(rel))
    fieldF("0/U", "volVectorField", "[0 1 -1 0 0 0 0]", "uniform (0 0 0)", {
        "inlet": "type mmsFsiVelocity;\nvalue uniform (0 0 0);",
        "outlet": OUTLET_U,
        "walls": ("type fixedValueCorrected;\nvalue uniform (0 0 0);" if a.uCorrected
                  else "type fixedValue;\nvalue uniform (0 0 0);"),
        # pimpleFluid: the moving-wall velocity of the partitioned couplings, so that the
        # interface flux equals the mesh flux of the prescribed (exact) interface motion
        "interface": ("type newMovingWallVelocity;\nvalue uniform (0 0 0);" if PIMPLE
                      else "type mmsFsiVelocityCorrected;\nvalue uniform (0 0 0);" if a.uCorrected
                      else "type mmsFsiVelocity;\nvalue uniform (0 0 0);")})
    FFP = "type fixedFluxPressure;\nvalue uniform 0;"
    fieldF("0/p", "volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0", {
        # pimpleFluid: fixedFluxPressure on the inlet, walls and interface, as
        # in the partitioned couplings
        "inlet": FFP if PIMPLE else "type mmsFsiPressureGradient;\nvalue uniform 0;",
        "outlet": "type mmsFsiPressure;\nvalue uniform 0;",
        "walls": FFP if PIMPLE else
                 ("type fixedGradientCorrected;\ngradient uniform 0;\nvalue uniform 0;" if a.pCorrected
                  else "type fixedGradient;\ngradient uniform 0;"),
        "interface": (("type movingWallPressure;\nvalue uniform 0;" if a.interfaceP == "movingWallPressure"
                       else "type zeroGradient;" if a.interfaceP == "zeroGradient" else FFP) if PIMPLE else
                     ("type fixedGradientCorrected;\ngradient uniform 0;\nvalue uniform 0;" if a.pCorrected
                      else "type fixedGradient;\ngradient uniform 0;"))})
    ifc = """type codedFixedValue;
value uniform (0 0 0);
name exactInterfaceMotion;
codeInclude #{ #include "mmsFsiExact.H" #};
codeOptions #{ -I%s #};
code
#{
    const scalar t = this->db().time().value();
    const pointField& pts = this->patch().localPoints();
    vectorField disp(pts.size());
    forAll(pts, i)
    {
        const scalar X = pts[i].x() - (*this)[i].x();
        double x[2], tr[2], n[2];
        mmsFsi::interface(X, t, x, tr, n);
        disp[i] = vector(x[0] - X, x[1] - mmsFsi::param_H, 0);
    }
    this->operator==(disp);
#};""" % GEN
    fieldF("0/pointDisplacement", "pointVectorField", "[0 1 0 0 0 0 0]", "uniform (0 0 0)", {
        "inlet": endMotion(a.variant), "outlet": endMotion(a.variant),
        "walls": "type fixedValue;\nvalue uniform (0 0 0);", "interface": ifc})
    with open(os.path.join(a.dir, "Allrun"), "w") as fh:
        fh.write("#!/bin/bash\ncd \"${0%/*}\" || exit\nsource ENVF > /dev/null\nblockMesh > log.blockMesh 2>&1\n"
                 "solids4Foam > log.solids4Foam 2>&1\n".replace("ENVF", a.env))
    os.chmod(os.path.join(a.dir, "Allrun"), 0o755)
    print("written", a.dir)
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# constant
# ---------------------------------------------------------------------------
w("constant/physicsProperties", "type fluidSolidInteraction;\n")
w("constant/g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n", "uniformDimensionedVectorField")
w("constant/solid/g", "dimensions [0 1 -2 0 0 0 0];\nvalue (0 0 0);\n", "uniformDimensionedVectorField")
w("constant/fluid/transportProperties", """
transportModel  Newtonian;
nu              nu [0 2 -1 0 0 0 0] {NU};
rho             rho [1 -3 0 0 0 0 0] 1;
""".replace("{NU}", "0.1" if a.variant == "mu10" else "0.01"))
w("constant/fluid/turbulenceProperties", "simulationType laminar;\n")
w("constant/solid/mechanicalProperties", """
planeStress     no;
mechanical
(
    layer
    {
        type            neoHookeanElastic;
        rho             rho [1 -3 0 0 0 0 0] 1;
        E               E [1 -1 -2 0 0 0 0] 10;
        nu              nu [0 0 0 0 0 0 0] 0.3;
    }
);
""")
w("constant/solid/dynamicMeshDict", "dynamicFvMesh staticFvMesh;\n")
w("constant/fluid/fvOptions", """
mmsSource
{
    type        mmsFsiSource;
    region      fluid;
    fields      (U);
}
""")
w("constant/solid/fvOptions", f"""
mmsSource
{{
    type        mmsFsiSource;
    region      solid;
    fields      ({'U' if QM else 'D'});
}}
""")
w("constant/fluid/dynamicMeshDict", f"""
dynamicFvMesh   dynamicMotionSolverFvMesh;
motionSolver    {a.motion};
diffusivity     {a.diffusivity};
""")

if QM:
    w("constant/fsiProperties", """
fluidSolidInterface NewtonQuasiMonolithic;
NewtonQuasiMonolithicCoeffs
{
    solidPatch interface;
    fluidPatch interface;
    coupled         yes;
    interfaceTransferMethod directMap;
    passViscousStress yes;
    interfaceFluxDiagnostics yes;
    fluidSystemScaleFactor 1e+08;
}
""")
    w("constant/fluid/fluidProperties", f"""
fluidModel    newtonIcoFluid;
newtonIcoFluidCoeffs
{{
    forceImplicitFlux yes;{RC}
    fluidFluxExtrapolationAlgorithm1 yes;
    addDivPhiUDamping no;
    pressureScaleFactor 1.0;
    stabilisation
    {{
        type RhieChow;
        scaleFactor {a.pStab};
    }}
}}
""")
    w("constant/solid/solidProperties", f"""
solidModel      nonLinearGeometryTotalLagrangianVelocity;
nonLinearGeometryTotalLagrangianVelocityCoeffs
{{{SSG}
    stabilisation
    {{
        momentum
        {{
            type            diffStencilLaplacian;
            scaleFactor     {a.solidStab};
        }}
    }}
}}
""")
else:
    coupling = {"iqnils": "IQNILS", "aitken": "Aitken", "robin": "fixedRelaxation"}[a.coupling]
    extra = {
        "IQNILS": "relaxationFactor %s;\n    couplingReuse 0;\n    relMinSignificant 1e-2;\n    qrSolveTolerance 1e-3;\n    reorthogonalizeCouplingColumns yes;\n    predictSolid no;\n" % a.relax,
        "Aitken": "relaxationFactor %s;\n" % a.relax,
        "fixedRelaxation": "relaxationFactor 1;\n    predictSolid yes;\n    robinPressureTolerance 1e-6;\n    robinFluxTolerance 1e-5;\n",
    }[coupling]
    w("constant/fsiProperties", f"""
fluidSolidInterface {coupling};
{coupling}Coeffs
{{
    solidPatch interface;
    fluidPatch interface;
    {extra}
    outerCorrTolerance {a.outerTol};
    nOuterCorr 200;
    interfaceTransferMethod directMap;
    coupled yes;
    predictor yes;
    writeResidualsToFile yes;
    {a.fsiExtra}
}}
""")
    w("constant/fluid/fluidProperties", """
fluidModel    pimpleFluid;
pimpleFluidCoeffs {}
""")
    w("constant/solid/solidProperties", f"""
solidModel      nonLinearGeometryTotalLagrangianTotalDisplacement;
nonLinearGeometryTotalLagrangianTotalDisplacementCoeffs
{{
    solutionAlgorithm    PETScSNES;
    predictor            yes;
    stabilisation
    {{
        momentum
        {{
            type            diffStencilLaplacian;
            scaleFactor     {a.solidStab};
        }}
    }}
}}
""")

# ---------------------------------------------------------------------------
# system
# ---------------------------------------------------------------------------
wi = a.writeInterval if a.writeInterval else a.end
LIB = {"full": "mmsFsi", "flat": "mmsFsiFlat", "mu10": "mmsFsiMu10"}[a.variant]
w("system/controlDict", f"""
application     solids4Foam;
libs            ({LIB});
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {a.end};
deltaT          {a.dt};
writeControl    adjustableRunTime;
writeInterval   {wi};
purgeWrite      0;
writeFormat     ascii;
writePrecision  12;
writeCompression off;
timeFormat      general;
timePrecision   10;
runTimeModifiable no;
adjustTimeStep  no;
functions
{{
    mmsErrors
    {{
        type            mmsFsiErrors;
        libs            ({LIB});
        writeFields     true;
    }}
}}
""")

if QM:
    w("system/fvSolution", f"""
solvers
{{
    UpU
    {{
        solver    petsc;
        options
        {{
            snes_type newtonls;
            snes_linesearch_type basic;
            snes_rtol "{a.snesRtol}";
            snes_atol "1e-12";
            snes_stol "0";
            snes_max_it "60";
            snes_monitor;
            snes_converged_reason;
            snes_max_linear_solve_fail "100";
            snes_mf_operator;
            ksp_type lgmres;
            ksp_gmres_restart "200";
            ksp_rtol "{a.kspRtol}";
            ksp_max_it "1000";
            pc_type fieldsplit;
            pc_fieldsplit_type multiplicative;
            fieldsplit_solid_ksp_type preonly;
            fieldsplit_solid_pc_type redundant;
            fieldsplit_solid_redundant_pc_type lu;
            fieldsplit_fluid_ksp_type preonly;
            fieldsplit_fluid_pc_type fieldsplit;
            fieldsplit_fluid_pc_fieldsplit_type schur;
            fieldsplit_fluid_pc_fieldsplit_schur_fact_type upper;
            fieldsplit_fluid_pc_fieldsplit_schur_precondition selfp;
            fieldsplit_fluid_pc_fieldsplit_block_size "3";
            fieldsplit_fluid_pc_fieldsplit_0_fields "0,1";
            fieldsplit_fluid_pc_fieldsplit_1_fields "2";
            fieldsplit_fluid_fieldsplit_0_ksp_type preonly;
            fieldsplit_fluid_fieldsplit_0_pc_type bjacobi;
            fieldsplit_fluid_fieldsplit_0_sub_pc_type ilu;
            fieldsplit_fluid_fieldsplit_0_sub_pc_factor_levels "0";
            fieldsplit_fluid_fieldsplit_1_ksp_type preonly;
            fieldsplit_fluid_fieldsplit_1_pc_type hypre;
            fieldsplit_fluid_fieldsplit_1_pc_hypre_type boomeramg;
        }}
    }}
}}
""")
    w("system/fluid/fvSchemes", """
d2dt2Schemes { default none; }
ddtSchemes
{
    default            none;
    ddt(U)             backward;
    ddt(dU)            backward;
    ddt(p)             backward;
}
gradSchemes
{
    default            GRADSCHEME;
}
divSchemes
{
    default            none;
    jacobian-div(phi,U) Gauss upwind;
    div(phi,U)         Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
    div((nuEff*dev(T(grad(U))))) Gauss linear;
    div(U)             Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default INTERP; }
snGradSchemes { default corrected; }
fluxRequired { default no; p; }
""".replace("INTERP", a.interp).replace("GRADSCHEME", "leastSquares" if a.grad == "leastSquares"
             else "leastSquares;\n    grad(U) %s;\n    grad(p) %s" % (a.grad, a.grad)))
    w("system/solid/fvSchemes", """
d2dt2Schemes { default backward; }
ddtSchemes { default backward; }
gradSchemes { default leastSquaresS4f; }
divSchemes { default Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
snGradSchemes { default corrected; }
interpolationSchemes { default linear; }
""")
    w("system/solid/fvSolution", "solvers {}\n")
else:
    w("system/fluid/fvSchemes", PIMPLE_SCHEMES)
    w("system/solid/fvSchemes", """
d2dt2Schemes { default backward; }
ddtSchemes { default backward; }
gradSchemes { default leastSquaresS4f; }
divSchemes { default Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
snGradSchemes { default corrected; }
interpolationSchemes { default linear; }
""")
    w("system/solid/fvSolution", """
solvers
{
    D
    {
        solver    petsc;
        options
        {
            snes_type newtonls;
            snes_linesearch_type basic;
            snes_rtol "1e-12";
            snes_atol "1e-14";
            snes_stol "0";
            snes_max_it "50";
            snes_mf_operator;
            ksp_type lgmres;
            ksp_gmres_restart "100";
            ksp_rtol "1e-8";
            pc_type lu;
        }
    }
}
""")

w("system/fluid/fvSolution", PIMPLE_SOLUTION)


# ---------------------------------------------------------------------------
# 0/
# ---------------------------------------------------------------------------
def field(rel, cls, dims, internal, bcs):
    body = f"dimensions {dims};\ninternalField {internal};\nboundaryField\n{{\n"
    for name, bc in bcs.items():
        body += f"    {name}\n    {{\n" + "".join(f"        {l}\n" for l in bc.strip().split("\n")) + "    }\n"
    body += "    frontAndBack { type empty; }\n}\n"
    w(rel, body, cls, os.path.basename(rel))


V0 = "uniform (0 0 0)"
if QM:
    Uint, pint = ("type fixedValueCorrected;\nvalue uniform (0 0 0);" if a.uCorrected
                  else "type fixedValue;\nvalue uniform (0 0 0);"), \
        ("type mmsFsiPressureOracle;\nmode %s;\ngradient uniform 0;\nvalue uniform 0;" % a.pOracle if a.pOracle
         else "type fixedGradientCorrected;\ngradient uniform 0;\nvalue uniform 0;" if a.pCorrected
         else "type fixedGradient;\ngradient uniform 0;")
elif a.coupling == "robin":
    Uint, pint = "type elasticWallVelocity;\nvalue uniform (0 0 0);", \
        "type elasticWallPressure;\nvalue uniform 0;"
else:
    Uint, pint = "type newMovingWallVelocity;\nvalue uniform (0 0 0);", \
        ("type zeroGradient;" if a.interfaceP == "zeroGradient"
         else "type movingWallPressure;\nvalue uniform 0;" if a.interfaceP == "movingWallPressure"
         else "type fixedFluxPressure;\nvalue uniform 0;")
if QM:
    pwall = (("type fixedGradient;\ngradient uniform 0;" if a.pOracle else pint)
             if a.wallP == "fixedGradient" else "type zeroGradient;")
else:
    pwall = "type zeroGradient;" if a.wallP == "zeroGradient" else "type fixedFluxPressure;\nvalue uniform 0;"
field("0/fluid/U", "volVectorField", "[0 1 -1 0 0 0 0]", V0, {
    "inlet": "type mmsFsiVelocity;\nvalue uniform (0 0 0);",
    "outlet": OUTLET_U,
    "walls": ("type fixedValueCorrected;\nvalue uniform (0 0 0);" if (QM and a.uCorrected)
              else "type fixedValue;\nvalue uniform (0 0 0);"),
    "interface": Uint})
field("0/fluid/p", "volScalarField", "[0 2 -2 0 0 0 0]", "uniform 0", {
    # QM: exact Neumann p; pimpleFluid: fixedFluxPressure (an exact Neumann p
    # with fixedValue U makes the pressure-corrected inlet flux inconsistent)
    "inlet": ("type mmsFsiPressureGradient;\nvalue uniform 0;" if QM
              else "type fixedFluxPressure;\nvalue uniform 0;"),
    "outlet": "type mmsFsiPressure;\nvalue uniform 0;",
    "walls": pwall,
    "interface": pint})
ends = endMotion(a.variant)
if a.motion.startswith("displacement"):
    field("0/fluid/pointDisplacement", "pointVectorField", "[0 1 0 0 0 0 0]", V0, {
        "inlet": ends, "outlet": ends, "walls": "type fixedValue;\nvalue uniform (0 0 0);",
        "interface": "type fixedValue;\nvalue uniform (0 0 0);"})
else:
    field("0/fluid/pointMotionU", "pointVectorField", "[0 1 -1 0 0 0 0]", V0, {
        p: "type fixedValue;\nvalue uniform (0 0 0);" for p in ("inlet", "outlet", "walls", "interface")})

clamp = "type fixedDisplacement;\nvalue uniform (0 0 0);"
side = clamp if a.variant != "flat" else "type mmsFsiSolidDisplacement;\nvalue uniform (0 0 0);"
field("0/solid/D", "volVectorField", "[0 1 0 0 0 0 0]", V0, {
    "interface": "type solidTraction;\ntraction uniform (0 0 0);\npressure uniform 0;\nvalue uniform (0 0 0);",
    "top": clamp, "left": side, "right": side})
field("0/solid/pointD", "pointVectorField", "[0 1 0 0 0 0 0]", V0, {
    p: "type calculated;\nvalue uniform (0 0 0);" for p in ("interface", "top", "left", "right")})
if QM:
    fv0 = "type fixedValue;\nvalue uniform (0 0 0);"
    sideU = fv0 if a.variant != "flat" else "type mmsFsiSolidVelocity;\nvalue uniform (0 0 0);"
    field("0/solid/U", "volVectorField", "[0 1 -1 0 0 0 0]", V0, {
        "interface": "type zeroGradient;", "top": fv0, "left": sideU, "right": sideU})

if a.np > 1:
    for reg in ("", "fluid/", "solid/"):
        w(f"system/{reg}decomposeParDict", f"""
numberOfSubdomains {a.np};
method          simple;
coeffs {{ n ({a.np} 1 1); }}
simpleCoeffs {{ n ({a.np} 1 1); delta 0.001; }}
""")
RUN = ("solids4Foam > log.solids4Foam 2>&1\n" if a.np == 1 else
       "decomposePar -region fluid > log.decomposePar.fluid 2>&1\n"
       "decomposePar -region solid > log.decomposePar.solid 2>&1\n"
       f"mpirun -np {a.np} solids4Foam -parallel > log.solids4Foam 2>&1\n")
with open(os.path.join(a.dir, "Allrun"), "w") as fh:
    fh.write("""#!/bin/bash
cd "${0%/*}" || exit
# DYLD_* variables do not survive a new /bin/bash on macOS: source here
source ENVFILE > /dev/null
blockMesh -region fluid > log.blockMesh.fluid 2>&1
blockMesh -region solid > log.blockMesh.solid 2>&1
RUNCMD""".replace("ENVFILE", a.env).replace("RUNCMD", RUN))
os.chmod(os.path.join(a.dir, "Allrun"), 0o755)
print("written", a.dir)
