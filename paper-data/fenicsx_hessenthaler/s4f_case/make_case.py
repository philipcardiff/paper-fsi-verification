"""Create the isolated 'ours-s4f' feF1 case (Addendum A6/A6b deltas only) from the pinned tutorial extraction. usage: make_case.py OUTDIR"""
import sys, re, shutil, pathlib
src = pathlib.Path(__import__("os").environ["S4F_ROOT"]) / "tut/tutorials/fluidSolidInteraction/hessenthalerFsi"; out = pathlib.Path(sys.argv[1])
if out.exists(): sys.exit("exists")
shutil.copytree(src, out)
def sub(path, pat, rep, count=1, flags=re.M):
    p = out / path; t = p.read_text(); t2, n = re.subn(pat, rep, t, count=count, flags=flags)
    assert n >= 1, (path, pat); p.write_text(t2)
# (1) outlet velocity
p = out / "0/fluid/U"; t = p.read_text()
i = t.index("    outlet"); j = t.index("}", t.index("{", i)) + 1
t = t[:i] + "    outlet\n    {\n        type            backflowStabilisedOutletVelocity;\n        beta            0.2;\n        smoothing       0.01;\n        value           uniform (0 0 0);\n    }" + t[j:]
p.write_text(t)
# (2) schemes
sub("system/fluid/fvSchemes", r"^(\s*default\s+)Gauss linear;", r"\1leastSquares;")
sub("system/fluid/fvSchemes", r"^\s*grad\(U\)\s+cellLimited Gauss linear 1;", "    grad(U)            leastSquares;")
sub("system/fluid/fvSchemes", r"div\(phi,U\)\s+Gauss linearUpwind grad\(U\);", "div(phi,U)         Gauss LUST grad(U);")
# (3) PIMPLE
sub("system/fluid/fvSolution", r"nOuterCorrectors\s+1;", "nOuterCorrectors    20;\n    outerCorrectorResidualControl\n    {\n        U { tolerance 1e-5; relTol 0; }\n        p { tolerance 1e-5; relTol 0; }\n    }")
# (4) FSI coupling
sub("constant/fsiProperties", r"nOuterCorr\s+30;", "nOuterCorr 60;")
sub("constant/fsiProperties", r"couplingReuse\s+5;", "couplingReuse 2;")
sub("constant/fsiProperties", r"(outerCorrTolerance\s+1e-4;)", r"\1\n    requireAllResidualMeasures yes;\n    allowUnconvergedCoupling yes;")
# (5) material
sub("constant/solid/mechanicalProperties", r"(\bmu\s+mu \[1 -1 -2 0 0 0 0\]\s+)58450;", r"\g<1>64200;")
sub("constant/solid/mechanicalProperties", r"(\bK\s+K \[1 -1 -2 0 0 0 0\]\s+)565017;", r"\g<1>620600;")
print("case ready", out)
