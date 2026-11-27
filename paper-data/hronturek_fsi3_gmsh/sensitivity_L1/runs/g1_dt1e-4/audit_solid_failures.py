"""Audit solid PETSc SNES/KSP failures in a solids4Foam FSI log.
For each failed solid solve (the SNES after "Solving the momentum equation for D using PETSc SNES"),
record time step and coupling iteration, then check whether a later coupling iteration of the same
time step had a converged solid solve and whether the step met the coupling tolerance (combined
residual <= tol) or was accepted unconverged.
Usage: audit_solid_failures.py TOL LOG..."""
import re, sys
tol = float(sys.argv[1])
fails = []          # (time, iter, reason)
steps = {}          # time -> dict(last_iter, last_res, solid_ok_after_fail)
t = it = None; block = None; cur_ok = True
re_t = re.compile(r"^Time = ([0-9.eE+-]+), iteration: (\d+)")
for log in sys.argv[2:]:
    with open(log, errors="replace") as fh:
        for line in fh:
            m = re_t.match(line)
            if m:
                t, it = m.group(1), int(m.group(2)); block = None
                s = steps.setdefault(t, {"iters": [], "res": None, "fail_iters": [], "ok_after_fail": []})
                s["iters"].append(it); continue
            if line.startswith("Solving the momentum equation for D using PETSc SNES"):
                block = "solid"; cur_ok = True; continue
            if line.startswith("Evolving the mesh motion") or line.startswith("Evolving fluid model"):
                block = "other"; continue
            if block == "solid":
                if "did not converge" in line or "DIVERGED" in line:
                    if cur_ok:
                        fails.append((t, it, line.strip())); steps[t]["fail_iters"].append(it)
                    cur_ok = False
                elif "Nonlinear solve converged" in line:
                    if t in steps and steps[t]["fail_iters"]: steps[t]["ok_after_fail"].append(it)
                    block = None
            if t is not None and line.startswith("FSI combined residual norm"):
                steps[t]["res"] = float(line.split(":")[-1])
print(f"solid SNES/KSP failures: {len(fails)}")
bad = []
for (ft, fi, why) in fails:
    s = steps[ft]
    later_ok = [i for i in s["ok_after_fail"] if i > fi]
    conv = s["res"] is not None and s["res"] <= tol
    flag = "OK" if later_ok and conv else "ACCEPTED_WITH_FAILED_SOLID"
    if flag != "OK": bad.append(ft)
    print(f"t={ft} iter={fi}: {why} | later converged solid iters {later_ok[:5]} | final combined residual {s['res']} -> {flag}")
print(f"steps accepted with a failed solid solve (no later converged solid iteration and converged coupling): {len(set(bad))} {sorted(set(bad), key=float)[:20]}")
