#!/usr/bin/env python3
"""Newton-iteration and termination statistics from the FEBio text logs (<tag>.log, plain or .xz).
usage: log_stats.py RAW_DIR OUT_DIR
Writes OUT_DIR/run_summary.csv (one row per run and per preserved failed attempt) and OUT_DIR/newton_iterations_per_step.csv (run, step, time, iterations, reformations, converged)
and OUT_DIR/rtc_m1_step2_residuals.csv (squared residual norm of every Newton iteration of the failing step of RTC_fsi_M1)."""
import sys, os, io, re, csv, lzma, glob


def read(path):
    if path.endswith(".xz"):
        return io.TextIOWrapper(lzma.open(path), errors="replace").read()
    return open(path, errors="replace").read()


def parse(text):
    steps = []          # dict(step, time, iters, refs, converged, residuals[list of (initial, current)])
    cur = None
    for ln in text.splitlines():
        m = re.match(r"===== beginning time step (\d+) : ([0-9.eE+-]+) =====", ln)
        if m:
            cur = dict(step=int(m.group(1)), time=float(m.group(2)), iters=None, refs=None, converged=None, residuals=[]); steps.append(cur); continue
        if cur is None:
            continue
        m = re.match(r"\s+residual\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)", ln)
        if m:
            cur["residuals"].append((float(m.group(1)), float(m.group(2))))
        m = re.match(r"\s+number of iterations\s+:\s+(\d+)", ln)
        if m:
            cur["iters"] = int(m.group(1))
        m = re.match(r"\s+number of reformations\s+:\s+(\d+)", ln)
        if m:
            cur["refs"] = int(m.group(1))
        if "------- converged at time" in ln:
            cur["converged"] = True
        if "failed to converge at time" in ln:
            cur["converged"] = False
    term = "NORMAL" if re.search(r"N O R M A L\s+T E R M I N A T I O N", text) else ("ERROR" if re.search(r"E R R O R\s+T E R M I N A T I O N", text) else "NONE")
    msgs = sorted(set(re.findall(r"\* ([^*\n]*(?:negative jacobian|Max\. nr\. of reformations|max\. nr\. of reformations|Problem is diverging)[^*\n]*?)\s*\*", text, flags=re.I)))
    m = re.search(r"Total number of equilibrium iterations \.+ : (\d+)", text)
    tot_it = int(m.group(1)) if m else None
    return steps, term, msgs, tot_it


def main():
    raw, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    rows, per = [], []
    for run in sorted(os.listdir(raw)):
        logs = [(run, "final", os.path.join(raw, run, run + ".log.xz"))]
        logs += [(run, os.path.basename(p)[:-7], p) for p in sorted(glob.glob(os.path.join(raw, run, "failed", "*.log.xz")))]
        for r, kind, p in logs:
            if not os.path.exists(p):
                continue
            steps, term, msgs, tot = parse(read(p))
            conv = [s for s in steps if s["converged"]]
            its = [s["iters"] for s in conv if s["iters"] is not None]
            rows.append(dict(run=r, attempt=kind, steps_started=len(steps), steps_converged=len(conv), last_time_started=steps[-1]["time"] if steps else "",
                             termination=term, mean_iterations_per_converged_step=round(sum(its) / len(its), 3) if its else "", max_iterations=max(its) if its else "",
                             total_equilibrium_iterations=tot if tot is not None else "", diagnostics="; ".join(msgs)))
            if kind == "final":
                for s in steps:
                    per.append(dict(run=r, step=s["step"], time=s["time"], iterations=s["iters"] if s["iters"] is not None else "", reformations=s["refs"] if s["refs"] is not None else "",
                                    converged=s["converged"]))
            if r == "RTC_fsi_M1" and kind == "final":
                bad = [s for s in steps if s["converged"] is False][-1]
                with open(os.path.join(out, "rtc_m1_step2_residuals.csv"), "w", newline="") as f:
                    w = csv.writer(f); w.writerow(["iteration_status_line", "squared_residual_initial", "squared_residual_current"])
                    for i, (a, b) in enumerate(bad["residuals"], 1):
                        w.writerow([i, a, b])
    with open(os.path.join(out, "run_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(os.path.join(out, "newton_iterations_per_step.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(per[0])); w.writeheader(); w.writerows(per)
    for r in rows:
        print(r["run"], r["attempt"], r["termination"], r["steps_converged"], r["mean_iterations_per_converged_step"], r["max_iterations"], r["diagnostics"])


if __name__ == "__main__":
    main()
