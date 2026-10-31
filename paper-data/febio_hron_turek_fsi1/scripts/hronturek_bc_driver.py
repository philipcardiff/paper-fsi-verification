#!/usr/bin/env python3
"""Sequential, budget-guarded submission of the Hron-Turek FSI1 Stage B/C jobs (one FEBio job at a time).
usage: hronturek_bc_driver.py CORES EXPECTED_CPUH TIMELIMIT RUN_DIR [RUN_DIR ...]      (RUN_DIR = runs/hronturekBC/<case>, contains <case>.feb and model_meta.json)
For every run dir: pre-registers the ledger row (status 'submitting') if absent, checks squeue/ps (no other FEBio job of this workflow running), checks the budget
(cumulative sacct CPU-h + EXPECTED_CPUH must stay <= BUDGET=3.0), submits ONLY through tools/submit_guarded.sh, waits, appends the sacct-based budget row,
and updates the ledger from the LOG (never from the exit code). Stops at the first failure to submit; the first job with a non-normal termination is logged and the driver
continues only if CONTINUE_ON_FAIL=1 (default: stop)."""
import sys, os, json, csv, subprocess, time, hashlib, re
ROOT = os.environ.get("HRON_PROJECT_ROOT", os.getcwd())   # site-specific driver (Slurm + ledger helpers of the original workspace); NOT needed for the reproduction checks; BUDGET = 3.0
BIN_SHA = "f5ae1ed6f347f875c392f93d847f2369f437345a9b785296d46b0bd6fbb72f02"
BT = os.path.join(ROOT, "runs/hronturekBC/budget.csv")
def sh(c): return subprocess.run(c, capture_output=True, text=True)
def cum():
    if not os.path.exists(BT): return 0.0
    return sum(float(r["alloc_cpuh"]) for r in csv.DictReader(open(BT)))
def main():
    cores, exp, tl = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3]; dirs = sys.argv[4:]
    if not os.path.exists(BT): open(BT, "w").write("jobid,case,cores,elapsed_s,alloc_cpuh,state,cumulative_cpuh,termination\n")
    for d in dirs:
        d = os.path.realpath(d); case = os.path.basename(d); meta = json.load(open(os.path.join(d, "model_meta.json")))
        feb = case + ".feb"; h = hashlib.sha256(open(os.path.join(d, feb), "rb").read()).hexdigest()
        rel = os.path.relpath(d, ROOT)
        q = sh(["squeue", "-u", os.environ["USER"], "-h"]).stdout.strip(); ps = sh(["bash", "-c", "ps -eo pid,user,args | grep '[f]ebio4 -i'"]).stdout.strip()
        mine = [l for l in ps.splitlines() if os.environ["USER"] in l]
        if q or mine: print("STOP: queue/process not empty:", q, mine); sys.exit(2)
        c0 = cum()
        print(f"=== {case}: purpose={meta.get('purpose', 'see pre-registration')} cores={cores} expected={exp} CPU-h cumulative={c0:.4f} of {BUDGET}; sha256 {h[:12]}", flush=True)
        if c0 + exp > BUDGET: print("STOP: projected cumulative cost would exceed the budget"); sys.exit(3)
        rows = [r for r in csv.reader(open(os.path.join(ROOT, "MULTIBENCHMARK_JOB_LEDGER.csv"))) if len(r) > 11 and r[11] == rel]
        if not rows:
            subprocess.run([os.path.join(ROOT, "tools/ledger_add.sh"), "HronTurekFSI1", case, json.dumps({k: meta[k] for k in meta if k in ("E", "g", "root", "elem", "nt", "k", "level")}),
                            json.dumps({k: meta[k] for k in meta if k in ("n_elem", "n_nodes", "steps", "dt", "t_end")}), h, BIN_SHA, meta.get("mesh_sha256", h),
                            meta.get("purpose", "Hron-Turek FSI1 stage B/C (pre-registered in HRON_TUREK_FSI1_STAGE_BC_RESULTS.md)"), "submitting", "", str(cores), rel, "", "pre-registered"], check=True)
        r = sh([os.path.join(ROOT, "tools/submit_guarded.sh"), d, feb, str(cores), "main", tl, "hbc_" + case[:20], "12G"]); print(r.stdout.strip(), r.stderr.strip(), flush=True)
        m = re.search(r"SUBMITTED (\d+)", r.stdout)
        if not m: print("STOP: submission refused/failed"); sys.exit(4)
        jid = m.group(1)
        subprocess.run([os.path.join(ROOT, "tools/ledger_update.py"), rel, f"job_id_or_pid={jid}", "status=running"])
        while True:
            time.sleep(10)
            s = sh(["sacct", "-j", jid, "-n", "-P", "-X", "-o", "JobID,State,ElapsedRaw,AllocCPUS"]).stdout.strip()
            if s:
                f = s.splitlines()[0].split("|")
                if f[1].split()[0] not in ("PENDING", "RUNNING", "COMPLETING", "REQUEUED"): break
        jid_, state, el, ac = f[0], f[1], int(f[2]), int(f[3]); ch = el * ac / 3600
        lf = os.path.join(d, case + ".log"); log = open(lf, errors="replace").read() if os.path.exists(lf) else ""
        term = "NORMAL" if "N O R M A L   T E R M I N A T I O N" in log else ("ERROR" if "E R R O R   T E R M I N A T I O N" in log else "UNKNOWN")
        c1 = cum() + ch
        open(BT, "a").write(f"{jid},{case},{ac},{el},{ch:.5f},{state},{c1:.5f},{term}\n")
        st = "completed-normal" if term == "NORMAL" else "failed-error-termination" if term == "ERROR" else "unknown-termination"
        subprocess.run([os.path.join(ROOT, "tools/ledger_update.py"), rel, f"status={st}", f"validation_status=termination {term} (log inspected); analysis pending"])
        print(f"--- {case}: job {jid} {state} elapsed {el}s alloc {ch:.4f} CPU-h cumulative {c1:.4f}; termination {term}", flush=True)
        if term != "NORMAL" and os.environ.get("CONTINUE_ON_FAIL") != "1": print("STOP: non-normal termination; inspect the log"); sys.exit(5)
if __name__ == "__main__":
    main()
