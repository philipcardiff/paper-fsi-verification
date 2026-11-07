#!/usr/bin/env python3
"""Regenerate every FEBio input from the generators in scripts/ and compare the sha256 with the value recorded at run time (inputs/INPUTS.csv, results/raw/*/model_meta.json),
and the shipped input files with their recorded file hash.  Package-only (python3 + numpy).  Usage from the package root:  python3 -I scripts/check_inputs.py
Expected: 'identical' for every row except the two early failed inlet attempts (B0_fluid_f1, B0_fluid_straight), whose generator output differs from the shipped file only in the
diagnostic monitor node sets (see provenance/INCONSISTENCIES.md I-9)."""
import csv, json, os, sys, subprocess, tempfile, hashlib, gzip, shutil
def h(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
rows = list(csv.DictReader(open("inputs/INPUTS.csv"))); tmp = tempfile.mkdtemp(); bad = 0
def regenerate(r):
    gen = os.path.join("scripts", r["generator"]); args = json.loads(r["args"]); tag = r["tag"]
    if tag.endswith("_generator_output"): tag = {"T3_straight_generator_output": "T3_straight", "T4_fsi_generator_output": "T4_fsi"}[tag]
    out = os.path.join(tmp, tag); os.makedirs(out, exist_ok=True)
    if r["generator"] == "make_cavity_plate_feb.py": cmd = [sys.executable, "-I", gen] + args + ["--out", out, "--tag", tag]
    elif r["generator"] == "make_cavity_plate_hex8_multi.py":
        if "@C3_fsi_M2_STAGED" in args:   # stage the M2 interface log + input exactly as the original run directory
            st = os.path.join(tmp, "M2stage"); os.makedirs(st, exist_ok=True)
            shutil.copy("results/raw/coupled/C3_fsi_M2/model_meta.json", st); shutil.copy("inputs/C3_fsi_M2.feb", os.path.join(st, "C3_fsi_M2.feb"))
            with gzip.open("results/raw/coupled/C3_fsi_M2/Mon_iface.txt.gz", "rb") as fi, open(os.path.join(st, "Mon_iface.txt"), "wb") as fo: shutil.copyfileobj(fi, fo)
            args = [st if a == "@C3_fsi_M2_STAGED" else a for a in args]
        cmd = [sys.executable, "-I", gen] + args + ["--out", out, "--tag", tag]
    else: cmd = [sys.executable, "-I", gen, "--out", out, "--tag", tag] + args
    subprocess.run(cmd, check=True, capture_output=True); return os.path.join(out, tag + ".feb")
for r in rows:
    f = regenerate(r); a = h(f); ok = a == r["sha256"]
    sf = r["shipped_file"]; sh = ("shipped file hash " + ("OK" if h(sf) == r["file_sha256"] else "MISMATCH")) if sf else "not shipped (>3 MB or regenerated only)"
    print(f"{'identical' if ok else 'DIFFERENT'}  {r['tag']:42s} {sh}"); bad += (not ok)
print("regenerated inputs differing from the recorded hash:", bad, "(expected 2: early failed attempts B0_fluid_f1, B0_fluid_straight)")
