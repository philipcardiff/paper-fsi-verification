#!/usr/bin/env python3
"""Regenerate the FEBio inputs with make_rtr_feb.py from the stored model_meta.json files and compare with inputs/INPUT_HASHES.csv.
usage: regen_inputs.py PACKAGE_ROOT OUT_DIR      (writes OUT_DIR/<tag>/{<tag>.feb, mesh_rtr.npz, model_meta.json}; exit 1 on an unexpected mismatch)"""
import csv, hashlib, json, os, subprocess, sys
root, out = sys.argv[1], sys.argv[2]
gen = os.path.join(root, "scripts", "make_rtr_feb.py")
flags = (("dt", "--dt"), ("t_end", "--t-end"), ("amp", "--amp"), ("k", "--k"), ("e_mesh", "--e-mesh"), ("rhoi", "--rhoi"), ("max_refs", "--max-refs"), ("three_field", "--three-field"),
         ("clamp_outer_mesh", "--clamp-outer-mesh"), ("no_solid", "--no-solid"), ("mu_f", "--mu-f"), ("solid", "--solid"), ("bc_time_shift", "--bc-time-shift"), ("min_res", "--min-res"))
bad = 0
for row in csv.DictReader(open(os.path.join(root, "inputs", "INPUT_HASHES.csv"))):
    t = row["tag"]
    mp = os.path.join(root, "results", "raw", t, "model_meta.json")
    if row["regenerable_bit_identical"].startswith("no") or not os.path.exists(mp):
        print(f"  {t}: not regenerable by the final generator ({row['regenerable_bit_identical']}); shipped as run" if row["regenerable_bit_identical"].startswith("no") else f"  {t}: no meta"); continue
    m = json.load(open(mp))
    a = [sys.executable, "-I", gen, "--mode", m["mode"], "--out", os.path.join(out, t), "--tag", t]
    a += ["--mesh", m["mesh"]] if m.get("mesh") else ["--nth", str(m["nth"]), "--nrs", str(m["nrs"]), "--nrf", str(m["nrf"])]
    for k, f in flags:
        if k in m and m[k] is not None:
            a += [f, str(m[k])]
    subprocess.run(a, check=True, capture_output=True)
    b = open(os.path.join(out, t, t + ".feb"), "rb").read()
    h = hashlib.sha256(b).hexdigest()
    if h == row["feb_sha256"]:
        print(f"  {t}: IDENTICAL (sha256 {h[:16]})")
    elif row["regenerable_bit_identical"].startswith("modulo"):
        shipped = os.path.join(root, "inputs", "final", t + ".feb")
        same = sorted(open(shipped, "rb").read().split(b"\n")) == sorted(b.split(b"\n"))
        print(f"  {t}: hash differs, identical modulo line order (NodeSet order): {same}")
        bad += 0 if same else 1
    else:
        print(f"  {t}: MISMATCH {h[:16]} vs {row['feb_sha256'][:16]}"); bad += 1
sys.exit(1 if bad else 0)
