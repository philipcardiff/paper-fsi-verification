#!/usr/bin/env python3
"""FSI1 CFD1 refinement gate (preregistration B/B2): table from stageC_analysis.json files."""
import json, os, math
R = "<project-root>/"; O = R + "OVERNIGHT2/runs/"
runs = [("L1 (1.25) pinned [stage C]", R + "runs/hronturekBC/C1p_L1_f1.25_K1e6"), ("L1.5 (1.875) pinned (rule)", O + "C2p_f1.875"), ("L2 (2.5) pinned (rule)", O + "C2p_f2.5"),
        ("L1.5 pins = L1 values", O + "C2pin_f1.875_P1"), ("L1.5 pins = 0", O + "C2pin_f1.875_P0"), ("L1.5 pins = 2x rule", O + "C2pin_f1.875_P2"),
        ("L1.5 unpinned", O + "C2u_f1.875"), ("L2 unpinned", O + "C2u_f2.5")]
rows = []
for lab, d in runs:
    p = d + "/stageC_analysis.json"
    if not os.path.exists(p): continue
    j = json.load(open(p)); J = j.get("junction", {})
    rows.append(dict(run=lab, cells=j["n_cells"], F1_drag=j["F1"][0], F1_lift=j["F1"][1], dragerr_F1=100 * j["drag_err_F1"], lifterr_F1=100 * j["lift_err_F1"],
                     dragerr_F2=100 * j["drag_err_F2"], lifterr_F2=100 * j["lift_err_F2"], dragerr_F3=100 * j["drag_err_F3"], lifterr_F3=100 * j["lift_err_F3"],
                     at_top=j["parity_top"]["a_t"], at_bot=j["parity_bottom"]["a_t"], at_cyl=j["parity_cyl_at"], mass=j["mass_imbalance"], spread_drag=j["spread_drag"],
                     pins=J.get("p_junction")))
json.dump(rows, open(R + "OVERNIGHT2/b_gate_table.json", "w"), indent=1)
md = "| run | cells | F1 drag err % | F1 lift err % | F2 drag/lift err % | F3 drag/lift err % | a_t top/bottom/cyl (Pa) | mass imb | drag spread | junction p (Pa) |\n|---|---|---|---|---|---|---|---|---|---|\n"
for x in rows:
    md += "| %s | %d | %+.2f | %+.1f | %+.2f / %+.2f | %+.2f / %+.1f | %.2f / %.2f / %.2f | %.1e | %.1e | %s |\n" % (x["run"], x["cells"], x["dragerr_F1"], x["lifterr_F1"], x["dragerr_F2"], x["lifterr_F2"], x["dragerr_F3"], x["lifterr_F3"], x["at_top"], x["at_bot"], x["at_cyl"], x["mass"], x["spread_drag"], x["pins"])
open(R + "OVERNIGHT2/b_gate_table.md", "w").write(md); print(md)
