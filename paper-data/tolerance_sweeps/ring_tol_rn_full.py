import sys, json, csv
from pathlib import Path
sys.path.insert(0, '/scratch/pcardiff/mwp721/ring/tut/ringAddedMass/verification/scripts')
import ring_added_mass_verification as V
OMEGA_DRY = 2.7125759973  # dry ring, standard solid, mesh 2, 100 steps (Sonic run)
RATIO_EX = 0.768654490979269 / 2.55504219581415
w = csv.writer(sys.stdout)
w.writerow(["run","coupling","tol","omega","omega_dry","ratio_error_pct","loss_per_period","mean_iterations","max_iterations"])
for c in sys.argv[1:]:
    case = Path(c); spec = json.loads((case/'spec.json').read_text())
    t, v = V.read_histories(case, spec); r = V.frequency_and_damping(t, v, spec['periods'])
    it, mx = V.coupling_iterations(case)
    w.writerow([case.name, "rn", spec['extra']['tol'], f"{r['omega']:.10f}", OMEGA_DRY,
                f"{100*(r['omega']/OMEGA_DRY/RATIO_EX-1):.6f}", f"{r['loss_per_period']:.7f}", f"{it:.3f}", mx])
