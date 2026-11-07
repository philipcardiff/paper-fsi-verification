"""Mesh/time convergence post-processing: observed order, Richardson extrapolation, GCI (Roache; Celik et al. 2008 ASME V&V 20 procedure).
Rules enforced (a result is only 'accepted' if all hold):
  - same formulation on all levels (caller's responsibility; recorded in the report),
  - >= 3 levels, monotone convergence (successive differences have the same sign and decrease in magnitude),
  - observed order p in [0.5, 3] (outside => 'not in asymptotic range', no extrapolation reported as accepted),
  - GCI uses Fs = 1.25 for >=3 levels.
Levels are given coarse -> fine with h (representative size) or cell counts (h ~ N^{-1/dim})."""
import math


def h_from_cells(n_cells, dim=2): return float(n_cells) ** (-1.0 / dim)


def observed_order(f, h, tol=1e-12, maxit=200):
    """f, h: three values coarse->fine (f3 coarsest... NOTE: here ordered coarse, medium, fine). Returns (p, status)."""
    if len(f) != 3 or len(h) != 3: raise ValueError("need exactly 3 levels")
    f3, f2, f1 = f; h3, h2, h1 = h          # 3 coarse, 2 medium, 1 fine
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 == 0: return float("nan"), "zero difference"
    r21, r32 = h2 / h1, h3 / h2
    s = 1.0 if e32 / e21 > 0 else -1.0
    if s < 0: oscillatory = True
    p = abs(math.log(abs(e32 / e21))) / math.log(r21)
    for _ in range(maxit):
        q = math.log((r21**p - s) / (r32**p - s))
        pn = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(pn - p) < tol: p = pn; break
        p = pn
    return p, ("oscillatory" if s < 0 else "monotone")


def richardson(f_fine, f_med, p, r21):
    return f_fine + (f_fine - f_med) / (r21**p - 1.0)


def gci(f_fine, f_med, p, r21, Fs=1.25):
    ea = abs((f_fine - f_med) / f_fine) if f_fine != 0 else float("nan")
    return Fs * ea / (r21**p - 1.0)


def analyse(f, h, label="", p_min=0.5, p_max=3.0, Fs=1.25):
    """f, h coarse->fine, 3 or more levels.  Uses the finest three for the order; also reports every consecutive triple."""
    n = len(f); out = dict(label=label, levels=n, values=list(f), h=list(h), triples=[])
    if n < 3: out.update(accepted=False, reason="fewer than 3 levels"); return out
    for k in range(n - 2):
        p, st = observed_order(f[k:k + 3], h[k:k + 3]); out["triples"].append(dict(idx=(k, k + 1, k + 2), p=p, status=st))
    f3, f2, f1 = f[-3:]; h3, h2, h1 = h[-3:]
    p, st = out["triples"][-1]["p"], out["triples"][-1]["status"]; r21 = h2 / h1
    d21, d32 = f2 - f1, f3 - f2
    monotone = (d21 * d32 > 0) and abs(d21) < abs(d32)
    in_range = (not math.isnan(p)) and p_min <= p <= p_max
    out.update(p=p, status=st, monotone=bool(monotone), p_in_range=bool(in_range), r21=r21)
    if monotone and in_range:
        fe = richardson(f1, f2, p, r21); out.update(accepted=True, f_extrap=fe, gci_fine=gci(f1, f2, p, r21, Fs), rel_diff_fine_vs_extrap=(f1 - fe) / fe if fe else float("nan"),
                                                    asymptotic_ratio=None)
        # asymptotic range check using the two GCIs on the coarser pair (should be ~1)
        g23 = gci(f2, f3, p, h3 / h2, Fs)
        out["asymptotic_ratio"] = g23 / (r21**p * out["gci_fine"]) if out["gci_fine"] else float("nan")
    else:
        out.update(accepted=False, reason="; ".join(x for x in ("not monotone" if not monotone else "", "p=%.3g outside [%.1f,%.1f]" % (p, p_min, p_max) if not in_range else "") if x))
    return out


def format_analysis(a):
    if not a.get("accepted"): return "%s: NOT ACCEPTED (%s); p per triple = %s" % (a["label"], a.get("reason"), ", ".join("%.3g(%s)" % (t["p"], t["status"]) for t in a["triples"]))
    return "%s: p=%.3f (finest triple), f_extrap=%.6g, GCI_fine=%.3g%%, rel(fine-extrap)=%.3g%%, asymptotic ratio=%.3f" % (
        a["label"], a["p"], a["f_extrap"], 100 * a["gci_fine"], 100 * a["rel_diff_fine_vs_extrap"], a["asymptotic_ratio"])
