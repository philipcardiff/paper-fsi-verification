"""Figures from the CSVs in ../data (numpy not needed; matplotlib).  usage: python3 -B -I scripts/make_figures.py"""
import csv, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "..", "data"); F = os.path.join(HERE, "..", "figures")
rd = lambda n: list(csv.DictReader(open(os.path.join(D, n))))
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
META = {"CreationDate": None, "Creator": None, "Producer": None}
fac = rd("factorial_V1.csv"); ps = rd("time_limit_p_sensitivity.csv"); off = rd("same_dt_offset_to_oomph.csv"); gap = rd("mf1_gap.csv")
col = {"S1": "#8c8c8c", "S2": "#1f4e79", "S3": "#b05a00"}; mk = {"S1": "s", "S2": "o", "S3": "^"}
# Figure 1: trough value at x = 10 against dt, and the same-dt offset to oomph
fig, ax = plt.subplots(1, 2, figsize=(6.8, 2.7))
a = ax[0]
for S in ("S1", "S2", "S3"):
    r = [q for q in fac if q["x_m"] == "10" and q["mesh"] == S]
    a.plot([float(q["dt_s"]) for q in r], [float(q["FE_V1_pctP"]) for q in r], mk[S] + "-", color=col[S], ms=4, lw=1, label="FE " + S)
r = [q for q in fac if q["x_m"] == "10" and q["mesh"] == "S2"]
a.plot([float(q["dt_s"]) for q in r], [float(q["oomph_same_dt_V1_pctP"]) for q in r], "x--", color="k", ms=5, lw=0.8, label="oomph, same dt")
a.axhline(float(r[0]["oomph_reference_V1_pctP"]), color="k", lw=0.8, ls=":", label="oomph reference (U_ref 0.30 % P)")
lo = min(float(q["limit_V1_pctP"]) for q in ps if q["mesh"] == "S2" and q["x_m"] == "10"); hi = max(float(q["limit_V1_pctP"]) for q in ps if q["mesh"] == "S2" and q["x_m"] == "10")
a.axhspan(lo, hi, color="#1f4e79", alpha=0.15, lw=0, label="FE S2 time limit, p-band")
a.set_xscale("log"); a.set_xticks([0.025, 0.0125, 0.00625, 0.003125]); a.set_xticklabels(["0.025", "0.0125", "0.00625", "0.003125"], rotation=30); a.minorticks_off(); a.invert_xaxis(); a.set_xlabel(r"$\Delta t$ (s)"); a.set_ylabel(r"first trough, $x=10$ (% of $P$)"); a.set_ylim(-101, -98.3)
a.legend(fontsize=6.5, frameon=False, loc="lower left")
b = ax[1]
for S in ("S1", "S2", "S3"):
    r = [q for q in off if q["x_m"] == "10" and q["mesh"] == S]
    b.plot([float(q["dt_s"]) for q in r], [float(q["FE_minus_oomph_same_dt_V1_pctP"]) for q in r], mk[S] + "-", color=col[S], ms=4, lw=1, label="FE %s minus oomph" % S)
b.set_xscale("log"); b.set_xticks([0.025, 0.0125, 0.00625, 0.003125]); b.set_xticklabels(["0.025", "0.0125", "0.00625", "0.003125"], rotation=30); b.minorticks_off(); b.invert_xaxis(); b.set_xlabel(r"$\Delta t$ (s)"); b.set_ylabel(r"same-$\Delta t$ trough difference (% of $P$)"); b.set_ylim(-0.25, 0.1); b.axhline(0, color="k", lw=0.5)
b.legend(fontsize=7, frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(F, "fenicsx_collapsible_trough_time.pdf"), metadata=META); fig.savefig(os.path.join(F, "fenicsx_collapsible_trough_time.png"), dpi=150)
# Figure 2: open-loop beam-minus-continuum gap (MF1) against the FE-minus-oomph trough offset, per x
fig, ax = plt.subplots(figsize=(4.4, 2.8))
xs = ["7.5", "10", "12.5"]; w = 0.26
mf = {q["x_m"] if q["x_m"] != "10.0" else "10": q for q in gap if q["beam_variant"] == "ne160_d2dt2"}
cont_minus_beam = [-float(mf[x]["dV1_beam_minus_continuum_pctP"]) for x in xs]
o2 = [float([q for q in off if q["mesh"] == "S2" and q["x_m"] == x and q["dt_s"] == "0.003125"][0]["FE_minus_oomph_same_dt_V1_pctP"]) for x in xs]
o3 = [float([q for q in off if q["mesh"] == "S3" and q["x_m"] == x and q["dt_s"] == "0.003125"][0]["FE_minus_oomph_same_dt_V1_pctP"]) for x in xs]
ax.bar([i - w for i in range(3)], o2, w, color=col["S2"], label="FE S2 minus oomph (same dt 0.003125)")
ax.bar([i for i in range(3)], o3, w, color=col["S3"], label="FE S3 minus oomph (same dt 0.003125)")
ax.bar([i + w for i in range(3)], cont_minus_beam, w, color="#4d7d2f", label="open loop: continuum minus beam (MF1)")
ax.axhline(0, color="k", lw=0.5); ax.set_xticks(range(3)); ax.set_xticklabels(["x = 7.5", "x = 10", "x = 12.5"]); ax.set_ylabel(r"first trough difference (% of $P$)")
h_, l_ = ax.get_legend_handles_labels(); fig.legend(h_, l_, fontsize=6.5, frameon=False, loc="lower center"); fig.tight_layout(rect=(0, 0.17, 1, 1))
fig.savefig(os.path.join(F, "fenicsx_collapsible_model_form.pdf"), metadata=META); fig.savefig(os.path.join(F, "fenicsx_collapsible_model_form.png"), dpi=150)
print("ok")
