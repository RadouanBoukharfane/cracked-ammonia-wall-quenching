"""TODO item 7: 'planar HOQ at the local mixture' model tested quench event by quench event (Schneider I Figs. 10/12).

  python item7_local_model.py  -> item7_local_model.csv, figs/item7_local_model.pdf

Curve Phi_q,1D(phi_loc) per gamma from the planar HOQ at the nominal mixture and at the fresh mixtures reproducing the
q10, q50, q90, mean and p95 crest mixtures (PhaseB/runs/HOQ_crest_*): peak flux in W/m2 divided by the q_ref of the
nominal mixture, against the target phi_loc on the T_track isotherm. Each local event of analysis/local_events.csv
(phi_loc of the crest tip 1 tau_F before its flux maximum) gets the interpolated planar prediction.
"""
import glob, json, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "paper"))
from make_figs import W, INK2, COL, MK, tag, PLANAR_MIX  # noqa: E402
R = os.path.join(ROOT, "PhaseB", "runs")

def peak(run):
    m = json.load(open(os.path.join(run, "case.json"))); d = pd.read_csv(os.path.join(run, "hoq_timeseries.csv"))
    d = d[d.t_tauF > 0.5 * m.get("dist_lF", 30.0)]
    z = np.load(os.path.join(run, "hoq_wallflux.npz"))
    return d.Phi_max_star.max() * float(z["q_ref"]), m

rows, curves = [], {}
for g in (0.2, 0.43, 0.75, 1.0):
    nom = os.path.join(R, f"HOQ_g{g:.2f}_phi0.50_P10.0_T750_planar")
    q0, _ = peak(nom); qref = float(np.load(os.path.join(nom, "hoq_wallflux.npz"))["q_ref"])
    pts = [(PLANAR_MIX[g][0], q0 / qref, "nominal")]
    for r in sorted(glob.glob(os.path.join(R, f"HOQ_crest_g{g:.2f}_phi0.50_P10.0_T750_*_planar"))):
        q, m = peak(r); pts.append((m["phi_loc_target"], q / qref, os.path.basename(r).split("_")[-2]))
    pts.sort(); curves[g] = pts
    print(f"gamma {g}: " + ", ".join(f"{lab} phi_loc {p:.4f} -> Phi/Phi_nom {f / pts[[x[2] for x in pts].index('nominal')][1]:.3f}" for p, f, lab in pts))
SUF = sys.argv[1] if len(sys.argv) > 1 else ""
ev = pd.read_csv(os.path.join(ROOT, "analysis", f"local_events{SUF}.csv"))
ev = ev[(ev.Tw == 750) & ev.phi_loc.notna()].copy()
ev["Phi_pred"] = [np.interp(p, [a for a, _, _ in curves[g]], [b for _, b, _ in curves[g]]) for p, g in zip(ev.phi_loc, ev.gamma)]
ev["extrap"] = [not (curves[g][0][0] <= p <= curves[g][-1][0]) for p, g in zip(ev.phi_loc, ev.gamma)]
ev["ratio"] = ev.Phi_q / ev.Phi_pred
ev.to_csv(os.path.join(ROOT, "analysis", f"item7_local_model{SUF}.csv"), index=False)
for g, q in ev.groupby("gamma"):
    print(f"gamma {g}: n {len(q)}, inside curve range {100 * (~q.extrap).mean():.0f}% | measured/predicted median {q.ratio.median():.3f} "
          f"IQR {q.ratio.quantile(.25):.3f}-{q.ratio.quantile(.75):.3f} | top-10% events: {q[q.Phi_q > q.Phi_q.quantile(.9)].ratio.median():.3f}"
          f" | corr(Phi_q, Phi_pred) {q[['Phi_q', 'Phi_pred']].corr().iloc[0, 1]:+.2f}")
fig, axs = plt.subplots(1, 2, figsize=(W, 2.9), gridspec_kw=dict(wspace=0.3, left=0.09, right=0.98, bottom=0.17, top=0.93))
for g, pts in curves.items():
    q = ev[ev.gamma == g]
    axs[0].scatter(q.phi_loc, q.Phi_q, s=2, color=COL[g], alpha=0.2, lw=0, rasterized=True)
    axs[0].plot([a for a, _, _ in pts], [b for _, b, _ in pts], MK[g] + "-", color=COL[g], mec="white", ms=6)
    v = np.sort(q.ratio.values); axs[1].plot(v, np.linspace(0, 1, len(v)), color=COL[g], lw=1.6, label=rf"$\gamma = {g:g}$")
axs[0].set_xlabel(r"$\phi_{loc}$ of the crest before contact"); axs[0].set_ylabel(r"$\Phi_q^*$ (nominal normalisation)")
axs[1].axvline(1, color=INK2, lw=0.8, ls=":"); axs[1].set_xlabel(r"measured / planar HOQ at $\phi_{loc}$")
axs[1].set_ylabel("cumulative fraction"); axs[1].legend(frameon=False, fontsize=7.5); axs[1].set_xlim(0.4, 1.8)
tag(axs[0], 0); tag(axs[1], 1)
fig.savefig(os.path.join(ROOT, "analysis", "figs", f"item7_local_model{SUF}.pdf"), dpi=300)
