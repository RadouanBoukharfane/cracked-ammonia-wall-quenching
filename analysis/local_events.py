"""TODO item 1: statistics of the local quenching events of the unstable HOQ (Schneider I Fig. 8, Traut Figs. 6-7).

  python local_events.py NPROC [LEAD]  -> local_events[_leadX].csv, figs/local_events[_leadX].pdf (LEAD in tau_F, default 1)

For every unstable run and every STRIDE-th wall column x (hoq_wallflux.npz with per-column y_hrr from post_hoq.py):
  Phi_q(x) = max_t Phi_w*(x, t) (t > 5 tau_F), t*(x) its instant; Pe_D(x) = y_HRR(x, t*)/l_D (diffusive thickness);
  curvature, phi_loc, gamma_loc of the crest tip above x on the T_track isotherm LEAD tau_F before t*(x)
  (PhaseB/crest_at_quench.tip). Normalised by the planar HOQ of the same gamma and T_w (Phi_q and Pe_D at its peak).
"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "paper"))
sys.path.insert(0, os.path.join(ROOT, "PhaseB"))
from make_figs import W, INK, INK2, HALO, COL, MK, tag  # noqa: E402
from crest_at_quench import tip, plt_list  # noqa: E402

OUT = os.path.join(ROOT, "analysis")
STRIDE = 8
LEAD = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
SUF = "" if LEAD == 1.0 else f"_lead{LEAD:g}"
A0 = os.path.join(ROOT, "A0", "results", "mei2021", "cases")


def lD_of(m):
    c = m["case"]
    return json.load(open(os.path.join(A0, f"g{c['gamma']:.2f}_phi{c['phi']:.2f}_P{c['P_atm']:.1f}_T{c['T_u']:.0f}.json")))["l_D"]


def planar_ref(g, tw):
    r = os.path.join(ROOT, "PhaseB", "runs", f"HOQ_g{g:.2f}_phi0.50_P10.0_T750_planar" + ("" if tw == 750 else f"_Tw{tw:.0f}"))
    m = json.load(open(os.path.join(r, "case.json")))
    z = np.load(os.path.join(r, "hoq_wallflux.npz"))
    t = z["t"] / m["tau_F"]; o = np.argsort(t)
    phi = z["phi_w"][o].mean(1) / z["q_ref"]; yh = z["y_hrr"][o].mean(1)
    k = t[o] > 0.5 * m.get("dist_lF", 30.0)
    i = np.where(k)[0][np.argmax(phi[k])]
    return phi[i], yh[i] / lD_of(m)


def main(nproc):
    from multiprocessing import Pool
    runs = [r for r in sorted(glob.glob(os.path.join(ROOT, "PhaseB", "runs", "HOQ_g*_t??.0*")))
            if os.path.exists(os.path.join(r, "hoq_wallflux.npz")) and "y_hrr" in np.load(os.path.join(r, "hoq_wallflux.npz"))]
    jobs, meta = [], {}
    for run in runs:
        m = json.load(open(os.path.join(run, "case.json")))
        z = np.load(os.path.join(run, "hoq_wallflux.npz"))
        t = z["t"] / m["tau_F"]; o = np.argsort(t); t = t[o]
        phi = z["phi_w"][o] / z["q_ref"]; yh = z["y_hrr"][o]
        k = t > 5.0
        cols = np.arange(0, phi.shape[1], STRIDE)
        iq = np.argmax(phi[k][:, cols], axis=0)
        tq = t[k][iq]
        Phi = phi[k][iq, cols]
        Pe = yh[k][iq, cols] / lD_of(m)
        xs = (cols + 0.5) * m["Lx"] / phi.shape[1]
        plts = plt_list(run, m["tau_F"]); tp = np.array([a for a, _ in plts])
        req = {}
        for j, (tt, x) in enumerate(zip(tq, xs)):
            p = plts[int(np.argmin(np.abs(tp - (tt - LEAD))))][1]
            req.setdefault(p, []).append(j)
        meta[run] = (m, xs, tq, Phi, Pe, req)
        jobs += [(p, m, [xs[j] for j in jj]) for p, jj in req.items()]
    with Pool(nproc) as pool:
        res = dict(pool.map(tip, jobs))
    rows = []
    for run, (m, xs, tq, Phi, Pe, req) in meta.items():
        g, tw = round(m["case"]["gamma"], 2), m["T_wall"]
        Pp, Pep = planar_ref(g, tw)
        kap = np.full(len(xs), np.nan); ph = kap.copy(); ga = kap.copy()
        for p, jj in req.items():
            for j, (y_, k_, f_, g_) in zip(jj, res[p]):
                kap[j], ph[j], ga[j] = k_, f_, g_
        for j in range(len(xs)):
            rows.append(dict(run=os.path.basename(run), gamma=g, Tw=tw, x_lF=xs[j] / m["l_F"], t_q=tq[j], Phi_q=Phi[j],
                             Pe_D=Pe[j], Phi_rel=Phi[j] / Pp, Pe_rel=Pe[j] / Pep, product=Phi[j] * Pe[j],
                             kappa_lF=kap[j], phi_loc=ph[j], gamma_loc=ga[j]))
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(OUT, f"local_events{SUF}.csv"), index=False)
    for (g, tw), q in d.groupby(["gamma", "Tw"]):
        print(f"gamma {g} Tw {tw:.0f}: n={len(q)}  Phi_rel median {q.Phi_rel.median():.2f} p90 {q.Phi_rel.quantile(.9):.2f} "
              f"| frac Phi_rel>1 {(q.Phi_rel > 1).mean():.2f} | Pe_rel median {q.Pe_rel.median():.2f} | Phi*Pe_D median "
              f"{q['product'].median():.2f} | corr(Phi_rel, kappa) {q[['Phi_rel', 'kappa_lF']].corr().iloc[0, 1]:+.2f}")

    fig = plt.figure(figsize=(W, 3.0))
    gs = fig.add_gridspec(1, 3, wspace=0.42, left=0.07, right=0.97, bottom=0.17, top=0.92)
    ax, bx, cx = (fig.add_subplot(gs[i]) for i in range(3))
    for g, q in d[d.Tw == 750].groupby("gamma"):
        for a_, c_ in ((ax, "Phi_rel"), (bx, "Pe_rel")):
            v = np.sort(q[c_].dropna().values)
            a_.plot(v, np.linspace(0, 1, len(v)), color=COL[g], lw=1.6, label=rf"$\gamma = {g:g}$")
    for a_, xl in ((ax, r"$\Phi_q/\Phi_{q,\mathrm{planar}}$"), (bx, r"$\mathrm{Pe}_q/\mathrm{Pe}_{q,\mathrm{planar}}$")):
        a_.axvline(1, color=INK2, lw=0.7, ls=(0, (2, 2))); a_.set_xlabel(xl); a_.set_ylabel("cumulative fraction of wall")
    ax.legend(frameon=False, fontsize=7.5)
    tag(ax, 0, x=-0.2); tag(bx, 1, x=-0.2)
    q = d[(d.Tw == 750)]
    sc = cx.scatter(q.Pe_rel, q.Phi_rel, c=q.kappa_lF, cmap="RdBu_r", vmin=-0.6, vmax=0.6, s=4, lw=0, rasterized=True)
    xx = np.linspace(0.4, 1.6, 50)
    cx.plot(xx, 1 / xx, color=INK, lw=1, ls=(0, (3, 2)))
    cx.text(1.25, 1 / 1.25 + 0.06, r"$\Phi_q \propto 1/\mathrm{Pe}_q$", fontsize=8, color=INK)
    cx.set_xlabel(r"$\mathrm{Pe}_q/\mathrm{Pe}_{q,\mathrm{planar}}$"); cx.set_ylabel(r"$\Phi_q/\Phi_{q,\mathrm{planar}}$")
    cb = fig.colorbar(sc, ax=cx, fraction=0.05, pad=0.02); cb.set_label(r"$\kappa\,l_F$ before contact")
    tag(cx, 2, x=-0.22)
    fig.savefig(os.path.join(OUT, "figs", f"local_events{SUF}.pdf"), dpi=300); plt.close(fig)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
