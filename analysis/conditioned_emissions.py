"""TODO item 4: near-wall emissions and NO/N2O reaction budgets of the unstable HOQ, conditioned on crests and troughs.

  python conditioned_emissions.py NPROC
     -> conditioned_emissions.csv, conditioned_budget_NO.csv, conditioned_budget_N2O.csv, figs/conditioned_emissions.pdf

Classification of the wall columns by the instant t*(x) of their local flux maximum (hoq_wallflux.npz): first third
reached = leading crests ("crest"), last third = troughs ("trough"), the rest "middle". t_Q as in PhaseB/emissions_tq.py
(hrr_total.csv). For each plotfile in [t_Q - 2, t_Q + 10] tau_F: per column, mole fractions averaged over 0 < y < 2 l_F
and the NO and N2O net production of every reaction integrated over 0 < y < 2 l_F; class averages; time integrals of the
budgets (kmol m^-2, duplicates summed). Emissions at t_Q + 5 and + 10 tau_F per class.
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
from make_figs import W, INK, INK2, HALO, COL, tag  # noqa: E402

OUT = os.path.join(ROOT, "analysis")
MECH = os.path.join(os.environ.get("PELE_HOME", "PeleLMeX"), "Submodules/PelePhysics/Mechanisms/Mei2021NH3/mechanism.yaml")
CLASSES = ("crest", "middle", "trough")
SP = ("NH3", "NO", "N2O")


def t_quench(run, tau):
    h = pd.read_csv(os.path.join(run, "hrr_total.csv"))
    tt = h.t.values / tau
    h0 = np.median(h.hrr_tot.values[:5]); k = np.argmax(h.hrr_tot.values < 0.5 * h0)
    a, b = h.iloc[k - 1], h.iloc[k]
    return tt[k - 1] + (a.hrr_tot - 0.5 * h0) / (a.hrr_tot - b.hrr_tot) * (tt[k] - tt[k - 1])


def one(args):
    p, m, cls = args
    import cantera as ct
    import yt
    yt.set_log_level(50)
    gas = ct.Solution(MECH)
    ds = yt.load(p)
    lev = ds.index.max_level
    n = (ds.domain_dimensions * ds.refine_by**lev).astype(int)
    dy = (ds.domain_width / n).d[1]
    ny = int(np.ceil(2.0 * m["l_F"] / dy))
    cg = ds.covering_grid(lev, left_edge=ds.domain_left_edge, dims=[n[0], ny, 1])
    T = cg["boxlib", "temp"].d[:, :, 0]
    Y = np.stack([np.clip(cg["boxlib", f"Y({s})"].d[:, :, 0], 0, None) for s in gas.species_names], -1)
    Y /= Y.sum(-1, keepdims=True)
    nx = T.shape[0]
    st = ct.SolutionArray(gas, nx * ny)
    st.TPY = T.ravel(), m["P_Pa"], Y.reshape(-1, Y.shape[-1])
    X = st.X.reshape(nx, ny, -1)
    q = st.net_rates_of_progress.reshape(nx, ny, -1)
    xs = (np.arange(nx) + 0.5) / nx                     # fraction of the domain width -> class of the column
    ci = cls[np.minimum((xs * len(cls)).astype(int), len(cls) - 1)]
    out = {"t": float(ds.current_time)}
    for s in SP:
        Xs = X[:, :, gas.species_index(s)].mean(axis=1)
        for c in range(3):
            out[f"X_{s}_{CLASSES[c]}"] = Xs[ci == c].mean()
    for s in ("NO", "N2O"):
        i = gas.species_index(s)
        nu = gas.product_stoich_coeffs[i, :] - gas.reactant_stoich_coeffs[i, :]
        w = (q * nu).sum(axis=1) * dy                    # (nx, nreac) kmol m^-2 s^-1
        for c in range(3):
            out[f"w{s}_{CLASSES[c]}"] = w[ci == c].mean(axis=0)
    return out


def main(nproc):
    from multiprocessing import Pool
    import cantera as ct
    gas = ct.Solution(MECH)
    eqs = [r.equation for r in gas.reactions()]
    runs = [r for r in sorted(glob.glob(os.path.join(ROOT, "PhaseB", "runs", "HOQ_g*_t??.0*")))
            if os.path.exists(os.path.join(r, "hrr_total.csv")) and os.path.exists(os.path.join(r, "hoq_wallflux.npz"))
            and "gamma" not in r and "1.00" not in os.path.basename(r)]
    tasks, meta = [], {}
    for run in runs:
        m = json.load(open(os.path.join(run, "case.json")))
        z = np.load(os.path.join(run, "hoq_wallflux.npz"))
        t = z["t"] / m["tau_F"]; o = np.argsort(t)
        phi = z["phi_w"][o]; k = t[o] > 5.0
        tstar = t[o][k][np.argmax(phi[k], axis=0)]
        q1, q2 = np.quantile(tstar, [1 / 3, 2 / 3])
        cls = np.where(tstar <= q1, 0, np.where(tstar >= q2, 2, 1))
        tQ = t_quench(run, m["tau_F"])
        plts = []
        for p in glob.glob(os.path.join(run, "plt*")):
            if not re.fullmatch(r"plt\d+", os.path.basename(p)):
                continue
            L = open(os.path.join(p, "Header")).read().split("\n")
            tt = float(L[2 + int(L[1]) + 1]) / m["tau_F"]
            if tQ - 2.01 <= tt <= tQ + 10.01:
                plts.append(p)
        meta[run] = (m, tQ, len(plts))
        tasks += [(p, m, cls) for p in plts]
    with Pool(nproc) as pool:
        res = pool.map(one, tasks)
    E, BN, BO = [], [], []
    i = 0
    for run in runs:
        m, tQ, n = meta[run]
        rr = sorted(res[i:i + n], key=lambda a: a["t"]); i += n
        tt = np.array([a["t"] for a in rr]) / m["tau_F"]
        ts = (tt - tQ) * m["tau_F"]
        g, tw = round(m["case"]["gamma"], 2), m["T_wall"]
        for c in CLASSES:
            row = dict(run=os.path.basename(run), gamma=g, Tw=tw, cls=c)
            for dt in (5.0, 10.0):
                j = int(np.argmin(np.abs(tt - (tQ + dt))))
                for s in SP:
                    row[f"{s}_{dt:g}"] = rr[j][f"X_{s}_{c}"]
            E.append(row)
            for s, B in (("NO", BN), ("N2O", BO)):
                wv = np.array([a[f"w{s}_{c}"] for a in rr])
                integ = pd.Series(np.trapezoid(wv, ts, axis=0) * 1e9, index=eqs).groupby(level=0).sum()
                B.append(pd.DataFrame(dict(run=os.path.basename(run), gamma=g, Tw=tw, cls=c, reaction=integ.index,
                                           value=integ.values)))
    E = pd.DataFrame(E)
    E.to_csv(os.path.join(OUT, "conditioned_emissions.csv"), index=False)
    pd.concat(BN).to_csv(os.path.join(OUT, "conditioned_budget_NO.csv"), index=False)
    pd.concat(BO).to_csv(os.path.join(OUT, "conditioned_budget_N2O.csv"), index=False)
    pl = pd.read_csv(os.path.join(ROOT, "PhaseB", "emissions_tq.csv")); pl["gamma"] = pl.gamma.round(2)
    pl = pl[~pl.run.str.contains("crest|mech")]          # nominal-mixture Mei runs only
    pl = pl[pl.unstable == 0]
    print("near-wall mole fractions at t_Q + 10 tau_F, class mean over realisations, relative to the planar HOQ")
    for (g, tw), q in E.groupby(["gamma", "Tw"]):
        ref = pl[(pl.gamma == g) & (pl.Tw == tw)].iloc[0]
        s = " | ".join(f"{c}: " + " ".join(f"{sp} x{q[q.cls == c][f'{sp}_10'].mean() / ref[f'{sp}_10']:.2f}" for sp in SP)
                       for c in ("crest", "trough"))
        print(f"  gamma {g} Tw {tw:.0f}: {s}")

    fig = plt.figure(figsize=(W, 3.0))
    gs = fig.add_gridspec(1, 3, wspace=0.38, left=0.08, right=0.98, bottom=0.2, top=0.9)
    for j, sp in enumerate(SP):
        ax = fig.add_subplot(gs[j])
        groups = [(0.2, 750.0), (0.43, 750.0), (0.75, 750.0), (0.43, 500.0)]
        for x, (g, tw) in enumerate(groups):
            q = E[(E.gamma == g) & (E.Tw == tw)]
            if not len(q):
                continue
            ref = pl[(pl.gamma == g) & (pl.Tw == tw)].iloc[0][f"{sp}_10"]
            for dx, c, mk in ((-0.15, "crest", "^"), (0.15, "trough", "v")):
                v = q[q.cls == c].groupby("run")[f"{sp}_10"].mean() / ref
                ax.plot([x + dx] * len(v), v, mk, color=COL[g], ms=4, alpha=0.45, mew=0)
                ax.plot(x + dx, v.mean(), mk, color=COL[g], ms=8, mec="white", mew=0.8)
        ax.axhline(1, color=INK2, lw=0.7); ax.set_yscale("log")
        ax.set_xticks(range(4), ["0.2", "0.43", "0.75", "0.43\n500 K"]); ax.set_xlabel(r"$\gamma$")
        ax.set_title(sp.replace("NH3", "NH$_3$").replace("N2O", "N$_2$O") + " near the wall, unstable/planar",
                     fontsize=8.5)
        tag(ax, j, x=-0.17)
    fig.legend([plt.Line2D([], [], marker="^", color=INK2, ls="", ms=7), plt.Line2D([], [], marker="v", color=INK2,
                ls="", ms=7)], ["beneath the first-quenching crests", "beneath the last-quenching troughs"],
               loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.savefig(os.path.join(OUT, "figs", "conditioned_emissions.pdf")); plt.close(fig)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
