"""TODO item 5: stability database of the Phase A dispersion relations (waves 1 and 2).

  python stability_database.py   -> stability_database.csv, figs/stability_vs_X.pdf, figs/stability_maps.pdf

For every Phase A case with a dispersion relation (the CFL-0.3 rerun is used when it exists; runs flagged noisy without
a rerun are excluded): omega_max, k at the maximum, cut-off k_c; A0 properties of the unstretched flame (Ze, Le_eff on the
unburnt and burnt sides, diffusivity-weighted, sigma, l_D/l_F). Correlation of omega_max with X_u = Ze(1 - Le_eff,u) and
X_b = Ze(1 - Le_eff,b) (linear fits, R^2), and maps omega_max(gamma, phi) for each (P, T_u).
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "paper"))
from make_figs import W, INK, INK2, MUTED, HALO, tag  # noqa: E402  (same figure standard as the paper)

OUT = os.path.join(ROOT, "analysis")


def complete(run):
    """Linear window closed (max A > 50 A0) or run reached amr.stop_time."""
    import re
    m = json.load(open(os.path.join(run, "case.json")))
    fr = np.load(os.path.join(run, "front.npz"))
    yf = fr["yf"]
    A = 2 * np.abs(np.fft.rfft(yf - np.nanmean(yf, 1, keepdims=True), axis=1))[:, 1:m["nmodes"] + 1] / yf.shape[1]
    if (np.nanmax(A, 1) > 50 * m["A0"]).any():
        return True
    stop = float(re.search(r"amr.stop_time = (\S+)", open(os.path.join(run, "input.2d")).read()).group(1))
    return fr["t"][-1] >= 0.98 * stop


def collect():
    c = pd.read_csv(os.path.join(ROOT, "PhaseA", "cases", "phaseA_cases.csv"))
    rows = []
    for _, r in c.iterrows():
        base = os.path.join(ROOT, "PhaseA", "runs", r["name"])
        run = base + "_cfl03" if os.path.exists(base + "_cfl03/dispersion.csv") else base
        f = os.path.join(run, "dispersion.csv")
        if not os.path.exists(f):
            continue
        q = os.path.join(run, "quality.json")
        noisy = json.load(open(q)).get("noisy", False) if os.path.exists(q) else False
        if noisy:
            continue
        if not complete(run):          # run still going and linear window not closed: provisional dispersion
            continue
        d = pd.read_csv(f)
        d = d[(d.r2 > 0.9) | (d.omega_bar > 0)]
        i = d.omega_bar.idxmax()
        if d["mode"][i] >= d["mode"].max() - 2:      # maximum at the highest resolved modes: numerical, not physical
            print("excluded (peak at the k edge):", r["name"])
            continue
        pos = d[d.omega_bar > 0]
        a = json.load(open(os.path.join(ROOT, "A0", "results", "mei2021", "cases", r["name"][2:] + ".json")))
        rows.append(dict(name=r["name"], wave=r["wave"], gamma=r["gamma"], phi=r["phi"], P=r["P_atm"], Tu=r["T_u"],
                         omega_max=d.omega_bar[i], k_max=d.k_bar[i], k_c=pos.k_bar.max() if len(pos) else np.nan,
                         Ze=a["Ze"], Le_u=a["Le_eff"]["u"]["diff"], Le_b=a["Le_eff"]["b"]["diff"], sigma=a["sigma"],
                         lD_lF=a["l_D"] / a["delta_T"], sL=a["sL"], cfl03=run.endswith("_cfl03")))
    d = pd.DataFrame(rows)
    d["X_u"] = d.Ze * (1 - d.Le_u)
    d["X_b"] = d.Ze * (1 - d.Le_b)
    return d


def fit(x, y):
    b, a = np.polyfit(x, y, 1)
    r2 = 1 - ((y - (a + b * x))**2).sum() / ((y - y.mean())**2).sum()
    return a, b, r2


def main():
    d = collect()
    d.to_csv(os.path.join(OUT, "stability_database.csv"), index=False)
    print(f"{len(d)} cases (wave 1: {(d.wave == 1).sum()}, wave 2: {(d.wave == 2).sum()})")
    for xc in ("X_u", "X_b"):
        a, b, r2 = fit(d[xc].values, d.omega_max.values)
        print(f"omega_max = {a:.3f} + {b:.4f} {xc}   R2 = {r2:.3f}")
    # multiple regression with the hydrodynamic part (sigma)
    A = np.column_stack([np.ones(len(d)), d.X_b, d.sigma])
    coef, *_ = np.linalg.lstsq(A, d.omega_max, rcond=None)
    pred = A @ coef
    r2 = 1 - ((d.omega_max - pred)**2).sum() / ((d.omega_max - d.omega_max.mean())**2).sum()
    print(f"omega_max = {coef[0]:.3f} + {coef[1]:.4f} X_b + {coef[2]:.4f} sigma   R2 = {r2:.3f}")

    # ---- figure 1: omega_max against X_u and X_b
    mk = {(1.0, 500.0): "o", (1.0, 750.0): "s", (5.0, 500.0): "^", (5.0, 750.0): "D", (10.0, 500.0): "v",
          (10.0, 750.0): "P"}
    cmap = plt.get_cmap("viridis")
    fig = plt.figure(figsize=(W, 3.1))
    gs = fig.add_gridspec(1, 2, wspace=0.28, left=0.085, right=0.88, bottom=0.15, top=0.9)
    for j, (xc, xl) in enumerate([("X_u", r"$Ze\,(1 - Le_{eff})$, unburnt side"),
                                  ("X_b", r"$Ze\,(1 - Le_{eff})$, burnt side")]):
        ax = fig.add_subplot(gs[j])
        for key, q in d.groupby(["P", "Tu"]):
            sc = ax.scatter(q[xc], q.omega_max, c=q.gamma, cmap=cmap, vmin=0, vmax=1, marker=mk.get(key, "o"), s=28,
                            edgecolor="white", linewidth=0.4, label=rf"{key[0]:g} atm, {key[1]:.0f} K")
        a, b, r2 = fit(d[xc].values, d.omega_max.values)
        xx = np.linspace(d[xc].min(), d[xc].max(), 10)
        ax.plot(xx, a + b * xx, color=INK2, lw=1, ls=(0, (4, 2)))
        ax.text(0.03, 0.95, rf"linear fit, $R^2 = {r2:.2f}$", transform=ax.transAxes, va="top", fontsize=8, color=INK2)
        ax.set_xlabel(xl); ax.set_ylabel(r"$\bar\omega_{\max}$")
        tag(ax, j, x=-0.16)
    cax = fig.add_axes([0.9, 0.2, 0.015, 0.65])
    cb = fig.colorbar(sc, cax=cax); cb.set_label(r"cracking ratio $\gamma$")
    present = sorted(set(zip(d.P, d.Tu)))
    fig.legend([plt.Line2D([], [], marker=mk.get(k, "o"), color=INK2, ls="", ms=6) for k in present],
               [rf"{k[0]:g} atm, {k[1]:.0f} K" for k in present], loc="upper center", ncol=6, frameon=False,
               fontsize=7.5, bbox_to_anchor=(0.48, 1.0))
    fig.savefig(os.path.join(OUT, "figs", "stability_vs_X.pdf")); plt.close(fig)

    # ---- figure 2: maps omega_max(gamma, phi) per (P, T_u)
    Ps, Ts = sorted(d.P.unique()), sorted(d.Tu.unique())
    fig, axs = plt.subplots(len(Ts), len(Ps), figsize=(W, 1.9 * len(Ts) + 0.6), squeeze=False, sharex=True, sharey=True)
    vmax = d.omega_max.max()
    for a_, T in enumerate(Ts):
        for b_, P in enumerate(Ps):
            ax = axs[a_, b_]
            q = d[(d.P == P) & (d.Tu == T)]
            if len(q):
                pv = q.pivot_table(index="phi", columns="gamma", values="omega_max")
                im = ax.imshow(pv.values, origin="lower", aspect="auto", cmap="magma", vmin=0, vmax=vmax,
                               extent=[-0.5, len(pv.columns) - 0.5, -0.5, len(pv.index) - 0.5])
                ax.set_xticks(range(len(pv.columns)), [f"{g:g}" for g in pv.columns], fontsize=6.5)
                ax.set_yticks(range(len(pv.index)), [f"{p:g}" for p in pv.index], fontsize=7)
                for (ii, ph) in enumerate(pv.index):
                    row = pv.loc[ph]
                    if row.notna().any():
                        jj = int(np.nanargmax(row.values))
                        ax.plot(jj, ii, "o", ms=3.5, mfc="none", mec="white", mew=0.8)
            ax.set_title(rf"{P:g} atm, $T_u = {T:.0f}$ K", fontsize=8)
            ax.grid(False)
            if a_ == len(Ts) - 1:
                ax.set_xlabel(r"$\gamma$")
            if b_ == 0:
                ax.set_ylabel(r"$\phi$")
    fig.subplots_adjust(left=0.08, right=0.88, bottom=0.12, top=0.93, wspace=0.08, hspace=0.3)
    cax = fig.add_axes([0.9, 0.2, 0.015, 0.65])
    cb = fig.colorbar(im, cax=cax); cb.set_label(r"$\bar\omega_{\max}$ (circles: maximum in $\gamma$)")
    fig.savefig(os.path.join(OUT, "figs", "stability_maps.pdf")); plt.close(fig)


if __name__ == "__main__":
    main()
