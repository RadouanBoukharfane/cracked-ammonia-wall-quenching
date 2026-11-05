"""TODO item 6: hydrodynamic theory of Matalon et al. (coefficients of Altantzis et al., as in Lehmann et al. 2025,
supplementary material, Eqs. 1-3) against the computed dispersion relations, and an effective numerical coefficient.

  python matalon.py   (needs stability_database.csv)  -> matalon.csv, figs/matalon.pdf

Theory, in units of l_D = D_th/s_L and l_D/s_L:  w = w0 k - Omega k^2,  Omega = B1 + Ze (Le_eff - 1) B2 + Pr (sigma - 1) B3,
with w0 = (sqrt(sigma^3 + sigma^2 - sigma) - sigma)/(sigma + 1) and B1, B2, B3 integrals of lambda(x)/lambda_u over
x = T/T_u in [1, sigma], evaluated from the A0 one-dimensional profiles (lambda from Cantera along the profile).
In the Phase A units (l_F, tau_F = l_F/s_L):  w_bar = w0 k_bar - Omega (l_D/l_F) k_bar^2, w_bar_max = w0^2/(4 Omega l_D/l_F).
Omega is computed with Le_eff from both sides (unburnt, burnt; diffusivity-weighted). The numerical Omega_num is fitted
to the computed w_bar/k_bar = w0 - Omega_num (l_D/l_F) k_bar for k_bar <= k_max.
"""
import json
import os

import numpy as np
import pandas as pd
import cantera as ct
from scipy.integrate import quad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(ROOT, "paper"))
from make_figs import W, INK, INK2, HALO, tag  # noqa: E402

MECH = os.path.join(ROOT, "mechanisms", "mei2021", "mei2021.yaml")
OUT = os.path.join(ROOT, "analysis")


def coefficients(name, gas):
    a = json.load(open(os.path.join(ROOT, "A0", "results", "mei2021", "cases", name + ".json")))
    p = pd.read_csv(os.path.join(ROOT, "A0", "results", "mei2021", "profiles", name + ".csv"))
    P = a["P_atm"] * ct.one_atm
    sp = [c[2:] for c in p.columns if c.startswith("X_") and c[2:] in gas.species_names]
    lam = []
    for _, row in p.iterrows():
        X = {s: max(row["X_" + s], 0.0) for s in sp}
        gas.TPX = row["T"], P, X
        lam.append(gas.thermal_conductivity)
    x_prof = p["T"].values / a["Tu"]
    lt = np.array(lam) / lam[0]
    o = np.argsort(x_prof)
    x_prof, lt = x_prof[o], lt[o]
    sig = a["sigma"]
    lam_t = lambda x: np.interp(x, x_prof, lt)                      # constant extrapolation beyond the profile
    w0 = (np.sqrt(sig**3 + sig**2 - sig) - sig) / (sig + 1)
    I1 = quad(lambda x: lam_t(x) / x, 1, sig, limit=200)[0]
    I0 = quad(lam_t, 1, sig, limit=200)[0]
    I2 = quad(lambda x: np.log((sig - 1) / (x - 1)) * lam_t(x) / x, 1, sig, limit=400, points=[1 + 1e-6])[0]
    den = sig + (sig + 1) * w0
    B1 = sig / 2 * (sig * (2 * w0 + sig + 1) / ((sig - 1) * den) * I1 + I0 / den)
    B2 = sig / 2 * ((1 + w0) * (sig + w0) / ((sig - 1) * den) * I2)
    B3 = sig / 2 * (2 * (sig - 1) / den * lam_t(sig) - 2 / den * I0)
    return dict(w0=w0, B1=B1, B2=B2, B3=B3, Pr=a["Pr_u"], Ze=a["Ze"], sigma=sig)


def main():
    gas = ct.Solution(MECH)
    d = pd.read_csv(os.path.join(OUT, "stability_database.csv"))
    rows = []
    for _, r in d.iterrows():
        name = r["name"][2:]
        c = coefficients(name, gas)
        out = dict(name=r["name"], gamma=r.gamma, phi=r.phi, P=r.P, Tu=r.Tu, omega_max=r.omega_max, lD_lF=r.lD_lF, **c)
        for side, le in (("u", r.Le_u), ("b", r.Le_b)):
            Om = c["B1"] + c["Ze"] * (le - 1) * c["B2"] + c["Pr"] * (c["sigma"] - 1) * c["B3"]
            out[f"Omega_{side}"] = Om
            out[f"wmax_th_{side}"] = c["w0"]**2 / (4 * Om * r.lD_lF) if Om > 0 else np.inf
        base = os.path.join(ROOT, "PhaseA", "runs", r["name"])
        run = base + "_cfl03" if os.path.exists(base + "_cfl03/dispersion.csv") else base
        dd = pd.read_csv(os.path.join(run, "dispersion.csv"))
        dd = dd[(dd.k_bar <= r.k_max) & (dd.omega_bar > 0)]
        if len(dd) >= 3:
            s, i0 = np.polyfit(dd.k_bar, dd.omega_bar / dd.k_bar, 1)
            out["Omega_num"] = -s / r.lD_lF
            out["w0_num"] = i0
        rows.append(out)
        print(f"{r['name']:30s} w_max {r.omega_max:.3f}  theory(u) {out['wmax_th_u']:.3f}  theory(b) "
              f"{out['wmax_th_b']:.3f}  Omega_num {out.get('Omega_num', np.nan):.2f} Omega_u {out['Omega_u']:.2f} "
              f"Omega_b {out['Omega_b']:.2f}", flush=True)
    m = pd.DataFrame(rows)
    m.to_csv(os.path.join(OUT, "matalon.csv"), index=False)
    ok = np.isfinite(m.wmax_th_b) & np.isfinite(m.wmax_th_u)
    for s in ("u", "b"):
        e = (m[f"wmax_th_{s}"][ok] / m.omega_max[ok] - 1)
        print(f"theory ({s}) / numerics - 1: median {100 * e.median():+.0f} %, IQR {100 * e.quantile(0.25):+.0f} to "
              f"{100 * e.quantile(0.75):+.0f} %")

    fig = plt.figure(figsize=(W, 3.1))
    gs = fig.add_gridspec(1, 2, wspace=0.3, left=0.085, right=0.88, bottom=0.15, top=0.92)
    ax, bx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    lim = [0, max(1.2, m.omega_max.max() * 1.1)]
    for s, mk_, lab in (("u", "o", "unburnt-side $Le_{eff}$"), ("b", "s", "burnt-side $Le_{eff}$")):
        sc = ax.scatter(m.omega_max[ok], m[f"wmax_th_{s}"][ok], c=m.gamma[ok], cmap="viridis", vmin=0, vmax=1,
                        marker=mk_, s=22, edgecolor="white", linewidth=0.4, label=lab)
    ax.plot(lim, lim, color=INK2, lw=0.8)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel(r"computed $\bar\omega_{\max}$"); ax.set_ylabel(r"Matalon theory $\bar\omega_{\max}$")
    ax.legend(frameon=False, loc="upper left")
    tag(ax, 0, x=-0.15)
    from make_figs import COL as GC
    for g in (0.20, 0.43, 0.75, 1.00):
        nm = f"A_g{g:.2f}_phi0.50_P10.0_T750"
        row = m[m.name == nm]
        if not len(row):
            continue
        row = row.iloc[0]
        base = os.path.join(ROOT, "PhaseA", "runs", nm)
        run = base + "_cfl03" if os.path.exists(base + "_cfl03/dispersion.csv") else base
        dd = pd.read_csv(os.path.join(run, "dispersion.csv"))
        dd = dd[(dd.r2 > 0.9) | (dd.omega_bar > 0)]
        kk = np.linspace(0, 3, 200)
        bx.plot(dd.k_bar, dd.omega_bar, color=GC[g], lw=1.6)
        bx.plot(kk, row.w0 * kk - row.Omega_u * row.lD_lF * kk**2, color=GC[g], lw=1.0, ls=(0, (3, 2)))
        bx.text(dd.k_bar[dd.omega_bar.idxmax()], dd.omega_bar.max() + 0.04, rf"{g:g}", color=GC[g], fontsize=8,
                ha="center", path_effects=HALO)
    bx.plot([], [], color=INK2, lw=1.6, label="computed"); bx.plot([], [], color=INK2, lw=1.0, ls=(0, (3, 2)),
                                                                   label="Matalon, unburnt-side $Le_{eff}$")
    bx.legend(frameon=False, loc="lower left", fontsize=7.5)
    bx.axhline(0, color=INK2, lw=0.6); bx.set_xlim(0, 3); bx.set_ylim(-0.6, 1.0)
    bx.set_xlabel(r"$\bar k$"); bx.set_ylabel(r"$\bar\omega$ ($\phi = 0.5$, 10 atm, 750 K)")
    tag(bx, 1, x=-0.17)
    cax = fig.add_axes([0.9, 0.2, 0.015, 0.65])
    cb = fig.colorbar(sc, cax=cax); cb.set_label(r"cracking ratio $\gamma$")
    fig.savefig(os.path.join(OUT, "figs", "matalon.pdf")); plt.close(fig)


if __name__ == "__main__":
    main()
