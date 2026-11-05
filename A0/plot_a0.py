"""Figures of Phase A0 from results/<tag>/a0_summary.csv.

  python plot_a0.py --tag mei2021 [--compare mei2021_nosoret]

Small multiples: rows = phi, columns = T_u, x = gamma, one line per P
(sequential single-hue ramp, P is ordered). Cases with s_L < SL_MIN are
treated as near the flammability limit and drawn as open markers.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

P_COLORS = {1.0: "#86b6e0", 5.0: "#3a7dc0", 10.0: "#0b2f6b"}
SL_MIN = 0.03  # m/s

QUANTITIES = {
    "sL": ("$s_L$ [cm/s]", lambda d: d.sL * 100, "log"),
    "delta_T": (r"$\delta_T$ [µm]", lambda d: d.delta_T * 1e6, "log"),
    "Ze": ("Ze (N$_2$ dilution)", lambda d: d.Ze, "linear"),
    "Le_eff": (r"Le$_{eff}$ (diffusivity-weighted, unburnt)", lambda d: d["Le_eff.u.diff"], "linear"),
    "ZeLe": (r"Ze (Le$_{eff}$ − 1)", lambda d: d.Ze * (d["Le_eff.u.diff"] - 1), "linear"),
    "I0": (r"$I_0$ predicted (Schneider law, H$_2$)", lambda d: d.I0_schneider, "linear"),
    "tau_F": (r"$\tau_F$ [ms]", lambda d: d.tau_F * 1e3, "log"),
}


def small_multiples(df, key, out):
    label, fn, yscale = QUANTITIES[key]
    phis, Tus = sorted(df.phi.unique()), sorted(df.T_u.unique())
    fig, axes = plt.subplots(len(phis), len(Tus), figsize=(3.2 * len(Tus), 2.3 * len(phis)),
                             sharex=True, sharey="row" if yscale == "linear" else False,
                             squeeze=False)
    for i, phi in enumerate(phis):
        for j, Tu in enumerate(Tus):
            ax = axes[i, j]
            for P, c in P_COLORS.items():
                d = df[(df.phi == phi) & (df.T_u == Tu) & (df.P_atm == P)].sort_values("gamma")
                if d.empty:
                    continue
                y = fn(d)
                ok = d.sL >= SL_MIN
                ax.plot(d.gamma[ok], y[ok], "-o", color=c, lw=2, ms=4, label=f"{P:g} atm")
                ax.plot(d.gamma[~ok], y[~ok], "o", mfc="none", color=c, ms=4)
            ax.set_yscale(yscale)
            ax.grid(alpha=0.25, lw=0.5)
            if i == 0:
                ax.set_title(f"$T_u$ = {Tu:g} K", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"φ = {phi:g}", fontsize=10)
            if i == len(phis) - 1:
                ax.set_xlabel("γ")
            if key == "ZeLe":
                ax.axhline(0, color="0.5", lw=0.8)
            if key == "I0":
                ax.axhline(1, color="0.5", lw=0.8)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper right", ncol=3, frameon=False, fontsize=9)
    fig.suptitle(label + f"   (open: $s_L$ < {SL_MIN*100:g} cm/s)", x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=130)
    plt.close(fig)


def le_definitions(df, out, P=10.0, Tu=750.0):
    """Le_eff with the three fuel weightings, unburnt vs burnt gas."""
    phis = sorted(df.phi.unique())
    fig, axes = plt.subplots(1, len(phis), figsize=(3 * len(phis), 3), sharey=True)
    styles = {"diff": ("#0b2f6b", "diffusivity (Lehmann)"), "molar": ("#3a7dc0", "molar (Chabot)"),
              "hrr": ("#86b6e0", "heat release")}
    for ax, phi in zip(axes, phis):
        d = df[(df.phi == phi) & (df.T_u == Tu) & (df.P_atm == P)].sort_values("gamma")
        for w, (c, lab) in styles.items():
            ax.plot(d.gamma, d[f"Le_eff.u.{w}"], "-o", color=c, lw=2, ms=4, label=lab + ", unburnt")
            ax.plot(d.gamma, d[f"Le_eff.b.{w}"], "--", color=c, lw=1.5, label=lab + ", burnt")
        ax.set_title(f"φ = {phi:g}", fontsize=10)
        ax.set_xlabel("γ")
        ax.axhline(1, color="0.5", lw=0.8)
        ax.grid(alpha=0.25, lw=0.5)
    axes[0].set_ylabel(r"Le$_{eff}$")
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.suptitle(f"Le$_{{eff}}$ definitions, P = {P:g} atm, $T_u$ = {Tu:g} K", x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0.12, 1, 0.95))
    fig.savefig(out, dpi=130)
    plt.close(fig)


def soret_effect(df, dfn, out):
    m = df.merge(dfn, on=["gamma", "phi", "P_atm", "T_u"], suffixes=("", "_ns"))
    m = m[(m.status == "ok") & (m.status_ns == "ok") & (m.sL >= SL_MIN)]
    phis, Tus = sorted(m.phi.unique()), sorted(m.T_u.unique())
    fig, axes = plt.subplots(len(phis), len(Tus), figsize=(3.2 * len(Tus), 2.3 * len(phis)),
                             sharex=True, sharey=True, squeeze=False)
    for i, phi in enumerate(phis):
        for j, Tu in enumerate(Tus):
            ax = axes[i, j]
            for P, c in P_COLORS.items():
                d = m[(m.phi == phi) & (m.T_u == Tu) & (m.P_atm == P)].sort_values("gamma")
                ax.plot(d.gamma, 100 * (d.sL / d.sL_ns - 1), "-o", color=c, lw=2, ms=4, label=f"{P:g} atm")
            ax.axhline(0, color="0.5", lw=0.8)
            ax.grid(alpha=0.25, lw=0.5)
            if i == 0:
                ax.set_title(f"$T_u$ = {Tu:g} K", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"φ = {phi:g}", fontsize=10)
            if i == len(phis) - 1:
                ax.set_xlabel("γ")
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper right", ncol=3, frameon=False, fontsize=9)
    fig.suptitle("Soret effect on $s_L$ [%]  (Soret / no Soret − 1)", x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="mei2021")
    ap.add_argument("--compare", default=None)
    a = ap.parse_args()
    df = pd.read_csv(os.path.join("results", a.tag, "a0_summary.csv"))
    df = df[df.status == "ok"]
    figdir = os.path.join("results", a.tag, "figs")
    os.makedirs(figdir, exist_ok=True)
    for k in QUANTITIES:
        small_multiples(df, k, os.path.join(figdir, f"{k}.png"))
    le_definitions(df, os.path.join(figdir, "Le_eff_definitions.png"))
    if a.compare:
        dfn = pd.read_csv(os.path.join("results", a.compare, "a0_summary.csv"))
        soret_effect(df, dfn, os.path.join(figdir, "soret_effect_sL.png"))
    print("figures in", figdir)
