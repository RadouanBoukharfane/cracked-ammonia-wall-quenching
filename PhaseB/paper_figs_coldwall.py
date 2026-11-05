"""Cold-wall figures of the paper (printed width 6.46 in, 9 pt): ../paper/figs/coldwall_collapse.pdf, coldwall_NO_budget.pdf
Needs tw_collapse.csv, tw_collapse_R2.csv (collapse_tw.py) and budget_reactions_NO_y2.csv, budget_y2.csv (no_budget.py)."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.lines import Line2D

W = 6.46
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"
plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "mathtext.fontset": "cm", "font.family": "serif", "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.axisbelow": True, "axes.edgecolor": "#b5b3ab", "xtick.color": INK2, "ytick.color": INK2,
    "savefig.dpi": 300})
OUT = "../paper/figs/"
COL = {0.2: "#2a78d6", 0.43: "#eb6834", 0.75: "#1baf7a"}
MK = {0.2: "o", 0.43: "s", 0.75: "^"}
LS = {15.0: (0, (1, 1.5)), 30.0: "-", 60.0: (0, (5, 2))}
HALO = [pe.withStroke(linewidth=2.5, foreground="white")]


def tag(ax, i, x=-0.02):
    ax.text(x, 1.03, f"({chr(97 + i)})", transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right",
            va="bottom")


# ---------------------------------------------------------------------------------------------------------------
# Message: at fixed T_w the initial distance matters; the cold-layer enthalpy deficit H collapses it.
d = pd.read_csv("tw_collapse.csv")
R = pd.read_csv("tw_collapse_R2.csv", index_col=0)
tg = [("Pe_Q", r"$\mathrm{Pe}_Q$", "quenching distance"),
      ("NH3_frac", r"$X_{\mathrm{NH_3}}/X_{\mathrm{NH_3},u}$", "residual ammonia"),
      ("N2O_NO", r"$X_{\mathrm{N_2O}}/X_{\mathrm{NO}}$", "N$_2$O-to-NO ratio")]
fig = plt.figure(figsize=(W, 5.0))
gs = fig.add_gridspec(2, 3, hspace=0.38, wspace=0.36, left=0.145, right=0.985, top=0.86, bottom=0.09)
k_ = 0
for row, (xc, xl, rk) in enumerate([("Tw", r"$T_w$ (K)", "theta_w"), ("H", r"$\mathcal{H}$", "H")]):
    for col, (k, lab, title) in enumerate(tg):
        ax = fig.add_subplot(gs[row, col])
        if row == 0:
            q = d[d.gamma == 0.43].pivot(index="Tw", columns="dist", values=k).sort_index()
            ax.fill_between(q.index, q.min(axis=1), q.max(axis=1), color=COL[0.43], alpha=0.13, lw=0)
            ax.set_title(title, fontsize=9, color=INK, pad=8)
        for (g, dist), q in d.groupby(["gamma", "dist"]):
            q = q.sort_values(xc)
            ax.plot(q[xc], q[k], ls=LS[dist], lw=1.4, color=COL[g], marker=MK[g], ms=4.5, mec="white", mew=0.6)
        ypos = 0.05
        if k == "NH3_frac":
            ax.set_yscale("log")
        if k == "N2O_NO":
            ax.axhspan(1, 4, color="#f4f3f0", zorder=0)
            ax.axhline(1, color=INK2, lw=0.7, ls=(0, (2, 2)))
            ax.text(0.97 if row == 0 else 0.03, 0.95, "N$_2$O dominant", transform=ax.transAxes, va="top",
                    ha="right" if row == 0 else "left", fontsize=7.5, color=INK2)
            ax.set_ylim(0, 3.8)
            ypos = 0.62 if row == 0 else 0.78
        ax.text(0.97, ypos, rf"$R^2 = {R.loc[rk, k]:.2f}$", transform=ax.transAxes, ha="right", fontsize=8,
                color=INK, path_effects=HALO)
        ax.set_xlabel(xl); ax.set_ylabel(lab)
        tag(ax, k_, x=-0.17); k_ += 1
fig.text(0.008, 0.68, "against wall temperature", rotation=90, ha="left", va="center", fontsize=9, color=INK2)
fig.text(0.008, 0.29, "against cold-layer\nenthalpy deficit", rotation=90, ha="left", va="center", fontsize=9,
         color=INK2, multialignment="center")
from matplotlib.patches import Patch
h1 = [Line2D([], [], color=COL[g], marker=MK[g], mec="white", ms=6, lw=1.4, label=rf"$\gamma = {g:g}$")
      for g in (0.2, 0.43, 0.75)]
h2 = [Line2D([], [], color=INK2, ls=LS[x], lw=1.4, label=rf"$d_0 = {x:.0f}\,l_F$") for x in (15.0, 30.0, 60.0)]
h3 = [Patch(color=COL[0.43], alpha=0.2, label=r"$d_0$ spread")]
fig.legend(handles=h1 + h2 + h3, loc="upper center", ncol=7, frameon=False, bbox_to_anchor=(0.55, 1.0),
           handlelength=2.2, columnspacing=0.9, handletextpad=0.4, fontsize=7.5)
fig.savefig(OUT + "coldwall_collapse.pdf")
plt.close(fig)

# ---------------------------------------------------------------------------------------------------------------
# Message: NO falls near a cold wall because its radical-driven production collapses, not through a new sink.
P = pd.read_csv("budget_reactions_NO_y2.csv", index_col=0, header=[0, 1])
P.columns = pd.MultiIndex.from_tuples([(float(a), int(b)) for a, b in P.columns])
B = pd.read_csv("budget_y2.csv")
hot, cold = (0.43, 750), (0.43, 400)
top = P[[hot, cold]].abs().max(axis=1).sort_values(ascending=False).index[:9][::-1]


def chem(r):
    for a, b in [("<=>", r"$\rightarrow$"), ("NH2", "NH$_2$"), ("NO2", "NO$_2$"), ("HO2", "HO$_2$"), ("H2O", "H$_2$O"),
                 ("N2O", "N$_2$O"), ("O2", "O$_2$"), ("H2NO", "H$_2$NO"), ("N2", "N$_2$"), ("H2 ", "H$_2$ ")]:
        r = r.replace(a, b)
    return r.replace("H$_2$NO", "H$_2$NO")


fig = plt.figure(figsize=(W, 3.4))
gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.25], wspace=1.05, left=0.1, right=0.985, top=0.9, bottom=0.15)
ax, bx = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
s = B[(B.sp == "NO") & (B.gamma == 0.43)].sort_values("Tw", ascending=False)
x = np.arange(len(s))
for j, (k, c, lb) in enumerate([("form", "#2a78d6", "production"), ("cons", "#eb6834", "consumption"),
                                ("transport", "#1baf7a", "transport")]):
    ax.bar(x + (j - 1) * 0.27, s[k], width=0.26, color=c, label=lb, edgecolor="white", lw=0.5)
ax.set_xticks(x, [f"{v:.0f}" for v in s.Tw]); ax.set_xlabel(r"wall temperature $T_w$ (K)")
ax.axhline(0, color=INK2, lw=0.7)
ax.set_ylabel(r"NO budget, $0 < y < 2\,l_F$ ($10^{-9}$ kmol m$^{-2}$)")
ax.annotate(r"production $\div 16$", xy=(4 - 0.27, s.form.iloc[-1]), xytext=(1.6, 520), fontsize=8, color="#1c5cab",
            arrowprops=dict(arrowstyle="-|>", color="#1c5cab", lw=0.8, connectionstyle="arc3,rad=-0.25"))
ax.legend(frameon=False, loc="lower right", fontsize=7.5, handlelength=1.0, bbox_to_anchor=(1.05, 0.0))
ax.set_ylim(-800, 900)
tag(ax, 0, x=-0.24)
yy = np.arange(len(top))
bx.barh(yy + 0.2, P.loc[top, hot], height=0.38, color="#c3c2b7", label=r"$T_w = 750$ K ($= T_u$)")
bx.barh(yy - 0.2, P.loc[top, cold], height=0.38, color="#0d366b", label=r"$T_w = 400$ K")
bx.set_yticks(yy, [chem(r) for r in top], fontsize=7.8)
bx.axvline(0, color=INK2, lw=0.7)
bx.set_xlabel(r"NO per reaction ($10^{-9}$ kmol m$^{-2}$)")
bx.set_ylim(-0.7, len(top) + 0.4)
bx.legend(frameon=False, loc="lower left", fontsize=7.5)
bx.grid(axis="y", visible=False)
bx.text(0.98, 0.98, "formation", transform=bx.transAxes, ha="right", va="top", fontsize=8, color="#1c5cab")
bx.text(0.02, 0.98, "consumption", transform=bx.transAxes, ha="left", va="top", fontsize=8, color="#b2262a")
tag(bx, 1, x=-0.02)
fig.savefig(OUT + "coldwall_NO_budget.pdf")
plt.close(fig)
