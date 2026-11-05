"""Species budget (NO, N2O, NO2, NH3) in the near-wall layer 0 < y < LAYER l_F of the planar HOQ runs, per reaction.

  python no_budget.py [LAYER=2]   (needs RUN/coldwall.npz from coldwall_analysis.py)

Time integral from t_Q - 2 to t_Q + 10 tau_F (t_Q: integrated heat release < 50 % of the free flame, as in
collapse_tw.py) of the net production of every reaction (duplicate reactions summed), in 1e-9 kmol/m2;
formation = sum of the reactions with a positive integral, consumption = sum of the negative ones;
transport across y = LAYER l_F = change of the inventory - chemistry.
Output: budget_y{LAYER}.csv, budget_reactions_{SP}_y{LAYER}.csv, figs/no_budget_y{LAYER}.png.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LAYER = sys.argv[1] if len(sys.argv) > 1 else "2"
T0, T1 = -2.0, 10.0
SPECIES = ("NO", "N2O", "NO2", "NH3")
CASES = [(0.43, 750), (0.43, 650), (0.43, 600), (0.43, 500), (0.43, 400), (0.20, 750), (0.20, 500), (0.20, 400),
         (0.75, 750), (0.75, 500), (0.75, 400)]


def run_name(g, tw):
    return f"runs/HOQ_g{g:.2f}_phi0.50_P10.0_T750_planar" + ("" if tw == 750 else f"_Tw{tw}")


def integ(t, f, a, b):
    k = (t >= a - 1e-9) & (t <= b + 1e-9)
    return np.trapezoid(f[k], t[k], axis=0)


rows, per = [], {s: {} for s in SPECIES}
for g, tw in CASES:
    run = run_name(g, tw)
    m = json.load(open(os.path.join(run, "case.json")))
    z = np.load(os.path.join(run, "coldwall.npz"))
    I = z["hrr"].sum(1); I0 = np.median(I[:5]); k = np.argmax(I < 0.5 * I0)
    tQ = z["t_tauF"][k - 1] + (I[k - 1] - 0.5 * I0) / (I[k - 1] - I[k]) * (z["t_tauF"][k] - z["t_tauF"][k - 1])
    t = (z["t_tauF"] - tQ) * m["tau_F"]
    a, b = T0 * m["tau_F"], T1 * m["tau_F"]
    for s in SPECIES:
        W = pd.Series(integ(t, z[f"w{s}_y{LAYER}"], a, b) * 1e9, index=z["reactions"]).groupby(level=0).sum()
        inv = z[f"inv{s}_y{LAYER}"] * 1e9
        di = np.interp(b, t, inv) - np.interp(a, t, inv)
        rows.append(dict(sp=s, gamma=g, Tw=tw, inv_start=np.interp(a, t, inv), inv_end=np.interp(b, t, inv),
                         chem=W.sum(), form=W[W > 0].sum(), cons=W[W < 0].sum(), transport=di - W.sum()))
        per[s][(g, tw)] = W

B = pd.DataFrame(rows)
B.to_csv(f"budget_y{LAYER}.csv", index=False)
pd.set_option("display.width", 230)
print(f"budgets over 0 < y < {LAYER} l_F, t_Q{T0:+g} .. t_Q{T1:+g} tau_F  [1e-9 kmol/m2]")
for s in SPECIES:
    print(f"\n== {s}")
    print(B[B.sp == s].drop(columns="sp").set_index(["gamma", "Tw"]).round(2).to_string())
    P = pd.DataFrame(per[s])
    P.to_csv(f"budget_reactions_{s}_y{LAYER}.csv")
    top = P.abs().max(axis=1).sort_values(ascending=False).index[:10]
    print(P.loc[top].round(2).to_string())

# figure: NO per-reaction contributions and NO budget, gamma = 0.43
P = pd.DataFrame(per["NO"])
cols = [(0.43, 750), (0.43, 600), (0.43, 500), (0.43, 400)]
top = P[cols].abs().max(axis=1).sort_values(ascending=False).index[:10][::-1]
fig, axs = plt.subplots(1, 2, figsize=(14, 5.8), gridspec_kw=dict(width_ratios=[2.1, 1]))
shades = ["#9ec5f4", "#5598e7", "#256abf", "#0d366b"]
h = 0.8 / len(cols)
for i, c in enumerate(cols):
    axs[0].barh(np.arange(len(top)) + (i - (len(cols) - 1) / 2) * h, P.loc[top, c], height=h * 0.9,
                color=shades[i], label=f"$T_w$ = {c[1]} K")
axs[0].set_yticks(np.arange(len(top)), top, fontsize=8)
axs[0].axvline(0, color="#52514e", lw=0.8)
axs[0].set_xlabel(r"NO produit (> 0) / consommé (< 0) sur $y < 2\,l_F$, $t_Q - 2$ à $t_Q + 10\,\tau_F$"
                  r" [10$^{-9}$ kmol m$^{-2}$]", fontsize=9)
axs[0].legend(fontsize=8, frameon=False, loc="lower left")
axs[0].set_title("γ = 0,43 : principales réactions du NO près de la paroi", fontsize=10, loc="left")
s = B[(B.sp == "NO") & (B.gamma == 0.43)].sort_values("Tw", ascending=False)
x = np.arange(len(s))
for j, (k, c, lb) in enumerate([("form", "#2a78d6", "formation"), ("cons", "#eb6834", "consommation"),
                                ("transport", "#1baf7a", "transport (résidu)")]):
    axs[1].bar(x + (j - 1) * 0.27, s[k], width=0.25, color=c, label=lb)
axs[1].set_xticks(x, [f"{v:.0f} K" for v in s.Tw])
axs[1].axhline(0, color="#52514e", lw=0.8)
axs[1].set_xlabel("$T_w$", fontsize=9)
axs[1].set_title("γ = 0,43 : bilan du NO (y < 2 l$_F$)", fontsize=10, loc="left")
axs[1].legend(fontsize=8, frameon=False)
for ax in axs:
    ax.grid(color="#e6e5e0", lw=0.6, axis="x" if ax is axs[0] else "y"); ax.set_axisbelow(True)
    for sp_ in ("top", "right"):
        ax.spines[sp_].set_visible(False)
fig.tight_layout()
fig.savefig(f"figs/no_budget_y{LAYER}.png", dpi=130, facecolor="#fcfcfb")
