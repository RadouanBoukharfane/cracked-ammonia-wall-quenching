"""Wall-temperature series (planar HOQ): table + figure. Usage: python summary_tw.py [DT_AFTER_TQ=10]"""
import glob, json, os, re, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

DT = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
rows = []
for r in sorted(glob.glob("runs/HOQ_g*_phi0.50_P10.0_T750_planar*")):
    f = os.path.join(r, "hoq_timeseries.csv")
    if not os.path.exists(f): continue
    m = json.load(open(os.path.join(r, "case.json")))
    d = pd.read_csv(f)
    late = d[d.t_tauF > 0.5 * m.get("dist_lF", 30.0)]
    iq = late.Phi_max_star.idxmax(); tq = d.t_tauF[iq]
    if tq + DT > d.t_tauF.max() + 0.3: print("too short:", r); continue
    a = d.loc[(d.t_tauF - (tq + DT)).abs().idxmin()]
    q = d.loc[iq]
    rows.append(dict(gamma=m["case"]["gamma"], Tw=m["T_wall"], dist=m.get("dist_lF", 30.0), tq=tq,
                     Phi_q=q.Phi_max_star, Pe_q=q.Pe_min, NH3_q=q.X_NH3_y2 * 100,
                     NH3=a.X_NH3_y2 * 100, NO=a.X_NO_y2 * 1e6, N2O=a.X_N2O_y2 * 1e6, N2O_NO=a.X_N2O_y2 / a.X_NO_y2))
t = pd.DataFrame(rows).sort_values(["gamma", "dist", "Tw"])
t.to_csv("tw_summary.csv", index=False)
print(f"values at t_q + {DT:g} tau_F, y < 2 l_F")
print(t.to_string(index=False, float_format="%.3g"))

COL = {0.2: "#2a78d6", 0.43: "#eb6834", 0.75: "#1baf7a", 1.0: "#eda100"}
MK = {0.2: "o", 0.43: "s", 0.75: "^", 1.0: "D"}
LS = {15.0: ":", 30.0: "-", 60.0: "--"}
panels = [("Phi_q", r"$\Phi_q^*$ (flux pariétal max.)"), ("Pe_q", r"Pe au flux max."),
          ("NH3", rf"NH$_3$ résiduel à $t_q$+{DT:g} τ$_F$ [%]"), ("NO", rf"NO à $t_q$+{DT:g} τ$_F$ [ppm]"),
          ("N2O", rf"N$_2$O à $t_q$+{DT:g} τ$_F$ [ppm]"), ("N2O_NO", rf"N$_2$O/NO à $t_q$+{DT:g} τ$_F$")]
fig, axs = plt.subplots(2, 3, figsize=(13, 7.5))
for ax, (k, lab) in zip(axs.flat, panels):
    for (g, dist), s in t.groupby(["gamma", "dist"]):
        if g == 1.0 and k in ("NH3", "NO", "N2O", "N2O_NO"): continue   # no NH3 in the fuel
        ax.plot(s.Tw, s[k], ls=LS[dist], lw=2, color=COL[g], marker=MK[g], ms=8,
                mec="#fcfcfb", mew=1.5, label=f"γ = {g:.2f}" + ("" if dist == 30 else f", d = {dist:.0f} l$_F$"))
    if k in ("NH3",): ax.set_yscale("log")
    ax.set_title(lab, fontsize=10, color="#0b0b0b", loc="left")
    ax.set_xlabel(r"$T_w$ [K]  ($T_u$ = 750 K)", color="#52514e")
    ax.grid(color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): ax.spines[sp].set_color("#b5b3ab")
    ax.tick_params(colors="#52514e")
axs[0, 0].legend(fontsize=8, frameon=False, handlelength=3.5)
axs[1, 2].text(0.98, 0.98, "γ = 1 (sans NH$_3$) omis\ndes panneaux d'émission", transform=axs[1, 2].transAxes,
               ha="right", va="top", fontsize=8, color="#52514e")
fig.suptitle("HOQ plan, φ = 0,5, 10 atm : effet de la température de paroi", x=0.01, ha="left", fontsize=12)
fig.tight_layout(); fig.savefig("figs/tw_series.png", dpi=130, facecolor="#fcfcfb")
