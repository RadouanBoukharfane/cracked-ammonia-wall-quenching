"""Cold-wall series: search a parameter that collapses the three initial distances (gamma = 0.43).

  python collapse_tw.py [DT=10]   (needs tw_summary.csv from summary_tw.py and RUN/coldwall.npz)

With T_w < T_u the flame stalls 1.5-2 l_F from the wall several tau_F BEFORE the wall heat flux peaks (the peak is
carried by conduction from the burnt gas), so the timing is based on heat release, not on the flux:
  t_Q  = first time the integrated heat release int HRR dy drops below 50 % of its free-flame value
  y_Q  = position of the HRR peak at t_Q (Pe_Q = y_Q / l_F)
  emissions: x-averaged mole fractions over 0 < y < 2 l_F at t_Q + DT tau_F
Cold-layer metrics on the profile when the flame leading edge (T_track) is Y_REF l_F from the wall, fresh gas
y < y_f - PREHEAT l_F:
  theta_w = (T_u - T_w)/(T_ad - T_u)                                 (wall only)
  delta_c = int (T_u - T)/(T_u - T_w) dy / l_F                        (cold-layer thickness)
  H       = int (rho/rho_u)(T_u - T) dy / ((T_ad - T_u) l_F)          (enthalpy deficit / flame enthalpy per l_F)
  T_eff   = mean T over 0 < y < y_Q                                   (gas temperature where the flame will quench)
  Theta_Q = (T_u - T_eff)/(T_ad - T_u)
Collapse quality: R^2 of one quadratic fit through the three distances (gamma = 0.43, 15 points).
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DT = float(sys.argv[1]) if len(sys.argv) > 1 else 10.0
Y_REF, PREHEAT = 8.0, 2.0
NH3_IN = {0.2: 9.59, 0.43: 6.65, 0.75: 2.81}


def analyse(run):
    m = json.load(open(os.path.join(run, "case.json")))
    z = np.load(os.path.join(run, "coldwall.npz"))
    y, T, rho, hrr, X, t = z["y_lF"], z["T"], z["rho"], z["hrr"], z["X"], z["t_tauF"]
    sp = list(z["species"])
    dy = y[1] - y[0]
    Ihrr = hrr.sum(1) * dy
    I0 = np.median(Ihrr[:5])
    kq = np.argmax(Ihrr < 0.5 * I0)
    f = (Ihrr[kq - 1] - 0.5 * I0) / (Ihrr[kq - 1] - Ihrr[kq])
    tQ = t[kq - 1] + f * (t[kq] - t[kq - 1])
    yQ = y[np.argmax(hrr[kq - 1])] * (1 - f) + y[np.argmax(hrr[kq])] * f
    yf = np.array([y[np.argmax(Tk > m["T_track"])] for Tk in T])
    pre = np.where(t < tQ)[0]
    k = pre[np.argmin(np.abs(yf[pre] - Y_REF))]
    Tu, Tw, Tad = m["Tu"], m["T_wall"], m["Tb_eq"]
    sel = y < yf[k] - PREHEAT
    deficit = np.clip(Tu - T[k, sel], 0, None)
    rho_u = rho[k, sel][-1] * T[k, sel][-1] / Tu
    T_eff = T[k, y < yQ].mean()
    ke = np.argmin(np.abs(t - (tQ + DT)))
    near = y < 2.0
    Xn = X[ke][near].mean(0)
    i = {s: sp.index(s) for s in ("NH3", "NO", "N2O", "H2")}
    return dict(tQ=tQ, Pe_Q=yQ, yf_ref=yf[k], t_ref=t[k] - tQ, T_edge=T[k, sel][-1], t_emis=t[ke] - tQ,
                theta_w=(Tu - Tw) / (Tad - Tu), delta_c=deficit.sum() * dy / (Tu - Tw) if Tu > Tw else 0.0,
                H=(rho[k, sel] / rho_u * deficit).sum() * dy / (Tad - Tu), T_eff=T_eff,
                Theta_Q=(Tu - T_eff) / (Tad - Tu),
                NH3=Xn[i["NH3"]] * 100, NO=Xn[i["NO"]] * 1e6, N2O=Xn[i["N2O"]] * 1e6, N2O_NO=Xn[i["N2O"]] / Xn[i["NO"]])


def main():
    s = pd.read_csv("tw_summary.csv")
    s = s[s.gamma < 1]
    rows = []
    for _, r in s.iterrows():
        run = f"runs/HOQ_g{r.gamma:.2f}_phi0.50_P10.0_T750_planar"
        if not (r.Tw == 750 and r.dist == 30):
            run += f"_Tw{r.Tw:.0f}" + ("" if r.dist == 30 else f"_d{r.dist:.0f}")
        rows.append(dict(gamma=r.gamma, dist=r.dist, Tw=r.Tw, Phi_q=r.Phi_q, t_flux=r.tq, **analyse(run)))
    d = pd.DataFrame(rows)
    d["t_flux"] -= d.tQ                                       # delay of the flux peak after the HRR quench
    d["NH3_frac"] = d.NH3 / d.gamma.map(NH3_IN)
    d["NO_rel"] = d.NO / d.gamma.map(d[d.Tw == 750].groupby("gamma").NO.mean())
    d.to_csv("tw_collapse.csv", index=False)
    pd.set_option("display.width", 220)
    print(f"emissions at t_Q + {DT:g} tau_F, y < 2 l_F ; t_flux = delay of the wall-flux peak after t_Q")
    print(d[["gamma", "dist", "Tw", "tQ", "t_flux", "Pe_Q", "t_ref", "T_edge", "delta_c", "H", "T_eff", "Theta_Q",
             "Phi_q", "NH3_frac", "NO", "N2O", "N2O_NO"]].to_string(index=False, float_format="%.3g"))

    g = d[d.gamma == 0.43]
    targets = ["Phi_q", "Pe_Q", "NH3_frac", "NO_rel", "N2O_NO"]
    cands = ["theta_w", "delta_c", "H", "Theta_Q"]
    R = pd.DataFrame(index=cands, columns=targets, dtype=float)
    for c in cands:
        for k in targets:
            p = np.polyfit(g[c], g[k], 2)
            R.loc[c, k] = 1 - ((g[k] - np.polyval(p, g[c]))**2).sum() / ((g[k] - g[k].mean())**2).sum()
    print("\nR^2, one quadratic fit through d = 15, 30, 60 l_F (gamma = 0.43)")
    print(R.round(3).to_string())
    R.to_csv("tw_collapse_R2.csv")
    best = R.mean(axis=1).idxmax()
    print("best (mean R^2):", best)

    COL = {0.2: "#2a78d6", 0.43: "#eb6834", 0.75: "#1baf7a"}
    MK = {0.2: "o", 0.43: "s", 0.75: "^"}
    LS = {15.0: ":", 30.0: "-", 60.0: "--"}
    lab = {"Phi_q": r"$\Phi_q^*$", "Pe_Q": r"Pe$_Q$", "NH3_frac": r"NH$_3$ résiduel / entrant",
           "NO_rel": r"NO / NO$(T_w = T_u)$", "N2O_NO": r"N$_2$O/NO"}
    xl = {"theta_w": r"$\theta_w = (T_u - T_w)/(T_{ad} - T_u)$", "delta_c": r"$\delta_c / l_F$",
          "H": r"$H$", "Theta_Q": r"$\Theta_Q = (T_u - T_{eff})/(T_{ad} - T_u)$"}
    fig, axs = plt.subplots(2, len(targets), figsize=(17, 6.8))
    for row, xc in enumerate(["theta_w", best if best != "theta_w" else "Theta_Q"]):
        for ax, k in zip(axs[row], targets):
            for (gm, dist), q in d.groupby(["gamma", "dist"]):
                q = q.sort_values(xc)
                ax.plot(q[xc], q[k], ls=LS[dist], lw=2, color=COL[gm], marker=MK[gm], ms=7, mec="#fcfcfb", mew=1.2,
                        label=f"γ = {gm:.2f}" + ("" if dist == 30 else f", d = {dist:.0f} l$_F$"))
            if k == "NH3_frac":
                ax.set_yscale("log")
            ax.set_xlabel(xl[xc], fontsize=9, color="#52514e")
            ax.set_title(f"{lab[k]}   R² (γ = 0,43) = {R.loc[xc, k]:.2f}", fontsize=9, loc="left")
            ax.grid(color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
            for sp_ in ("top", "right"):
                ax.spines[sp_].set_visible(False)
    axs[0, 0].legend(fontsize=7, frameon=False, handlelength=3.5)
    fig.tight_layout()
    fig.savefig("figs/tw_collapse.png", dpi=130, facecolor="#fcfcfb")


if __name__ == "__main__":
    main()
