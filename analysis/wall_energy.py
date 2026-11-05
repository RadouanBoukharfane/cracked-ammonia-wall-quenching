"""Review point 6: cumulative wall heat per unit area E_w(t) = int <Phi_w>_x dt after the first contact, unstable vs
planar HOQ, and its sensitivity to the end of the integration window.

  python wall_energy.py  -> wall_energy.csv, figs/wall_energy.pdf (diagnostic)
E_w normalised by rho_u c_p,u s_L (T_ad - T_u) tau_F = q_ref tau_F (energy of one flame thickness of fresh gas heated
to T_ad, per unit wall area).
"""
import glob, json, os
import numpy as np, pandas as pd
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
R = os.path.join(ROOT, "PhaseB", "runs")
WIN = (5.0, 8.0, 10.0, 15.0)


def ew(run):
    m = json.load(open(os.path.join(run, "case.json"))); z = np.load(os.path.join(run, "hoq_wallflux.npz"))
    t = z["t"] / m["tau_F"]; o = np.argsort(t); t = t[o]
    phi = z["phi_w"][o].mean(1) / z["q_ref"]
    k = t > (5.0 if "planar" not in run else 0.5 * m.get("dist_lF", 30.0))
    tc = t[k][np.argmax(z["phi_w"][o][k].max(1) / z["q_ref"] > 0.05)]
    E = np.concatenate([[0], np.cumsum(0.5 * (phi[1:] + phi[:-1]) * np.diff(t))])
    E0 = np.interp(tc, t, E)
    return {f"E_{w:g}": (np.interp(tc + w, t, E) - E0) if tc + w <= t[-1] else np.nan for w in WIN} | {"t_end": t[-1] - tc}


rows = []
for g in (0.2, 0.43, 0.75, 1.0):
    for tw in (750, 500):
        sfx = "" if tw == 750 else "_Tw500"
        pl = os.path.join(R, f"HOQ_g{g:.2f}_phi0.50_P10.0_T750_planar{sfx}")
        un = sorted(r for r in glob.glob(os.path.join(R, f"HOQ_g{g:.2f}_phi0.50_P10.0_T750_t??.0{sfx}"))
                    if os.path.exists(os.path.join(r, "hoq_wallflux.npz")))
        if not un or not os.path.exists(os.path.join(pl, "hoq_wallflux.npz")):
            continue
        P = ew(pl)
        for r in un:
            U = ew(r)
            rows.append(dict(run=os.path.basename(r), gamma=g, Tw=tw, t_end=U["t_end"],
                             **{f"ratio_{w:g}": U[f"E_{w:g}"] / P[f"E_{w:g}"] for w in WIN},
                             **{f"Eu_{w:g}": U[f"E_{w:g}"] for w in WIN}, **{f"Ep_{w:g}": P[f"E_{w:g}"] for w in WIN}))
d = pd.DataFrame(rows); d.to_csv(os.path.join(ROOT, "analysis", "wall_energy.csv"), index=False)
pd.set_option("display.width", 200)
print(d.groupby(["gamma", "Tw"])[[f"ratio_{w:g}" for w in WIN] + ["t_end"]].agg(["mean", "std"]).round(3).to_string())
print(d.groupby(["gamma", "Tw"])[[f"Ep_{w:g}" for w in WIN]].mean().round(3).to_string())
