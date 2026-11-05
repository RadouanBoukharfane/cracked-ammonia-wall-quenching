"""Review point 7: sensitivity of the cold-layer collapse (H) to the evaluation instant (front Y_REF l_F from the wall)
and to the preheat cut-off (PREHEAT l_F before the front). R^2 of one quadratic fit through d0 = 15, 30, 60 l_F
(gamma = 0.43), as in collapse_tw.py.   Run from PhaseB/:  python ../analysis/h_sensitivity.py
"""
import sys
import numpy as np, pandas as pd
sys.argv = [sys.argv[0]]
sys.path.insert(0, ".")
import collapse_tw as C
s = pd.read_csv("tw_summary.csv"); s = s[s.gamma == 0.43]
out = []
for yref in (6.0, 8.0, 10.0, 12.0):
    for pre in (1.0, 2.0, 3.0):
        C.Y_REF, C.PREHEAT = yref, pre
        rows = []
        for _, r in s.iterrows():
            run = f"runs/HOQ_g{r.gamma:.2f}_phi0.50_P10.0_T750_planar"
            if not (r.Tw == 750 and r.dist == 30):
                run += f"_Tw{r.Tw:.0f}" + ("" if r.dist == 30 else f"_d{r.dist:.0f}")
            rows.append(dict(Phi_q=r.Phi_q, **C.analyse(run)))
        d = pd.DataFrame(rows); d["NH3_frac"] = d.NH3 / C.NH3_IN[0.43]
        d["NO_rel"] = d.NO / d.NO[d.theta_w == 0].mean()
        row = dict(Y_REF=yref, PREHEAT=pre)
        for k in ("Phi_q", "Pe_Q", "NH3_frac", "NO_rel", "N2O_NO"):
            for c in ("H", "theta_w"):
                p = np.polyfit(d[c], d[k], 2)
                row[f"{c}:{k}"] = 1 - ((d[k] - np.polyval(p, d[c]))**2).sum() / ((d[k] - d[k].mean())**2).sum()
        out.append(row)
o = pd.DataFrame(out)
o.to_csv("../analysis/h_sensitivity.csv", index=False)
pd.set_option("display.width", 250)
print(o[["Y_REF", "PREHEAT"] + [c for c in o.columns if c.startswith("H:")]].round(3).to_string(index=False))
print("theta_w (independent of Y_REF/PREHEAT):", o[[c for c in o.columns if c.startswith("theta_w:")]].iloc[0].round(3).to_dict())
