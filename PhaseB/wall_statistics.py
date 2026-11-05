"""Local quenching statistics of all HOQ runs (TODO items 2 and 3), from RUN/hoq_wallflux.npz (post_hoq.py with
per-column y_hrr / hrr_col).

  python wall_statistics.py      -> wall_events.csv (one row per run and wall column), wall_stats.csv (one row per run)

Per wall column x:
  Phi_q(x)  = max_t Phi_w*(x, t), start-up excluded (t > dist/2 tau_F for planar runs, t > 5 tau_F for unstable runs);
  Pe_q(x)   = min_t y_HRR(x, t)/l_F over the times at which the column still burns (peak HRR of the column above
              10 % of the free-flame peak), i.e. the closest approach of the reaction zone; removes the residual
              wall heat release after extinction that biased the earlier "Pe_min";
  Pe_at_q(x)= y_HRR/l_F at the instant of Phi_q(x) (when the column burns).
Per run (item 3): peak and mean wall heat flux in MW/m2, time-averaged and RMS flux over the quenching phase
(first contact to +15 tau_F), fraction of the wall whose local peak exceeds the planar peak.
"""
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
WINDOW = 15.0           # tau_F after the first contact, for time-averaged / RMS flux


def analyse(run):
    m = json.load(open(os.path.join(run, "case.json")))
    z = np.load(os.path.join(run, "hoq_wallflux.npz"))
    if "y_hrr" not in z:
        return None
    t = z["t"] / m["tau_F"]
    o = np.argsort(t)
    t, phi, yh, hc = t[o], z["phi_w"][o] / z["q_ref"], z["y_hrr"][o] / m["l_F"], z["hrr_col"][o]
    unstable = "planar" not in os.path.basename(run)
    t_start = 5.0 if unstable else 0.5 * m.get("dist_lF", 30.0)
    k = t > t_start
    hrr_ref = np.median(hc[:5].max(axis=1))
    burning = hc > 0.1 * hrr_ref
    phi_k, yh_k, b_k, t_k = phi[k], yh[k], burning[k], t[k]
    iq = np.argmax(phi_k, axis=0)
    nx = phi.shape[1]
    Phi_q = phi_k[iq, np.arange(nx)]
    Pe_q = np.where(b_k, yh_k, np.inf).min(axis=0)
    Pe_at = np.where(b_k[iq, np.arange(nx)], yh_k[iq, np.arange(nx)], np.nan)
    ev = pd.DataFrame(dict(x_lF=(np.arange(nx) + 0.5) * m["Lx"] / nx / m["l_F"], Phi_q=Phi_q, Pe_q=Pe_q,
                           Pe_at_q=Pe_at, t_q=t_k[iq]))
    ev.loc[~np.isfinite(ev.Pe_q), "Pe_q"] = np.nan
    gmax = phi_k.max(axis=1)
    t_contact = t_k[np.argmax(gmax > 0.05)]
    w = (t >= t_contact) & (t <= t_contact + WINDOW)
    qref = float(z["q_ref"]) / 1e6
    st = dict(run=os.path.basename(run), gamma=round(m["case"]["gamma"], 2), Tw=m["T_wall"],
              dist=m.get("dist_lF", np.nan), unstable=int(unstable), crest="crest" in run,
              q_ref_MW=qref, Phi_peak=Phi_q.max(), Phi_peak_MW=Phi_q.max() * qref,
              Phi_mean_window=phi[w].mean(), Phi_rms_window=phi[w].std(), Phi_local_mean=Phi_q.mean(),
              Phi_local_p95=np.percentile(Phi_q, 95), Pe_q_min=np.nanmin(ev.Pe_q), Pe_q_median=np.nanmedian(ev.Pe_q),
              t_contact=t_contact)
    return ev, st


def main():
    runs = sorted(r for r in glob.glob(os.path.join(HERE, "runs", "HOQ_*"))
                  if "_mech_" not in r and not r.endswith("_L3"))   # mechanism and grid checks: not in the statistics
    E, S = [], []
    for r in runs:
        if not os.path.exists(os.path.join(r, "hoq_wallflux.npz")):
            continue
        res = analyse(r)
        if res is None:
            continue
        ev, st = res
        ev["run"] = st["run"]
        for k in ("gamma", "Tw", "unstable", "crest"):
            ev[k] = st[k]
        E.append(ev); S.append(st)
    E, S = pd.concat(E), pd.DataFrame(S)
    # fraction of the wall above the planar peak (same gamma and T_w, d0 = 30)
    frac = []
    for _, s in S.iterrows():
        ref = S[(S.gamma == s.gamma) & (S.Tw == s.Tw) & (S.unstable == 0) & (~S.crest) & S.run.str.contains("_T750_") &
                (S.run.str.endswith("planar") | S.run.str.endswith(f"planar_Tw{s.Tw:.0f}"))]
        if s.unstable and len(ref):
            ev = E[E.run == s.run]
            frac.append((ev.Phi_q > ref.Phi_peak.values[0]).mean())
        else:
            frac.append(np.nan)
    S["frac_above_planar"] = frac
    E.to_csv(os.path.join(HERE, "wall_events.csv"), index=False)
    S.to_csv(os.path.join(HERE, "wall_stats.csv"), index=False)
    pd.set_option("display.width", 250)
    cols = ["run", "Phi_peak", "Phi_peak_MW", "Phi_local_mean", "Phi_mean_window", "Phi_rms_window", "Pe_q_min",
            "Pe_q_median", "frac_above_planar"]
    print(S[(S.unstable == 1) | S.run.str.endswith("planar")][cols].to_string(index=False, float_format="%.3g"))


if __name__ == "__main__":
    main()
