"""Near-wall emissions at fixed delays after a common, heat-release-based quenching instant (planar and unstable HOQ).

  python emissions_tq.py NPROC RUN [RUN ...]     (needs RUN/hoq_timeseries.csv from post_hoq.py)
  -> RUN/hrr_total.csv (t, integrated heat release) and emissions_tq.csv (one line per run, appended/replaced)

t_Q: first time the heat release integrated over the whole domain drops below 50 % of its value before any contact
with the wall (median of the first 5 plotfiles). Same definition for planar and unstable runs (for unstable runs the
quenching is spread over ~10 tau_F; t_Q marks the middle of the process). t_contact: first time Phi_w* > 0.05.
Emissions: x-averaged mole fractions over 0 < y < 2 l_F (post_hoq.py), interpolated at t_Q + 0, 5, 10 tau_F.
"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

DELAYS = (0.0, 5.0, 10.0)


def total_hrr(p):
    import yt
    yt.set_log_level(50)
    ds = yt.load(p)
    ad = ds.all_data()
    v = (ad["boxlib", "HeatRelease"] * ad["index", "dx"] * ad["index", "dy"]).sum()   # finest data only
    return p, float(ds.current_time), float(v)


def main(nproc, runs):
    from multiprocessing import Pool
    jobs = []
    for run in runs:
        if os.path.exists(os.path.join(run, "hrr_total.csv")):
            continue
        jobs += [p for p in glob.glob(os.path.join(run, "plt*")) if re.fullmatch(r"plt\d+", os.path.basename(p))]
    with Pool(nproc) as pool:
        res = pool.map(total_hrr, jobs)
    by = {}
    for p, t, v in res:
        by.setdefault(os.path.dirname(p), []).append((t, v))
    rows = []
    for run in runs:
        m = json.load(open(os.path.join(run, "case.json")))
        f = os.path.join(run, "hrr_total.csv")
        if run in by:
            pd.DataFrame(sorted(by[run]), columns=["t", "hrr_tot"]).to_csv(f, index=False)
        h = pd.read_csv(f)
        h["t_tauF"] = h.t / m["tau_F"]
        h0 = np.median(h.hrr_tot.values[:5])
        k = np.argmax(h.hrr_tot.values < 0.5 * h0)
        a, b = h.iloc[k - 1], h.iloc[k]
        tQ = a.t_tauF + (a.hrr_tot - 0.5 * h0) / (a.hrr_tot - b.hrr_tot) * (b.t_tauF - a.t_tauF)
        d = pd.read_csv(os.path.join(run, "hoq_timeseries.csv"))
        late = d[d.t_tauF > 0.5 * m.get("dist_lF", 10.0)]
        t_contact = late.t_tauF[late.Phi_max_star > 0.05].min()
        r = dict(run=os.path.basename(run), gamma=m["case"]["gamma"], Tw=m["T_wall"],
                 unstable=int("planar" not in run), t_contact=t_contact, tQ=tQ, dur=tQ - t_contact,
                 Phi_q=late.Phi_max_star.max())
        for dt in DELAYS:
            for s in ("NH3", "NO", "N2O", "H2"):
                r[f"{s}_{dt:g}"] = np.interp(tQ + dt, d.t_tauF, d[f"X_{s}_y2"])
        rows.append(r)
        print(f"{r['run']:46s} t_contact {t_contact:5.1f} t_Q {tQ:5.1f} | " + " | ".join(
            f"+{dt:g}: NH3 {r[f'NH3_{dt:g}']*100:.2f}% NO {r[f'NO_{dt:g}']*1e6:.0f} N2O {r[f'N2O_{dt:g}']*1e6:.0f}"
            for dt in DELAYS), flush=True)
    out = "emissions_tq.csv"
    new = pd.DataFrame(rows)
    if os.path.exists(out):
        old = pd.read_csv(out)
        new = pd.concat([old[~old.run.isin(new.run)], new])
    new.sort_values(["Tw", "gamma", "unstable", "run"]).to_csv(out, index=False)


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2:])
