"""Local mixture of the crests that set the wall heat flux, measured just before they touch the wall.

  python crest_at_quench.py NPROC RUN [RUN ...]   (needs RUN/hoq_wallflux.npz from post_hoq.py)

Events: every plotfile time t* at which the maximum wall heat flux exceeds FRAC x its global maximum; x* = position
of that maximum. For each event and each lead time dt in LEADS (tau_F), on the plotfile at t* - dt, the T_track
isoline (PhaseB0/local_mixture.py definitions: phi_loc = n_H/(2 n_O), gamma_loc = 2 X_H2/(3 X_NH3 + 2 X_H2)) is
sampled within |x - x*| < 1 l_F, and the point closest to the wall (the crest tip) is kept.
Output: RUN/crest_at_quench.csv and a summary line per run (mean and max of the tip mixture over the events).
"""
import glob
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "PhaseB0"))
from local_mixture import isoline_samples, MECH_ROOT  # noqa: E402

FRAC, LEADS, YMAX = 0.7, (4.0, 3.0, 2.0, 1.0), 14.0


def plt_list(run, tau):
    out = []
    for p in glob.glob(os.path.join(run, "plt*")):
        if not re.fullmatch(r"plt\d+", os.path.basename(p)):
            continue
        L = open(os.path.join(p, "Header")).read().split("\n")
        out.append((float(L[2 + int(L[1]) + 1]) / tau, p))
    return sorted(out)


def tip(args):
    p, m, xs = args
    import cantera as ct
    import yt
    yt.set_log_level(50)
    gas = ct.Solution(os.path.join(MECH_ROOT, m["mech"], "mechanism.yaml"))
    W = gas.molecular_weights
    nH = np.array([gas.n_atoms(s, "H") for s in gas.species_names]) / W
    nO = np.array([gas.n_atoms(s, "O") for s in gas.species_names]) / W
    iH2, iNH3 = gas.species_index("H2"), gas.species_index("NH3")
    lF = m["l_F"]
    ds = yt.load(p)
    lev = ds.index.max_level
    n = (ds.domain_dimensions * ds.refine_by**lev).astype(int)
    dx = (ds.domain_width / n).d[0]
    ny = int(np.ceil(YMAX * lF / dx))
    cg = ds.covering_grid(lev, left_edge=ds.domain_left_edge, dims=[n[0], ny, 1])
    T = cg["boxlib", "temp"].d[:, :, 0]
    Y = np.stack([cg["boxlib", f"Y({s})"].d[:, :, 0] for s in gas.species_names], axis=-1)
    phi = (Y @ nH) / (2 * (Y @ nO))
    X = Y / W
    X /= X.sum(axis=-1, keepdims=True)
    gam = 2 * X[..., iH2] / np.maximum(3 * X[..., iNH3] + 2 * X[..., iH2], 1e-30)
    s = isoline_samples(T, {"phi": phi, "gamma": gam}, m["T_track"], dx)
    res = []
    for x0 in xs:
        if not s:
            res.append((np.nan,) * 4); continue
        dxp = np.abs(((s["x"] - x0) + 0.5 * m["Lx"]) % m["Lx"] - 0.5 * m["Lx"])   # periodic distance
        k = np.where(dxp < lF)[0]
        if len(k) == 0:
            res.append((np.nan,) * 4); continue
        j = k[np.argmin(s["y"][k])]
        res.append((s["y"][j] / lF, s["kappa"][j] * lF, s["phi"][j], s["gamma"][j]))
    return p, res


def main(nproc, runs):
    from multiprocessing import Pool
    import pandas as pd
    jobs, meta = [], {}
    for run in runs:
        m = json.load(open(os.path.join(run, "case.json")))
        z = np.load(os.path.join(run, "hoq_wallflux.npz"))
        t, pw = z["t"] / m["tau_F"], z["phi_w"]
        nx = pw.shape[1]
        gmax = pw.max()
        ev = [(t[i], (np.argmax(pw[i]) + 0.5) * m["Lx"] / nx, pw[i].max() / z["q_ref"])
              for i in range(len(t)) if pw[i].max() >= FRAC * gmax]
        plts = plt_list(run, m["tau_F"])
        tp = np.array([a for a, _ in plts])
        req = {}
        for (ts, xs, f) in ev:
            for dt in LEADS:
                p = plts[int(np.argmin(np.abs(tp - (ts - dt))))][1]
                req.setdefault(p, []).append((ts, dt, xs, f))
        meta[run] = (m, req)
        jobs += [(p, m, [r[2] for r in rr]) for p, rr in req.items()]
    with Pool(nproc) as pool:
        out = dict(pool.map(tip, jobs))
    for run in runs:
        m, req = meta[run]
        rows = []
        for p, rr in req.items():
            for (ts, dt, xs, f), (yt_, kap, ph, ga) in zip(rr, out[p]):
                rows.append(dict(t_event=ts, lead=dt, x_lF=xs / m["l_F"], Phi_star=f, y_tip_lF=yt_, kappa_lF=kap,
                                 phi_loc=ph, gamma_loc=ga))
        d = pd.DataFrame(rows).sort_values(["lead", "t_event"])
        d.to_csv(os.path.join(run, "crest_at_quench.csv"), index=False)
        for dt in LEADS:
            q = d[d.lead == dt].dropna()
            print(f"{os.path.basename(run)} lead {dt:.1f} tau_F: {len(q)} events, y_tip {q.y_tip_lF.mean():.1f} l_F, "
                  f"kappa {q.kappa_lF.mean():+.2f} | phi_loc mean {q.phi_loc.mean():.3f} max {q.phi_loc.max():.3f} | "
                  f"gamma_loc mean {q.gamma_loc.mean():.3f} max {q.gamma_loc.max():.3f}", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2:])
