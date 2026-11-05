"""Cold-wall series (planar HOQ): near-wall profiles, cold-layer metrics and NO reaction budget.

  python coldwall_analysis.py NPROC RUN [RUN ...]   -> RUN/coldwall.npz

Per plotfile in [t_q - 22, t_q + 14] tau_F (planar runs: x-averaged profiles, finest level, y < 40 l_F):
- profiles T, rho, HRR, mole fractions;
- NO, N2O, NO2, NH3 net molar production rate of every reaction, integrated over 0 < y < 2 l_F and 0 < y < 5 l_F
  (kmol m^-2 s^-1), and the NO inventory int rho Y_NO dy over the same layers (transport = residual of the budget).
t_q = time of the maximum wall heat flux once the flame is close (t > dist/2 tau_F), from hoq_timeseries.csv.
"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

MECH_ROOT = os.path.join(os.environ.get("PELE_HOME", "PeleLMeX"), "Submodules/PelePhysics/Mechanisms")
Y_MAX = 40.0
LAYERS = (2.0, 5.0)
BUDGET = ("NO", "N2O", "NO2", "NH3")


def one(args):
    p, m = args
    import cantera as ct
    import yt
    yt.set_log_level(50)
    gas = ct.Solution(os.path.join(MECH_ROOT, m["mech"], "mechanism.yaml"))
    lF = m["l_F"]
    ds = yt.load(p)
    lev = ds.index.max_level
    n = (ds.domain_dimensions * ds.refine_by**lev).astype(int)
    dy = (ds.domain_width / n).d[1]
    ny = int(min(n[1], np.ceil(Y_MAX * lF / dy)))
    cg = ds.covering_grid(lev, left_edge=ds.domain_left_edge, dims=[n[0], ny, 1])
    T = cg["boxlib", "temp"].d[:, :, 0].mean(0)
    rho = cg["boxlib", "density"].d[:, :, 0].mean(0)
    hrr = cg["boxlib", "HeatRelease"].d[:, :, 0].mean(0)
    Y = np.stack([np.clip(cg["boxlib", f"Y({s})"].d[:, :, 0].mean(0), 0, None) for s in gas.species_names], -1)
    Y /= Y.sum(-1, keepdims=True)
    y = (np.arange(ny) + 0.5) * dy
    st = ct.SolutionArray(gas, ny)
    st.TPY = T, m["P_Pa"], Y
    X = st.X
    q = st.net_rates_of_progress                                                    # (ny, nreac) kmol/m3/s
    r = {"t": float(ds.current_time)}
    for s_ in BUDGET:
        i_ = gas.species_index(s_)
        nu = (gas.product_stoich_coeffs[i_, :] - gas.reactant_stoich_coeffs[i_, :])  # (nreac,)
        for lim in LAYERS:
            k = y < lim * lF
            r[f"w{s_}_y{lim:g}"] = (q[k] * nu).sum(0) * dy        # per reaction, kmol/m2/s
            r[f"inv{s_}_y{lim:g}"] = (rho[k] * Y[k, i_]).sum() * dy / gas.molecular_weights[i_]   # kmol/m2
    return r, dict(T=T, rho=rho, hrr=hrr, X=X, y=y)


def main(nproc, runs):
    from multiprocessing import Pool
    import cantera as ct
    tasks, meta = [], {}
    for run in runs:
        m = json.load(open(os.path.join(run, "case.json")))
        d = pd.read_csv(os.path.join(run, "hoq_timeseries.csv"))
        late = d[d.t_tauF > 0.5 * m.get("dist_lF", 30.0)]
        tq = float(d.t_tauF[late.Phi_max_star.idxmax()])
        plts = sorted((p for p in glob.glob(os.path.join(run, "plt*")) if re.fullmatch(r"plt\d+", os.path.basename(p))),
                      key=lambda p: int(os.path.basename(p)[3:]))
        sel = []
        for p in plts:
            L = open(os.path.join(p, "Header")).read().split("\n")
            t = float(L[2 + int(L[1]) + 1]) / m["tau_F"]
            if tq - 22.01 <= t <= tq + 14.01:
                sel.append(p)
        meta[run] = (m, tq, len(sel))
        tasks += [(p, m) for p in sel]
    with Pool(nproc) as pool:
        out = pool.map(one, tasks)
    i = 0
    gas = ct.Solution(os.path.join(MECH_ROOT, "Mei2021NH3", "mechanism.yaml"))
    for run in runs:
        m, tq, nsel = meta[run]
        res = out[i:i + nsel]; i += nsel
        o = sorted(res, key=lambda a: a[0]["t"])
        save = {"t_tauF": np.array([a[0]["t"] for a in o]) / m["tau_F"], "tq": tq, "y_lF": o[0][1]["y"] / m["l_F"],
                "species": np.array(gas.species_names), "reactions": np.array([R.equation for R in gas.reactions()])}
        for k in ("T", "rho", "hrr", "X"):
            save[k] = np.array([a[1][k] for a in o])
        for s_ in BUDGET:
            for lim in LAYERS:
                save[f"w{s_}_y{lim:g}"] = np.array([a[0][f"w{s_}_y{lim:g}"] for a in o])
                save[f"inv{s_}_y{lim:g}"] = np.array([a[0][f"inv{s_}_y{lim:g}"] for a in o])
        np.savez_compressed(os.path.join(run, "coldwall.npz"), **save)
        print(run, "t_q =", tq, "plotfiles:", nsel, flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2:])
