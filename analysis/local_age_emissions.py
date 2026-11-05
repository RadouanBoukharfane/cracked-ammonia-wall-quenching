"""Review point 5: near-wall NH3, NO, N2O of crests and troughs compared at the same local time after quenching,
t - t_q(x), instead of the same global time t_Q + 10 tau_F (conditioned_emissions.py).

  python local_age_emissions.py NPROC  -> local_age_emissions.csv, local_age_classes.csv

t_q(x): instant of the local wall-flux maximum (hoq_wallflux.npz, t > 5 tau_F). Classes as in conditioned_emissions.py
(first / last third of t_q). Per plotfile and column: mole fractions averaged over 0 < y < 2 l_F. For each column,
values interpolated at t_q(x) + DT, DT in DTS; planar reference at the same delay after its own flux maximum.
Also: curvature of the crest tip 2 tau_F before t_q per class (local_events_lead2.csv) to check that the temporal
classes are geometric crests and troughs.
"""
import glob, json, os, re, sys
import numpy as np, pandas as pd
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
R = os.path.join(ROOT, "PhaseB", "runs")
MECH = os.path.join(os.environ.get("PELE_HOME", "PeleLMeX"), "Submodules/PelePhysics/Mechanisms/Mei2021NH3/mechanism.yaml")
SP = ("NH3", "NO", "N2O")
DTS = (1.0, 2.0, 3.0, 4.0)


def one(args):
    p, lF = args
    import yt, cantera as ct
    yt.set_log_level(50)
    gas = ct.Solution(MECH)
    W = gas.molecular_weights
    ds = yt.load(p)
    lev = ds.index.max_level
    n = (ds.domain_dimensions * ds.refine_by**lev).astype(int)
    dy = (ds.domain_width / n).d[1]
    ny = int(np.ceil(2.0 * lF / dy))
    cg = ds.covering_grid(lev, left_edge=ds.domain_left_edge, dims=[n[0], ny, 1])
    Y = np.stack([np.clip(cg["boxlib", f"Y({s})"].d[:, :, 0], 0, None) for s in gas.species_names], -1)
    X = (Y / W) / (Y / W).sum(-1, keepdims=True)
    return float(ds.current_time), {s: X[:, :, gas.species_index(s)].mean(1) for s in SP}


def plts(run):
    out = []
    for p in glob.glob(os.path.join(run, "plt*")):
        if re.fullmatch(r"plt\d+", os.path.basename(p)) and "Y(NH3)" in open(os.path.join(p, "Header")).read():
            out.append(p)
    return sorted(out)


def series(run, pool):
    m = json.load(open(os.path.join(run, "case.json")))
    res = sorted(pool.map(one, [(p, m["l_F"]) for p in plts(run)]), key=lambda a: a[0])
    t = np.array([a[0] for a in res]) / m["tau_F"]
    return m, t, {s: np.array([a[1][s] for a in res]) for s in SP}


def main(nproc):
    from multiprocessing import Pool
    rows, cls_rows = [], []
    ev = pd.read_csv(os.path.join(ROOT, "analysis", "local_events_lead2.csv"))
    with Pool(nproc) as pool:
        for g, tw in ((0.2, 750), (0.43, 750), (0.75, 750), (0.43, 500)):
            sfx = "" if tw == 750 else "_Tw500"
            pr = os.path.join(R, f"HOQ_g{g:.2f}_phi0.50_P10.0_T750_planar{sfx}")
            m, t, X = series(pr, pool)
            z = np.load(os.path.join(pr, "hoq_wallflux.npz")); tz = z["t"] / m["tau_F"]
            k = tz > 0.5 * m.get("dist_lF", 30.0); tqp = tz[k][np.argmax(z["phi_w"][k].mean(1))]
            ref = {s: {dt: np.interp(tqp + dt, t, X[s].mean(1)) for dt in DTS} for s in SP}
            for run in sorted(glob.glob(os.path.join(R, f"HOQ_g{g:.2f}_phi0.50_P10.0_T750_t??.0{sfx}"))):
                if not os.path.exists(os.path.join(run, "hoq_wallflux.npz")):
                    continue
                m, t, X = series(run, pool)
                z = np.load(os.path.join(run, "hoq_wallflux.npz")); tz = z["t"] / m["tau_F"]; o = np.argsort(tz)
                tz, phi = tz[o], z["phi_w"][o]; k = tz > 5.0
                tq = tz[k][np.argmax(phi[k], axis=0)]
                nx = len(tq); nxp = X["NH3"].shape[1]
                tq = np.interp((np.arange(nxp) + 0.5) / nxp, (np.arange(nx) + 0.5) / nx, tq)   # to plotfile columns
                q1, q2 = np.quantile(tq, [1 / 3, 2 / 3])
                cls = np.where(tq <= q1, "crest", np.where(tq >= q2, "trough", "middle"))
                for dt in DTS:
                    ok = tq + dt <= t[-1]
                    for c in ("crest", "trough"):
                        sel = ok & (cls == c)
                        if sel.sum() < 0.5 * (cls == c).sum():
                            continue
                        row = dict(run=os.path.basename(run), gamma=g, Tw=tw, cls=c, dt=dt, n=int(sel.sum()))
                        for s in SP:
                            v = [np.interp(tq[j] + dt, t, X[s][:, j]) for j in np.where(sel)[0]]
                            row[s] = np.mean(v) / ref[s][dt]
                        rows.append(row)
                e = ev[ev.run == os.path.basename(run)]
                if len(e):
                    ce = np.interp(e.x_lF.values / (m["Lx"] / m["l_F"]), (np.arange(nxp) + 0.5) / nxp, tq)
                    eq1, eq2 = np.quantile(ce, [1 / 3, 2 / 3])
                    for c, msk in (("crest", ce <= eq1), ("trough", ce >= eq2)):
                        cls_rows.append(dict(run=os.path.basename(run), gamma=g, Tw=tw, cls=c,
                                             kappa_med=np.nanmedian(e.kappa_lF.values[msk]),
                                             frac_convex=np.nanmean(e.kappa_lF.values[msk] > 0),
                                             phi_loc_med=np.nanmedian(e.phi_loc.values[msk])))
                print(os.path.basename(run), "done", flush=True)
    d = pd.DataFrame(rows); d.to_csv(os.path.join(ROOT, "analysis", "local_age_emissions.csv"), index=False)
    c = pd.DataFrame(cls_rows); c.to_csv(os.path.join(ROOT, "analysis", "local_age_classes.csv"), index=False)
    pd.set_option("display.width", 200)
    print(d.groupby(["gamma", "Tw", "dt", "cls"])[list(SP)].mean().unstack("cls").round(2).to_string())
    print(c.groupby(["gamma", "Tw", "cls"])[["kappa_med", "frac_convex", "phi_loc_med"]].mean().round(3).to_string())


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
