"""Phase B post-processing of a head-on quenching run (Schneider et al. 2025 Part I quantities + emissions).

  python post_hoq.py RUN_DIR [NPROC]   -> RUN_DIR/hoq_timeseries.csv, RUN_DIR/hoq_wallflux.npz

Per plotfile (finest level, region y < Y_MAX l_F above the wall):
- wall heat flux Phi_w(x) = lambda_w dT/dy|_w, one-sided 2nd order on cell centres
  dT/dy = (-8 T_w + 9 T_1 - T_2) / (3 dy); lambda from Cantera (mixture-averaged) at the first cell state;
- quenching position: y of the maximum heat release per column -> min over x (leading point) and mean;
- near-wall emissions: x-averaged mole fractions of NH3, H2, NO, N2O, NO2 over y < 2 l_F and y < 5 l_F.
Normalisation: Phi* = Phi / (rho_u c_p,u s_L (T_ad - T_u)), y* = y / delta_T (Pe = y / delta_T).
hoq_wallflux.npz also stores, per column x and plotfile, y_hrr (position of the HRR peak) and hrr_col (its value).
"""
import glob
import json
import os
import re
import sys

import numpy as np

MECH_ROOT = os.path.join(os.environ.get("PELE_HOME", "PeleLMeX"), "Submodules/PelePhysics/Mechanisms")
Y_MAX = 12.0
SPECIES = ["NH3", "H2", "NO", "N2O", "NO2"]


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
    dx = (ds.domain_width / n).d[0]
    ny = int(min(n[1], np.ceil(Y_MAX * lF / dx)))
    cg = ds.covering_grid(lev, left_edge=ds.domain_left_edge, dims=[n[0], ny, 1])
    T = cg["boxlib", "temp"].d[:, :, 0]
    hrr = cg["boxlib", "HeatRelease"].d[:, :, 0]
    Y = {s: cg["boxlib", f"Y({s})"].d[:, :, 0] for s in gas.species_names}
    Tw = m["T_wall"]
    dTdy = (-8 * Tw + 9 * T[:, 0] - T[:, 1]) / (3 * dx)
    lam = np.empty(n[0])
    for i in range(n[0]):
        gas.TPY = T[i, 0], m["P_Pa"], {s: max(Y[s][i, 0], 0.0) for s in gas.species_names}
        lam[i] = gas.thermal_conductivity
    phi_w = lam * dTdy                                    # W/m2, > 0: heat into the wall (T increases away from it)
    y = (np.arange(ny) + 0.5) * dx
    jmax = np.argmax(hrr, axis=1)
    y_hrr = y[jmax]
    W = gas.molecular_weights
    Ymat = np.stack([Y[s] for s in gas.species_names], axis=-1)
    X = (Ymat / W) / (Ymat / W).sum(axis=-1, keepdims=True)
    r = {"t": float(ds.current_time), "phi_w_max": phi_w.max(), "phi_w_mean": phi_w.mean(),
         "y_hrr_min": y_hrr.min(), "y_hrr_mean": y_hrr.mean(), "hrr_max": hrr.max()}
    for lim in (2.0, 5.0):
        k = y < lim * lF
        for s in SPECIES:
            r[f"X_{s}_y{lim:g}"] = X[:, k, gas.species_index(s)].mean()
    return r, phi_w, y_hrr, hrr.max(axis=1)


def main(run, nproc=1):
    import cantera as ct
    from multiprocessing import Pool
    m = json.load(open(os.path.join(run, "case.json")))
    gas = ct.Solution(os.path.join(MECH_ROOT, m["mech"], "mechanism.yaml"))
    gas.TPY = m["Tu"], m["P_Pa"], m["Y_u"]
    q_ref = gas.density * gas.cp_mass * m["sL"] * (m["Tb_eq"] - m["Tu"])
    plts = sorted((p for p in glob.glob(os.path.join(run, "plt*")) if re.fullmatch(r"plt\d+", os.path.basename(p))),
                  key=lambda p: int(os.path.basename(p)[3:]))
    with Pool(nproc) as pool:
        out = pool.map(one, [(p, m) for p in plts])
    import pandas as pd
    df = pd.DataFrame([o[0] for o in out]).sort_values("t").drop_duplicates("t")
    df["t_tauF"] = df.t / m["tau_F"]
    df["Phi_max_star"] = df.phi_w_max / q_ref
    df["Phi_mean_star"] = df.phi_w_mean / q_ref
    df["Pe_min"] = df.y_hrr_min / m["l_F"]
    df["Pe_mean"] = df.y_hrr_mean / m["l_F"]
    df.to_csv(os.path.join(run, "hoq_timeseries.csv"), index=False)
    np.savez(os.path.join(run, "hoq_wallflux.npz"), t=np.array([o[0]["t"] for o in out]),
             phi_w=np.array([o[1] for o in out]), q_ref=q_ref,
             y_hrr=np.array([o[2] for o in out]), hrr_col=np.array([o[3] for o in out]))   # per column (x)
    # T_w < T_u: the wall/gas temperature jump gives a spurious flux peak at start-up -> search the quenching
    # peak only once the flame is close (t > dist/2 tau_F, planar runs)
    late = df[df.t_tauF > 0.5 * m.get("dist_lF", 10.0)]     # unstable runs (no dist_lF): skip the first 5 tau_F
    i = late.Phi_max_star.idxmax()
    print(f"{m['name']}: q_ref = {q_ref/1e6:.3f} MW/m2 ; Phi_q* = {df.Phi_max_star[i]:.3f} at t = "
          f"{df.t_tauF[i]:.2f} tau_F ; Pe_q (min over time of leading y_HRR/delta_T) = {df.Pe_min.min():.3f} ; "
          f"Pe at Phi_q: {df.Pe_min[i]:.3f}")
    print(df[["t_tauF", "Phi_max_star", "Phi_mean_star", "Pe_min", "Pe_mean", "X_NH3_y2", "X_NO_y2", "X_N2O_y2",
              "X_H2_y2"]].to_string(index=False, float_format="%.4g"))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1)
