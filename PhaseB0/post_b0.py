"""Phase B0 post-processing: reactivity factor I0(t) of the unstable free flame (Schneider et al. 2025).

  python post_b0.py RUN_DIR [NPROC]   -> RUN_DIR/b0_timeseries.csv

For each plotfile (finest level, whole domain):
  s_c = -int omega_O2 dA / (rho_u (Y_O2,u - Y_O2,b) L_x)      consumption speed based on O2 (valid for two fuels)
  A   = length of the T = T_track isoline (contourpy), A/L_x = flame-surface ratio
  I0  = (s_c / s_L) / (A / L_x)
plus the inflow velocity from AC_History.dat (global burning speed under active control).
Y_O2,b: HP equilibrium of the unburnt mixture (Cantera, the PelePhysics YAML).
"""
import glob
import json
import os
import re
import sys

import cantera as ct
import contourpy
import numpy as np
import pandas as pd

MECH_ROOT = os.path.join(os.environ.get("PELE_HOME", "PeleLMeX"), "Submodules/PelePhysics/Mechanisms")


def one_plotfile(args):
    p, m, rho_u, YO2_u, YO2_b = args
    import yt
    yt.set_log_level(50)
    sL, lF, tau, Lx = m["sL"], m["l_F"], m["tau_F"], m["Lx"]
    ds = yt.load(p)
    lev = ds.index.max_level
    n = (ds.domain_dimensions * ds.refine_by**lev).astype(int)
    cg = ds.covering_grid(lev, left_edge=ds.domain_left_edge, dims=n)
    T = cg["boxlib", "temp"].d[:, :, 0]
    wO2 = cg["boxlib", "I_R(O2)"].d[:, :, 0]
    dx = (ds.domain_width / n).d[0]
    sc = -wO2.sum() * dx * dx / (rho_u * (YO2_u - YO2_b) * Lx)
    lines = contourpy.contour_generator(z=T.T, line_type="Separate").lines(m["T_track"])
    A = sum(np.sum(np.hypot(*np.diff(l, axis=0).T)) for l in lines) * dx
    ys = np.concatenate([l[:, 1] for l in lines]) * dx if lines else np.array([np.nan])
    r = {"t_tauF": float(ds.current_time) / tau, "sc_sL": sc / sL, "A_Lx": A / Lx,
         "I0": (sc / sL) / (A / Lx), "y_mean_lF": np.mean(ys) / lF, "brush_lF": (ys.max() - ys.min()) / lF}
    print(f"{os.path.basename(p)}: t = {r['t_tauF']:5.1f} tau_F  s_c/s_L = {r['sc_sL']:.3f}  A/Lx = {r['A_Lx']:.3f}  "
          f"I0 = {r['I0']:.3f}  brush = {r['brush_lF']:.1f} l_F", flush=True)
    return r


def main(run, nproc=1):
    from multiprocessing import Pool
    m = json.load(open(os.path.join(run, "case.json")))
    gas = ct.Solution(os.path.join(MECH_ROOT, m["mech"], "mechanism.yaml"))
    gas.TPY = m["Tu"], m["P_Pa"], m["Y_u"]
    rho_u, YO2_u = gas.density, gas.Y[gas.species_index("O2")]
    gas.equilibrate("HP")
    YO2_b = gas.Y[gas.species_index("O2")]
    sL, lF, tau, Lx = m["sL"], m["l_F"], m["tau_F"], m["Lx"]
    plts = sorted((p for p in glob.glob(os.path.join(run, "plt*")) if re.fullmatch(r"plt\d+", os.path.basename(p))),
                  key=lambda p: int(os.path.basename(p)[3:]))
    with Pool(nproc) as pool:
        rows = pool.map(one_plotfile, [(p, m, rho_u, YO2_u, YO2_b) for p in plts])
    df = pd.DataFrame(rows)
    ac = os.path.join(run, "AC_History.dat")
    if os.path.exists(ac):
        h = np.loadtxt(ac)
        df["Vin_sL"] = np.interp(df.t_tauF * tau, h[:, 1], h[:, 2]) / sL
    df.to_csv(os.path.join(run, "b0_timeseries.csv"), index=False)
    print(df.to_string(index=False, float_format="%.3f"))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1)
