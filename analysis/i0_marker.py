"""Review point 1b: sensitivity of the reactivity factor I0 to the flame-surface marker (saturated B0 flames).

  python i0_marker.py NPROC   -> i0_marker.csv (per plotfile), printed means over 20-40 tau_F

A = length of the isoline of: T = T_track (paper), c_O2 = (Y_O2,u - Y_O2)/(Y_O2,u - Y_O2,b) = 0.5 and 0.9, and the
fuel progress variables of Schneider et al. / Howarth & Aspden, PV_H2 = 1 - Y_H2/Y_H2,u = 0.9 (gamma > 0) and
PV_NH3 = 1 - Y_NH3/Y_NH3,u = 0.9 (gamma < 1). s_c from the O2 consumption as in post_b0.py.
"""
import glob, json, os, re, sys
import numpy as np, pandas as pd, contourpy, cantera as ct
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
MECH_ROOT = os.path.join(os.environ.get("PELE_HOME", "PeleLMeX"), "Submodules/PelePhysics/Mechanisms")
RUNS = ["B0_g0.20_phi0.50_P10.0_T750", "B0_g0.43_phi0.50_P10.0_T750", "B0_g0.75_phi0.50_P10.0_T750",
        "B0_g1.00_phi0.50_P10.0_T750", "B0_g0.43_phi0.40_P10.0_T750", "B0_g0.43_phi0.50_P10.0_T500",
        "B0_g0.60_phi0.50_P10.0_T500"]


def length(z, lev, dx):
    lines = contourpy.contour_generator(z=z.T, line_type="Separate").lines(lev)
    return sum(np.sum(np.hypot(*np.diff(l, axis=0).T)) for l in lines) * dx


def one(args):
    p, m, ref = args
    import yt
    yt.set_log_level(50)
    ds = yt.load(p)
    lev = ds.index.max_level
    n = (ds.domain_dimensions * ds.refine_by**lev).astype(int)
    cg = ds.covering_grid(lev, left_edge=ds.domain_left_edge, dims=n)
    dx = (ds.domain_width / n).d[0]
    f = lambda k: cg["boxlib", k].d[:, :, 0]
    T, wO2, YO2 = f("temp"), f("I_R(O2)"), f("Y(O2)")
    sc = -wO2.sum() * dx * dx / (ref["rho_u"] * (ref["YO2_u"] - ref["YO2_b"]) * m["Lx"]) / m["sL"]
    r = {"t_tauF": float(ds.current_time) / m["tau_F"], "sc_sL": sc}
    cO2 = (ref["YO2_u"] - YO2) / (ref["YO2_u"] - ref["YO2_b"])
    A = {"T_track": length(T, m["T_track"], dx), "cO2_0.5": length(cO2, 0.5, dx), "cO2_0.9": length(cO2, 0.9, dx)}
    if ref["YH2_u"] > 1e-6:
        A["PV_H2_0.9"] = length(1 - f("Y(H2)") / ref["YH2_u"], 0.9, dx)
    if ref["YNH3_u"] > 1e-6:
        A["PV_NH3_0.9"] = length(1 - f("Y(NH3)") / ref["YNH3_u"], 0.9, dx)
    for k, v in A.items():
        r[f"A_{k}"] = v / m["Lx"]; r[f"I0_{k}"] = sc / (v / m["Lx"])
    return r


def main(nproc):
    from multiprocessing import Pool
    rows = []
    for run in RUNS:
        d = os.path.join(ROOT, "PhaseB0", "runs", run)
        m = json.load(open(os.path.join(d, "case.json")))
        gas = ct.Solution(os.path.join(MECH_ROOT, m["mech"], "mechanism.yaml"))
        gas.TPY = m["Tu"], m["P_Pa"], m["Y_u"]
        ref = dict(rho_u=gas.density, YO2_u=gas["O2"].Y[0], YH2_u=gas["H2"].Y[0], YNH3_u=gas["NH3"].Y[0])
        gas.equilibrate("HP"); ref["YO2_b"] = gas["O2"].Y[0]
        t0 = 30.0 if "T500" in run and "g0.43" in run else (20.0 if "T500" in run else (15.0 if "phi0.40" in run else 20.0))
        plts = []
        cand = glob.glob(os.path.join(d, "species", "from_chk*", "plt_sp*")) or \
            [p for p in glob.glob(os.path.join(d, "plt*")) if re.fullmatch(r"plt\d+", os.path.basename(p))]
        for p in cand:
            if "Y(O2)" not in open(os.path.join(p, "Header")).read():
                continue
            L = open(os.path.join(p, "Header")).read().split("\n")
            if float(L[2 + int(L[1]) + 1]) / m["tau_F"] >= t0 - 1e-6:
                plts.append(p)
        with Pool(nproc) as pool:
            res = pool.map(one, [(p, m, ref) for p in sorted(plts)])
        for r in res:
            r["run"] = run
        rows += res
        q = pd.DataFrame(res)
        print(run, f"n={len(q)} t>={t0}: " + "  ".join(f"{c[3:]}: {q[c].mean():.3f}" for c in q.columns if c.startswith("I0_")),
              flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(ROOT, "analysis", "i0_marker.csv"), index=False)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
