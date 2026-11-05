"""Sanity checks of the Phase A wave-2 runs (setup and results), read-only.

  python wave2_sanity.py NPROC   -> wave2_sanity_done.csv, wave2_sanity_running.csv, printed flags

Finished runs (dispersion.csv): noise flag (quality.json), fit quality at the peak, peak inside the resolved k range,
small-k limit against Darrieus-Landau, number of samples in the linear window, smoothness along gamma within each
(phi, P, T_u) row, residual against the omega_max ~ X_b regression of all wave-1 + wave-2 cases.
Running runs: temporals (max T against the 1D adiabatic value, integrated heat release drift), and the front on the
first, middle and last plotfiles (front inside the extraction band 7.2-12.8 l_F, mean drift, growth of the seeded
modes, unseeded-mode noise). Nothing is written in the run directories.
"""
import glob
import json
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
PA = os.path.join(ROOT, "PhaseA")
sys.path.insert(0, PA)
from extract_front import front_from_plotfile  # noqa: E402

OUT = os.path.join(ROOT, "analysis")


def wave2_names():
    return [l.split()[1] for l in open(os.path.join(PA, "cases", "wave2_jobs.txt"))]


def running():
    q = subprocess.run(["squeue", "-u", os.environ["USER"], "-h", "-t", "R", "-o", "%j"], capture_output=True, text=True)
    return set(q.stdout.split())


def done_checks(names):
    rows = []
    for n in names:
        r = os.path.join(PA, "runs", n)
        f = os.path.join(r, "dispersion.csv")
        if not os.path.exists(f):
            continue
        m = json.load(open(os.path.join(r, "case.json")))
        d = pd.read_csv(f)
        q = json.load(open(os.path.join(r, "quality.json"))) if os.path.exists(os.path.join(r, "quality.json")) else {}
        dd = d[(d.r2 > 0.9) | (d.omega_bar > 0)].reset_index(drop=True)
        i = dd.omega_bar.idxmax()
        fr = np.load(os.path.join(r, "front.npz"))
        t = fr["t"]; yf = fr["yf"]
        A = 2 * np.abs(np.fft.rfft(yf - np.nanmean(yf, 1, keepdims=True), axis=1))[:, 1:m["nmodes"] + 1] / yf.shape[1]
        above = np.where(np.nanmax(A, 1) > 50 * m["A0"])[0]
        tlin = t[above[0]] if len(above) else t[-1]
        nwin = int(((t >= m["tau_F"]) & (t <= tlin)).sum())
        c = m["case"]
        rows.append(dict(name=n, gamma=c["gamma"], phi=c["phi"], P=c["P_atm"], Tu=c["T_u"], omega_max=dd.omega_bar[i],
                         k_max=dd.k_bar[i], mode_max=int(dd["mode"][i]), r2_peak=dd.r2[i],
                         n_good=int(((d.r2 > 0.9) & (d.omega_bar > 0)).sum()),
                         om1_over_DL=d.omega_bar[0] / d.omega_DL_bar[0], noisy=q.get("noisy", None),
                         noise=q.get("noise_unseeded_2tau", np.nan), nwin=nwin, tlin_tau=tlin / m["tau_F"],
                         front_nan=float(np.isnan(yf).mean()), sL=m["sL"], sL_case=c["sL"]))
    D = pd.DataFrame(rows)
    # smoothness along gamma in each (phi, P, Tu) row: deviation from the mean of the two neighbours
    D["dev_neigh"] = np.nan
    for _, g in D.groupby(["phi", "P", "Tu"]):
        g = g.sort_values("gamma")
        v = g.omega_max.values
        for j in range(1, len(g) - 1):
            D.loc[g.index[j], "dev_neigh"] = v[j] / (0.5 * (v[j - 1] + v[j + 1])) - 1
    # regression residual on the whole database
    S = pd.read_csv(os.path.join(OUT, "stability_database.csv"))
    ok = S.omega_max > 0.05
    p = np.polyfit(np.log(S.X_b[ok]), np.log(S.omega_max[ok]), 1)
    S["res"] = np.log(S.omega_max) - np.polyval(p, np.log(S.X_b))
    mad = 1.4826 * np.median(np.abs(S.res[ok] - np.median(S.res[ok])))
    D = D.merge(S[["name", "X_b", "res"]], on="name", how="left")
    D["res_sigma"] = D.res / mad
    D["flags"] = ""
    def flag(mask, txt):
        D.loc[mask, "flags"] += txt + ";"
    flag(D.noisy == True, "noisy")  # noqa: E712
    flag(D.r2_peak < 0.9, "poor-fit-peak")
    flag((D.mode_max <= 1) | (D.mode_max >= 40), "peak-at-k-edge")
    flag(D.nwin < 5, "short-window")
    flag((D.om1_over_DL > 2) | (D.om1_over_DL < 0), "small-k-vs-DL")
    flag(D.front_nan > 0.01, "front-lost")
    flag(D.dev_neigh.abs() > 0.3, "not-smooth-in-gamma")
    flag(D.res_sigma.abs() > 3, "outlier-vs-database")
    flag((D.sL / D.sL_case - 1).abs() > 0.02, "sL-mismatch")
    return D


def run_check(args):
    n, plts = args
    r = os.path.join(PA, "runs", n)
    m = json.load(open(os.path.join(r, "case.json")))
    out = dict(name=n)
    try:
        te = pd.read_csv(os.path.join(r, "temporals", "tempExtremas"), usecols=["time", "max_temp"])
        ts = pd.read_csv(os.path.join(r, "temporals", "tempState"), usecols=["time", "heatRelease"])
        out["Tmax_over_Tad"] = te.max_temp.iloc[-1] / m["T_end_1d"]
        out["Tmax_peak_over_Tad"] = te.max_temp.max() / m["T_end_1d"]
        h0 = ts.heatRelease.iloc[: max(5, len(ts) // 50)].median()
        out["HR_last_over_0"] = ts.heatRelease.iloc[-1] / h0
        out["HR_max_over_0"] = ts.heatRelease.max() / h0
        out["t_tau"] = ts.time.iloc[-1] / m["tau_F"]
    except Exception as e:  # noqa: BLE001
        out["temporal_err"] = str(e)[:60]
    fronts = []
    for p in plts:
        try:
            t, x, yf = front_from_plotfile(p, m)
            fronts.append((t, yf))
        except Exception as e:  # noqa: BLE001
            out["front_err"] = str(e)[:60]
    if fronts:
        lF = m["l_F"]
        for lab, (t, yf) in zip(("first", "mid", "last"), fronts):
            A = 2 * np.abs(np.fft.rfft(np.nan_to_num(yf - np.nanmean(yf)))) / len(yf) / lF
            out[f"t_{lab}"] = t / m["tau_F"]
            out[f"ymean_{lab}"] = np.nanmean(yf) / lF
            out[f"nan_{lab}"] = float(np.isnan(yf).mean())
            out[f"Aseed_{lab}"] = A[1:m["nmodes"] + 1].max()
            out[f"Anoise_{lab}"] = A[m["nmodes"] + 1:101].max()
    return out


def main(nproc):
    names = wave2_names()
    D = done_checks(names)
    D.to_csv(os.path.join(OUT, "wave2_sanity_done.csv"), index=False)
    pd.set_option("display.width", 250)
    print(f"finished runs checked: {len(D)}; flagged: {(D["flags"] != '').sum()}")
    cols = ["name", "omega_max", "k_max", "mode_max", "r2_peak", "n_good", "om1_over_DL", "nwin", "dev_neigh", "res_sigma",
            "flags"]
    print(D[D["flags"] != ""][cols].to_string(index=False, float_format="%.3g"))
    print("summary finished: omega_max %.2f-%.2f, r2_peak min %.3f, nwin min %d, |dev_neigh| max %.2f, |res| max %.1f sigma"
          % (D.omega_max.min(), D.omega_max.max(), D.r2_peak.min(), D.nwin.min(), D.dev_neigh.abs().max(),
             D.res_sigma.abs().max()))
    run = sorted(running() & set(names))
    tasks = []
    for n in run:
        r = os.path.join(PA, "runs", n)
        pl = sorted((p for p in glob.glob(os.path.join(r, "plt*")) if re.fullmatch(r"plt\d+", os.path.basename(p))),
                    key=lambda p: int(os.path.basename(p)[3:]))
        if len(pl) >= 3:
            pl = [pl[0], pl[len(pl) // 2], pl[-2]]      # pl[-1] may still be being written
        tasks.append((n, pl))
    from multiprocessing import Pool
    with Pool(nproc) as pool:
        R = pd.DataFrame(pool.map(run_check, tasks))
    R["flags"] = ""
    def flag(mask, txt):
        R.loc[mask, "flags"] += txt + ";"
    flag((R.Tmax_peak_over_Tad > 1.05) | (R.Tmax_over_Tad < 0.95), "Tmax-off")
    flag((R.HR_last_over_0 < 0.8) | (R.HR_last_over_0 > 1.5), "HR-drift")
    for lab in ("first", "mid", "last"):
        if f"ymean_{lab}" in R:
            flag((R[f"ymean_{lab}"] < 8.0) | (R[f"ymean_{lab}"] > 12.0), f"front-near-band-edge-{lab}")
            flag(R[f"nan_{lab}"] > 0.01, f"front-lost-{lab}")
            flag(R[f"Anoise_{lab}"] > 1e-8 * 10, f"noise-{lab}")
    if "Aseed_last" in R:
        flag(R.Aseed_last < R.Aseed_first, "no-growth")
    R.to_csv(os.path.join(OUT, "wave2_sanity_running.csv"), index=False)
    print(f"\nrunning runs checked: {len(R)}; flagged: {(R["flags"] != '').sum()}")
    cols = ["name", "t_tau", "Tmax_over_Tad", "HR_last_over_0", "ymean_first", "ymean_last", "Aseed_first", "Aseed_last",
            "Anoise_last", "flags"]
    print(R[[c for c in cols if c in R]].to_string(index=False, float_format="%.3g"))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
