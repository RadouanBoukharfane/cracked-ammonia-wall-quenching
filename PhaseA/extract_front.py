"""Front tracking and dispersion relation of a Phase A run (Lehmann 2024, section 3.2).

  python extract_front.py RUN_DIR [--keep-plt] [--a-lin-factor 50]

1. For each plotfile: finest level in the refined band, front y_f(x) = first crossing of
   T = T_track (case.json) along y, linear interpolation -> front.npz (t, x, y_f).
2. A_i(t) = 2 |FFT(y_f)|_i / N for the NMODES perturbed modes, omega_i = slope of ln A_i over the
   linear window [t_start, t_lin], t_lin = first time max_i A_i > 50 A0 (= 5e-5 l_F for A0 = 1e-6 l_F: beyond,
   mode coupling feeds the damped modes, which then seem to grow; seen on the 10 bar validation),
   t_start = 1 tau_F (initial transient). Per mode, samples below the noise floor (A < 1e-8 l_F) are
   dropped; at least 4 samples are needed. Normalised: omega_bar = omega tau_F, k_bar = k l_F.
3. Darrieus-Landau reference omega_DL_bar = k_bar s/(s+1) (sqrt(s + 1 - 1/s) - 1), s = rho_u/rho_b.
-> dispersion.csv. Plotfiles are deleted after extraction unless --keep-plt (the first and last are kept).
"""
import argparse
import glob
import json
import os
import re
import shutil

import numpy as np

A_FLOOR = 1e-8  # noise floor of the front amplitude, in l_F
NOISE_MAX = 1e-8  # quality flag: max amplitude of unseeded modes 41..100 at 2 tau_F (clean runs: 1e-10 to 2e-9 l_F)


def noise_level(t, yf, meta):
    """Largest amplitude (in l_F) of the unseeded modes 41..100 near t = 2 tau_F: numerical forcing."""
    if t[-1] < 2.0 * meta["tau_F"]:
        return float("nan")
    i = int(np.argmin(np.abs(t - 2.0 * meta["tau_F"])))
    a = 2 * np.abs(np.fft.rfft(yf[i] - yf[i].mean())) / yf.shape[1] / meta["l_F"]
    return float(a[meta["nmodes"] + 1:101].max())


def front_from_plotfile(plt, meta):
    import yt
    yt.set_log_level(50)
    ds = yt.load(plt)
    lF = meta["l_F"]
    lev = ds.index.max_level
    dx = ds.index.get_smallest_dx().d
    y_lo, y_hi = 7.2 * lF, 12.8 * lF                       # inside the static band
    nx = int(round(meta["Lx"] / dx))
    ny = int((y_hi - y_lo) / dx)
    cg = ds.covering_grid(lev, left_edge=[0.0, y_lo, 0.0], dims=[nx, ny, 1])
    T = cg["boxlib", "temp"].d[:, :, 0]
    y = y_lo + (np.arange(ny) + 0.5) * dx
    Tt = meta["T_track"]
    yf = np.full(nx, np.nan)
    for i in range(nx):
        c = T[i]
        j = np.where((c[:-1] < Tt) & (c[1:] >= Tt))[0]
        if len(j):
            j = j[0]
            yf[i] = y[j] + (Tt - c[j]) / (c[j + 1] - c[j]) * dx
    x = (np.arange(nx) + 0.5) * dx
    return float(ds.current_time.d), x, yf


def dispersion(t, x, yf, meta, a_lin):
    lF, tau, n = meta["l_F"], meta["tau_F"], meta["nmodes"]
    N = yf.shape[1]
    A = 2 * np.abs(np.fft.rfft(yf - yf.mean(axis=1, keepdims=True), axis=1))[:, 1:n + 1] / N
    k = 2 * np.pi * np.arange(1, n + 1) / meta["Lx"]
    t_start = 1.0 * tau
    above = np.where(A.max(axis=1) > a_lin * meta["A0"])[0]
    t_lin = t[above[0]] if len(above) else t[-1]
    w = (t >= t_start) & (t <= t_lin)
    omega = np.full(n, np.nan)
    r2 = np.full(n, np.nan)
    if w.sum() >= 4:
        for i in range(n):
            wi = w & (A[:, i] > A_FLOOR * lF)
            if wi.sum() < 4:
                continue
            lnA = np.log(A[wi, i])
            p = np.polyfit(t[wi], lnA, 1)
            omega[i] = p[0]
            res = lnA - np.polyval(p, t[wi])
            r2[i] = 1 - res.var() / max(lnA.var(), 1e-300)
    s = meta["rho_u"] / meta["rho_b"]
    kb = k * lF
    om_dl = kb * s / (s + 1) * (np.sqrt(s + 1 - 1 / s) - 1)
    return k, A, omega, r2, om_dl, t_start, t_lin, int(w.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--keep-plt", action="store_true")
    ap.add_argument("--a-lin-factor", dest="a_lin", type=float, default=50.0,
                    help="end of the linear window when max A > factor * A0")
    a = ap.parse_args()
    meta = json.load(open(os.path.join(a.run, "case.json")))
    pmf = np.loadtxt(os.path.join(a.run, "pmf.dat"), skiprows=2)
    meta["rho_u"], meta["rho_b"] = pmf[0, 3], pmf[-1, 3]

    npz = os.path.join(a.run, "front.npz")
    done = dict(np.load(npz)) if os.path.exists(npz) else {"t": np.zeros(0), "yf": None, "x": None}
    # PeleLMeX renames a plotfile overwritten at restart to pltNNNNN.old.<id>: ignore those
    plts = sorted((p for p in glob.glob(os.path.join(a.run, "plt*")) if re.fullmatch(r"plt\d+", os.path.basename(p))),
                  key=lambda p: int(os.path.basename(p)[3:]))
    ts, yfs = list(done["t"]), [] if done["yf"] is None else list(done["yf"])
    x = done["x"]
    for p in plts:
        t, x, yf = front_from_plotfile(p, meta)
        if not np.any(np.isclose(ts, t, rtol=0, atol=1e-15)):
            ts.append(t)
            yfs.append(yf)
    order = np.argsort(ts)
    t, yf = np.array(ts)[order], np.array(yfs)[order]
    np.savez(npz, t=t, x=x, yf=yf)
    if not a.keep_plt and len(plts) > 2:
        for p in plts[1:-1]:
            shutil.rmtree(p)

    noise = noise_level(t, yf, meta)
    json.dump({"noise_unseeded_2tau": noise, "noisy": bool(noise > NOISE_MAX)},
              open(os.path.join(a.run, "quality.json"), "w"))
    if noise > NOISE_MAX:
        print(f"WARNING {meta['name']}: unseeded modes at {noise:.1e} l_F at 2 tau_F (clean < {NOISE_MAX:.0e}): "
              f"numerical forcing, dispersion relation not reliable")
    if np.isnan(yf).any():
        print(f"warning: {np.isnan(yf).sum()} front points not found (front left the band?)")
    k, A, om, r2, om_dl, t0, t1, nw = dispersion(t, x, np.nan_to_num(yf, nan=np.nanmean(yf)), meta, a.a_lin)
    lF, tau = meta["l_F"], meta["tau_F"]
    out = np.column_stack([np.arange(1, len(k) + 1), k * lF, om * tau, om_dl, r2, A[0] / lF, A[-1] / lF])
    np.savetxt(os.path.join(a.run, "dispersion.csv"), out, delimiter=",", fmt="%.6e",
               header="mode,k_bar,omega_bar,omega_DL_bar,r2,A0_over_lF,Aend_over_lF", comments="")
    om_ok = np.where(r2 > 0.9, om, np.nan)
    if np.all(np.isnan(om_ok)):
        print(f"{meta['name']}: {len(t)} snapshots up to {t[-1]/tau:.2f} tau_F, no usable linear window yet "
              f"(fit [{t0/tau:.2f}, {t1/tau:.2f}] tau_F, {nw} pts)")
        return
    i = np.nanargmax(om_ok)
    print(f"{meta['name']}: {len(t)} snapshots, fit window [{t0/tau:.2f}, {t1/tau:.2f}] tau_F ({nw} pts), "
          f"omega_bar_max = {om[i]*tau:.3f} at k_bar = {k[i]*lF:.3f} (DL there: {om_dl[i]:.3f})")


if __name__ == "__main__":
    main()
