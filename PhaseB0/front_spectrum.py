"""Spectrum of the leading-edge front y_f(x) (first crossing of T_track per column, level-0 grid) in the B0 runs.
  python front_spectrum.py RUN [RUN ...] -> per run: share of the front variance in modes n = 1, 2 (lambda = L_x, L_x/2),
  energy-weighted mean wavelength, at plotfiles with t >= 8 tau_F (every 2 tau_F)."""
import glob, json, os, sys
import numpy as np, yt
yt.set_log_level(50)
for run in sys.argv[1:]:
    m = json.load(open(os.path.join(run, "case.json")))
    rows = []
    for p in sorted(glob.glob(run + "/plt?????")):
        L = open(p + "/Header").read().split("\n"); t = float(L[2 + int(L[1]) + 1]) / m["tau_F"]
        if t < 8 or abs(t - 2 * round(t / 2)) > 0.1:
            continue
        ds = yt.load(p)
        g = ds.covering_grid(1, ds.domain_left_edge, ds.domain_dimensions * [2, 2, 1])
        T = g["boxlib", "temp"].d[:, :, 0]
        y = (np.arange(T.shape[1]) + .5) * ds.domain_width[1].d / T.shape[1] / m["l_F"]
        yf = np.array([y[np.argmax(c > m["T_track"])] for c in T])
        E = np.abs(np.fft.rfft(yf - yf.mean()))**2
        E[0] = 0
        n = np.arange(E.size)
        lam = m["Lx"] / m["l_F"] / np.maximum(n, 1)
        rows.append((t, E[1:3].sum() / E.sum(), (E[1:] * lam[1:]).sum() / E[1:].sum(), yf.std()))
    r = np.array(rows)
    print(f"{os.path.basename(run):32s} t={r[0,0]:.0f}-{r[-1,0]:.0f}: share n=1,2 = {np.round(r[:,1],2)}  "
          f"<lambda>/lF = {np.round(r[:,2],0)}  rms(y_f)/lF = {np.round(r[:,3],1)}", flush=True)
