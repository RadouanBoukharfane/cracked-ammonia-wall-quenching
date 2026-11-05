"""Quick-look PNGs of the unstable HOQ runs: T, heat release, Y(NH3) (latest plotfile), and a zoom with AMR levels."""
import glob, json, os, sys
import numpy as np, yt, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
yt.set_log_level(50)
runs = sys.argv[1:]
YMAX = 60.0   # l_F

def load(p, m, lev=2):
    ds = yt.load(p); lF = m["l_F"]
    nx = ds.domain_dimensions[0] * 2**lev
    dx = ds.domain_width[0].d / nx
    ny = int(YMAX * lF / dx)
    g = ds.covering_grid(lev, ds.domain_left_edge, [nx, ny, 1])
    f = {k: g["boxlib", k].d[:, :, 0].T for k in ["temp", "HeatRelease", "Y(NH3)", "Y(H2)"]}
    return ds, f, ds.domain_width[0].d / lF

fig, ax = plt.subplots(len(runs), 3, figsize=(14, 3.2 * len(runs)), squeeze=False)
for i, run in enumerate(runs):
    m = json.load(open(os.path.join(run, "case.json")))
    p = sorted(glob.glob(run + "/plt?????"))[-1]
    ds, f, Lx = load(p, m)
    ext = [0, Lx, 0, YMAX]
    t = float(ds.current_time) / m["tau_F"]
    for j, (k, cm, lab) in enumerate([("temp", "inferno", "T [K]"), ("HeatRelease", "magma", "HRR [W/m³]"),
                                      ("Y(NH3)", "viridis", "Y(NH3)")]):
        im = ax[i, j].imshow(f[k], origin="lower", extent=ext, cmap=cm, aspect="equal")
        fig.colorbar(im, ax=ax[i, j], shrink=0.8, label=lab)
        ax[i, j].contour(np.linspace(0, Lx, f["temp"].shape[1]), np.linspace(0, YMAX, f["temp"].shape[0]),
                         f["temp"], [m["T_track"]], colors="c", linewidths=0.5)
        ax[i, j].set_title(f"{m['name'][-5:]}  t = {t:.2f} τ_F  ({os.path.basename(p)})", fontsize=9)
        ax[i, j].set_xlabel("x / l_F"); ax[i, j].set_ylabel("y / l_F (paroi en y = 0)")
fig.tight_layout(); fig.savefig("figs/hoq_unstable_fields.png", dpi=110)

# zoom on the first run: flame structure + AMR level boxes
run = runs[0]; m = json.load(open(os.path.join(run, "case.json"))); lF = m["l_F"]
p = sorted(glob.glob(run + "/plt?????"))[-1]
ds, f, Lx = load(p, m)
ny, nx = f["temp"].shape
x = (np.arange(nx) + .5) * Lx / nx; y = (np.arange(ny) + .5) * YMAX / ny
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
for a, k, cm in [(ax[0], "HeatRelease", "magma"), (ax[1], "Y(H2)", "viridis")]:
    im = a.pcolormesh(x, y, f[k], cmap=cm, shading="auto"); fig.colorbar(im, ax=a, label=k)
    a.contour(x, y, f["temp"], [m["T_track"]], colors="c", linewidths=0.8)
    for g in ds.index.grids:
        if g.Level == 0: continue
        (x0, y0), (x1, y1) = g.LeftEdge.d[:2] / lF, g.RightEdge.d[:2] / lF
        a.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, lw=0.4,
                                  ec="w" if g.Level == 2 else "lime"))
    a.set_xlim(0, 30); a.set_ylim(20, 50); a.set_aspect("equal")
    a.set_xlabel("x / l_F"); a.set_ylabel("y / l_F")
    a.set_title(f"{m['name'][-5:]} zoom, {k} ; boîtes niveau 1 (vert), niveau 2 (blanc)", fontsize=9)
fig.tight_layout(); fig.savefig("figs/hoq_unstable_zoom_amr.png", dpi=110)

# fine resolution of flame thickness: cells across the HRR peak along columns
hr = f["HeatRelease"]; dy = YMAX / ny
w = [(hr[:, i] > 0.5 * hr[:, i].max()).sum() for i in range(0, nx, 16)]
print("cells per HRR FWHM (fine level): median", np.median(w), "min", np.min(w), "dy/lF", dy)
