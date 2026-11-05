"""Local mixture along the wrinkled front of the saturated B0 flame (two-fuel effect).

  python local_mixture.py RUN_DIR      (uses RUN_DIR/species/from_chk*/plt_sp*, else RUN_DIR/plt* at 20-38 tau_F)

- phi_loc = n_H / (2 n_O) from elemental moles (conserved by chemistry; 1 at stoichiometry for H2 and NH3),
  sampled on the isotherm T_track (reaction zone side, as in Schneider) and on T_pre = T_u + 0.1 (T_b - T_u).
- gamma_loc = 2 X_H2 / (3 X_NH3 + 2 X_H2) on T_pre (preheat zone: fuels separated by preferential diffusion,
  not yet consumed). Stoichiometric cracking conserves the fuel H/N ratio, so gamma has no elemental definition.
- curvature kappa of the isoline (positive = convex towards the unburnt gas), in 1/l_F.
Outputs: local_mixture.csv (all samples), local_mixture_stats.txt, figs/local_mixture_<name>.png
"""
import glob
import json
import os
import sys

import cantera as ct
import contourpy
import numpy as np

MECH_ROOT = os.path.join(os.environ.get("PELE_HOME", "PeleLMeX"), "Submodules/PelePhysics/Mechanisms")


def isoline_samples(T, fields, level, dx):
    """Points of the T = level isoline with bilinear-interpolated fields and the local curvature."""
    lines = contourpy.contour_generator(z=T.T, line_type="Separate").lines(level)
    out = []
    for l in lines:
        if len(l) < 7:
            continue
        x, y = l[:, 0], l[:, 1]              # in grid units (index space of T.T: x = column, y = row)
        # curvature from a smoothed parametrisation (5-point window)
        dxs, dys = np.gradient(x), np.gradient(y)
        ddx, ddy = np.gradient(dxs), np.gradient(dys)
        kap = (dxs * ddy - dys * ddx) / np.maximum((dxs**2 + dys**2) ** 1.5, 1e-12) / dx
        # orient: positive when convex towards unburnt (low T): check T gradient direction
        i0 = np.clip(np.floor(x).astype(int), 0, T.shape[0] - 2)
        j0 = np.clip(np.floor(y).astype(int), 0, T.shape[1] - 2)
        fx, fy = x - i0, y - j0

        def interp(F):
            return ((1 - fx) * (1 - fy) * F[i0, j0] + fx * (1 - fy) * F[i0 + 1, j0]
                    + (1 - fx) * fy * F[i0, j0 + 1] + fx * fy * F[i0 + 1, j0 + 1])
        gTx = interp(np.gradient(T, axis=0))
        gTy = interp(np.gradient(T, axis=1))
        # normal pointing to burnt = grad T; tangent (dxs, dys); sign so that kappa > 0 bulges into unburnt
        sign = np.sign(dxs * gTy - dys * gTx)
        kap = kap * sign
        row = {"x": x * dx, "y": y * dx, "kappa": kap}
        for k, F in fields.items():
            row[k] = interp(F)
        out.append(row)
    return {k: np.concatenate([r[k] for r in out]) for k in out[0]} if out else {}


def main(run):
    import yt
    yt.set_log_level(50)
    m = json.load(open(os.path.join(run, "case.json")))
    gas = ct.Solution(os.path.join(MECH_ROOT, m["mech"], "mechanism.yaml"))
    lF, Tu, Tb = m["l_F"], m["Tu"], m["Tb_eq"]
    T_pre = Tu + 0.1 * (Tb - Tu)
    W = gas.molecular_weights
    nH = np.array([gas.n_atoms(s, "H") for s in gas.species_names]) / W
    nO = np.array([gas.n_atoms(s, "O") for s in gas.species_names]) / W
    iH2, iNH3 = gas.species_index("H2"), gas.species_index("NH3")
    phi0 = 0.5
    g0 = m["case"]["gamma"]
    allrows, lines = [], []
    plts = sorted(glob.glob(os.path.join(run, "species", "from_chk*", "plt_sp*")))
    if not plts:                       # series 2: plotfiles already carry Y(...); states at 20/26/32/38 tau_F
        for p in sorted(glob.glob(os.path.join(run, "plt?????"))):
            L = open(os.path.join(p, "Header")).read().split("\n")
            t = float(L[2 + int(L[1]) + 1]) / m["tau_F"]
            if any(abs(t - tt) < 0.1 for tt in (20, 26, 32, 38)):
                plts.append(p)
    for p in plts:
        ds = yt.load(p)
        lev = ds.index.max_level
        n = (ds.domain_dimensions * ds.refine_by**lev).astype(int)
        cg = ds.covering_grid(lev, left_edge=ds.domain_left_edge, dims=n)
        dx = (ds.domain_width / n).d[0]
        T = cg["boxlib", "temp"].d[:, :, 0]
        Y = np.stack([cg["boxlib", f"Y({s})"].d[:, :, 0] for s in gas.species_names], axis=-1)
        H = Y @ nH
        O = Y @ nO
        phi = H / (2 * O)
        X = Y / W
        X /= X.sum(axis=-1, keepdims=True)
        gam = 2 * X[..., iH2] / np.maximum(3 * X[..., iNH3] + 2 * X[..., iH2], 1e-30)
        hrr = cg["boxlib", "HeatRelease"].d[:, :, 0]
        t = float(ds.current_time) / m["tau_F"]
        for tag, lvl in (("track", m["T_track"]), ("pre", T_pre)):
            s = isoline_samples(T, {"phi": phi, "gamma": gam, "hrr": hrr}, lvl, dx)
            if not s:
                continue
            s["kappa_lF"] = s["kappa"] * lF
            for i in range(len(s["x"])):
                allrows.append((t, tag, s["x"][i] / lF, s["y"][i] / lF, s["kappa_lF"][i], s["phi"][i],
                                s["gamma"][i], s["hrr"][i]))
            ph, ga, ka = s["phi"], s["gamma"], s["kappa_lF"]
            pos, neg = ka > 0.2, ka < -0.2
            lines.append(f"t = {t:5.1f} tau_F, {tag:5s} (T = {lvl:6.1f} K), {len(ph):5d} pts | "
                         f"phi_loc mean {ph.mean():.3f} std {ph.std():.3f} [p5 {np.percentile(ph, 5):.3f}, "
                         f"p95 {np.percentile(ph, 95):.3f}] | gamma_loc mean {ga.mean():.3f} std {ga.std():.3f} "
                         f"[p5 {np.percentile(ga, 5):.3f}, p95 {np.percentile(ga, 95):.3f}] | "
                         f"convex (k>0.2): phi {ph[pos].mean():.3f} gamma {ga[pos].mean():.3f} ; "
                         f"concave (k<-0.2): phi {ph[neg].mean():.3f} gamma {ga[neg].mean():.3f} | "
                         f"corr(phi,kappa) {np.corrcoef(ph, ka)[0, 1]:+.2f} corr(gamma,kappa) {np.corrcoef(ga, ka)[0, 1]:+.2f}")
            print(lines[-1], flush=True)
    hdr = f"# unburnt: phi = {phi0}, gamma = {g0}; T_track = {m['T_track']:.1f} K, T_pre = {T_pre:.1f} K\n"
    open(os.path.join(run, "local_mixture_stats.txt"), "w").write(hdr + "\n".join(lines) + "\n")
    np.savetxt(os.path.join(run, "local_mixture.csv"),
               np.array([r[:1] + (1.0 if r[1] == "track" else 0.0,) + r[2:] for r in allrows]), delimiter=",",
               header="t_tauF,on_track,x_lF,y_lF,kappa_lF,phi_loc,gamma_loc,hrr", comments="", fmt="%.6g")
    # figure: joint PDF on T_pre + phi vs kappa on T_track
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pre = np.array([r for r in allrows if r[1] == "pre"], dtype=object)
    trk = np.array([r for r in allrows if r[1] == "track"], dtype=object)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    h = ax[0].hist2d(pre[:, 5].astype(float), pre[:, 6].astype(float), bins=60, cmap="Blues", density=True)
    ax[0].plot(phi0, g0, "r+", ms=14, mew=2, label="unburnt mixture")
    ax[0].set_xlabel(r"$\varphi_{loc}$ (preheat isotherm)")
    ax[0].set_ylabel(r"$\gamma_{loc}$ (preheat isotherm)")
    ax[0].legend(fontsize=8, frameon=False)
    ax[1].scatter(trk[:, 4].astype(float), trk[:, 5].astype(float), s=1, alpha=0.2, color="#1f4fbf")
    ax[1].axhline(phi0, color="r", lw=1)
    ax[1].set_xlabel(r"$\kappa\, l_F$ (T_track isoline, > 0: convex to unburnt)")
    ax[1].set_ylabel(r"$\varphi_{loc}$ (reaction zone)")
    ax[1].set_xlim(-3, 3)
    ax[2].scatter(pre[:, 4].astype(float), pre[:, 6].astype(float), s=1, alpha=0.2, color="#c0392b")
    ax[2].axhline(g0, color="r", lw=1)
    ax[2].set_xlabel(r"$\kappa\, l_F$ (preheat isotherm)")
    ax[2].set_ylabel(r"$\gamma_{loc}$")
    ax[2].set_xlim(-3, 3)
    for a in ax:
        a.grid(alpha=0.25, lw=0.5)
    fig.suptitle(f"{m['name']}: local mixture along the saturated front (4 snapshots, 20–38 τ_F)", fontsize=10)
    fig.tight_layout()
    os.makedirs(os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs"), exist_ok=True)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figs", f"local_mixture_{m['name']}.png")
    fig.savefig(out, dpi=130)
    print(out)


if __name__ == "__main__":
    main(sys.argv[1])
