"""Planar HOQ at the local mixture of the crests of the unstable front (Schneider Part I argument).

  python make_crest_hoq.py --b0 ../PhaseB0/runs/B0_g0.43_phi0.50_P10.0_T750 --phi_loc 0.436 --gamma_loc 0.263 --tag cmean

The crests of the saturated unstable flame burn a locally enriched mixture (phi_loc, gamma_loc on the T_track isotherm,
PhaseB0/local_mixture.py). The fresh mixture (phi_u, gamma_u) whose unstretched 1D flame shows the same (phi_loc,
gamma_loc) on the same isotherm (absolute temperature T_track of the reference B0) is found by Newton iterations on
Cantera free flames (Mei mechanism, Soret, as Phase A). For gamma = 1 only phi_u is solved for.
A pseudo-B0 directory (pmf.dat, case.json, input.2d with the keys make_hoq.py reads) is written in crest/, then
make_hoq.py planar builds runs/HOQ_crest_<...>_planar. Its wall heat flux, in W/m2, is compared with the unstable HOQ.
"""
import argparse
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "PhaseA"))
import make_case as mc  # noqa: E402


def local_mix(gas, f, T_track):
    W = gas.molecular_weights
    nH = np.array([gas.n_atoms(s, "H") for s in gas.species_names]) / W
    nO = np.array([gas.n_atoms(s, "O") for s in gas.species_names]) / W
    k = np.argmax(f.T > T_track)
    w = (T_track - f.T[k - 1]) / (f.T[k] - f.T[k - 1])
    Y = (1 - w) * f.Y[:, k - 1] + w * f.Y[:, k]
    X = (1 - w) * f.X[:, k - 1] + w * f.X[:, k]
    phi = (Y @ nH) / (2 * (Y @ nO))
    xh, xn = X[gas.species_index("H2")], X[gas.species_index("NH3")]
    return phi, 2 * xh / (3 * xn + 2 * xh)


def solve(case):
    return mc.flame(case)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--b0", required=True)
    ap.add_argument("--phi_loc", type=float, required=True)
    ap.add_argument("--gamma_loc", type=float, default=None)
    ap.add_argument("--tag", required=True)
    a = ap.parse_args()
    m0 = json.load(open(os.path.join(a.b0, "case.json")))
    c0 = dict(m0["case"])
    T_track = m0["T_track"]
    g_fixed = c0["gamma"] >= 1.0 or a.gamma_loc is None
    x = np.array([c0["phi"], c0["gamma"]])
    tgt = np.array([a.phi_loc, a.gamma_loc if a.gamma_loc is not None else np.nan])

    def F(v):
        c = dict(c0, phi=float(v[0]), gamma=1.0 if g_fixed else float(v[1]))
        gas, f, *_ = solve(c)
        return np.array(local_mix(gas, f, T_track))

    for it in range(8):
        r = F(x) - tgt
        res = abs(r[0]) if g_fixed else np.abs(r).max()
        print(f"iter {it}: phi_u = {x[0]:.4f} gamma_u = {x[1]:.4f} -> local {r + tgt} (target {tgt})", flush=True)
        if res < 2e-4:
            break
        h = 0.005
        if g_fixed:
            d = (F(x + [h, 0])[0] - (r[0] + tgt[0])) / h
            x[0] -= r[0] / d
        else:
            J = np.empty((2, 2))
            for j in range(2):
                e = np.zeros(2); e[j] = h
                J[:, j] = (F(x + e) - (r + tgt)) / h
            x -= np.linalg.solve(J, r)
    case = dict(c0, phi=float(x[0]), gamma=1.0 if g_fixed else float(x[1]))
    gas, f, P, pkey, Tb, lF, x_fl = solve(case)
    sL, Tu = f.velocity[0], f.T[0]
    name = f"B0_crest_{os.path.basename(a.b0)[3:]}_{a.tag}"
    out = os.path.join(HERE, "crest", name)
    os.makedirs(out, exist_ok=True)
    mc.write_pmf(os.path.join(out, "pmf.dat"), gas, f)
    Lx = 100 * lF
    meta = {"name": name, "case": case, "P_Pa": P, "pdir": mc.PDIRS[pkey], "sL": sL, "Tu": Tu, "Tb_eq": Tb,
            "l_F": lF, "tau_F": lF / sL, "T_track": Tu + mc.C_TRACK * (Tb - Tu), "Lx": Lx, "Ly": Lx,
            "dx_fine": Lx / 512 / 4, "mech": mc.MECH, "sigma_1d": float(f.density[0] / f.density[-1]),
            "rho_u": float(f.density[0]), "Y_u": dict(zip(gas.species_names, map(float, f.Y[:, 0]))),
            "crest_of": os.path.abspath(a.b0), "phi_loc_target": a.phi_loc, "gamma_loc_target": a.gamma_loc,
            "local_mix_1d": list(local_mix(gas, f, T_track))}
    json.dump(meta, open(os.path.join(out, "case.json"), "w"), indent=1, default=float)
    open(os.path.join(out, "input.2d"), "w").write(
        f"amr.gradT.adjacent_difference_greater = {0.02 * (Tb - Tu):.1f}\nprob.standoff = {30.0 * lF - x_fl:.10e}\n")
    print(f"{name}: phi_u = {case['phi']:.4f}, gamma_u = {case['gamma']:.4f}, sL = {sL:.4f} m/s "
          f"(ref {m0['sL']:.4f}, ratio {sL / m0['sL']:.3f}), l_F = {lF * 1e6:.1f} um", flush=True)
    subprocess.run([sys.executable, os.path.join(HERE, "make_hoq.py"), "planar", "--b0", out], check=True)


if __name__ == "__main__":
    main()
