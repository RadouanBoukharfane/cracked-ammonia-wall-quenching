"""TODO item 9: planar HOQ with another kinetic mechanism (Zhang et al. 2021, 10 bar CEPTR build) for robustness.

  python make_mech_hoq.py --mech Zhang2021NH3 --pdir 9_869atm --P_atm 9.869 --gamma 0.43 [--tw 750 500]

Writes a pseudo-B0 directory crest/B0_mech_<...>/ (pmf.dat in the species order of the chosen mechanism, case.json,
input.2d keys read by make_hoq.py) from a Cantera free flame of that mechanism (phi = 0.5, T_u = 750 K), then calls
make_hoq.py planar for each wall temperature. Pressure of the Zhang build: 10 bar (9.869 atm), i.e. 1.3 % below the
10 atm of the Mei runs; normalised quantities are compared.
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "PhaseA"))
import make_case as mc  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mech", default="Zhang2021NH3")
    ap.add_argument("--pdir", default="9_869atm")
    ap.add_argument("--P_atm", type=float, default=9.869)
    ap.add_argument("--gamma", type=float, default=0.43)
    ap.add_argument("--phi", type=float, default=0.5)
    ap.add_argument("--tw", type=float, nargs="+", default=[750.0, 500.0])
    a = ap.parse_args()
    mc.MECH, mc.PELE_MECH = a.mech, os.path.join(mc.MECH_ROOT, a.mech)
    case = dict(gamma=a.gamma, phi=a.phi, P_atm=a.P_atm, T_u=750.0)
    gas, f, P, pkey, Tb, lF, x_fl = mc.flame(case)
    sL, Tu = f.velocity[0], f.T[0]
    name = f"B0_mech_{a.mech}_g{a.gamma:.2f}_phi{a.phi:.2f}_P{a.P_atm:.1f}_T750"
    out = os.path.join(HERE, "crest", name)
    os.makedirs(out, exist_ok=True)
    mc.write_pmf(os.path.join(out, "pmf.dat"), gas, f)
    Lx = 100 * lF
    meta = {"name": name, "case": case, "P_Pa": P, "pdir": a.pdir, "sL": sL, "Tu": Tu, "Tb_eq": Tb, "l_F": lF,
            "tau_F": lF / sL, "T_track": Tu + mc.C_TRACK * (Tb - Tu), "Lx": Lx, "Ly": Lx, "dx_fine": Lx / 512 / 4,
            "mech": a.mech, "sigma_1d": float(f.density[0] / f.density[-1]), "rho_u": float(f.density[0]),
            "Y_u": dict(zip(gas.species_names, map(float, f.Y[:, 0])))}
    json.dump(meta, open(os.path.join(out, "case.json"), "w"), indent=1, default=float)
    open(os.path.join(out, "input.2d"), "w").write(
        f"amr.gradT.adjacent_difference_greater = {0.02 * (Tb - Tu):.1f}\nprob.standoff = {30.0 * lF - x_fl:.10e}\n")
    print(f"{name}: sL = {sL:.4f} m/s, l_F = {lF * 1e6:.1f} um, T_ad = {Tb:.0f} K", flush=True)
    for tw in a.tw:
        cmd = [sys.executable, os.path.join(HERE, "make_hoq.py"), "planar", "--b0", out]
        if tw != 750.0:
            cmd += ["--tw", f"{tw:g}"]
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
