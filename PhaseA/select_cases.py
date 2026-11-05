"""Phase A case selection from the A0 results (Mei 2021, MA + Soret).

Rules (proposal.md section 4.4/5, updated with A0):
  - T_u in {500, 750} K, phi in {0.4, 0.5, 0.6, 0.7, 0.8}, P in {1, 5, 10} atm, all gamma;
  - kept if the A0 flame converged and s_L >= SL_MIN (otherwise at the flammability limit);
  - wave 1 = phi in {0.4, 0.5, 0.6} and P in {1, 10} atm (the cases B0/B will be chosen from),
    wave 2 = the rest;
  - plus the Lehmann 2024 / D'Alessio validation cases (no N2, see VALIDATION).

Per case the Lehmann setup is derived: domain 80 l_F x 20 l_F, base dx = l_F/12.8,
2 static levels in y in [7, 13] l_F (dx_fine = l_F/51.2), front at y0 = 10 l_F,
duration N_TAU tau_F. Cost is a rough estimate to be replaced by the scaling test.

  python select_cases.py
"""
import os

import numpy as np
import pandas as pd

A0 = "../A0/results/mei2021/a0_summary.csv"
SL_MIN = 0.03          # m/s
T_U = [500.0, 750.0]
PHI = [0.4, 0.5, 0.6, 0.7, 0.8]
WAVE1 = {"phi": [0.4, 0.5, 0.6], "P_atm": [1.0, 10.0]}
N_TAU = 8              # tau_F simulated (Lehmann: plateau of omega within ~6 tau_F)
CFL = 0.5
US_PER_CELL_STEP = 100e-6   # core-s per cell and step (guess for ~31 species, to be measured)
C_TRACK = 0.4          # front isotherm T = T_u + C_TRACK (T_b - T_u); Lehmann's 1000 K at 298 K ~ 0.38-0.45

# Lehmann 2024, Fig. 3a (with D'Alessio 2024): NH3/H2 blend without N2
VALIDATION = [{"name": "val_L24_XH0.5_phi0.5_T500_P1", "X_H2_F": 0.5, "phi": 0.5, "T_u": 500.0, "P_bar": 1.0},
              {"name": "val_L24_XH0.5_phi0.5_T500_P10", "X_H2_F": 0.5, "phi": 0.5, "T_u": 500.0, "P_bar": 10.0}]


def lehmann_setup(row):
    lF = row.delta_T
    nx_base, ny_base = 1024, 256                     # 80 l_F / (l_F/12.8), 20 l_F / (l_F/12.8)
    n_fine = int(80 * 51.2) * int(6 * 51.2)          # level 2 band, 6 l_F high
    n_l1 = int(80 * 25.6) * int(8 * 25.6)            # level 1 band, ~8 l_F high
    cells = nx_base * ny_base + n_l1 + n_fine
    steps = N_TAU * 51.2 * row.sigma / CFL           # dt = CFL dx_fine / (sigma s_L), no subcycling
    return pd.Series({
        "l_F_um": lF * 1e6, "Lx_mm": 80 * lF * 1e3, "Ly_mm": 20 * lF * 1e3,
        "dx_fine_um": lF / 51.2 * 1e6, "T_track": row.T_u + C_TRACK * (row.Tb - row.T_u),
        "t_end_ms": N_TAU * row.tau_F * 1e3, "cells_M": cells / 1e6, "steps": steps,
        "core_h_est": cells * steps * US_PER_CELL_STEP / 3600})


def main():
    df = pd.read_csv(A0)
    df["ZeLe"] = df.Ze * (df["Le_eff.u.diff"] - 1)
    sel = df[df.T_u.isin(T_U) & df.phi.isin(PHI)].copy()
    n_all = len(sel)
    drop = sel[(sel.status != "ok") | (sel.sL < SL_MIN)]
    sel = sel[(sel.status == "ok") & (sel.sL >= SL_MIN)].copy()
    sel["wave"] = np.where(sel.phi.isin(WAVE1["phi"]) & sel.P_atm.isin(WAVE1["P_atm"]), 1, 2)
    # the most unstable gamma of each (T_u, phi, P) by the A0 parameter Ze (Le_eff - 1)
    imin = sel.groupby(["T_u", "phi", "P_atm"]).ZeLe.idxmin()
    sel["ZeLe_min_of_row"] = sel.index.isin(imin)
    sel = pd.concat([sel, sel.apply(lehmann_setup, axis=1)], axis=1)
    sel["name"] = [f"A_g{g:.2f}_phi{p:.2f}_P{P:04.1f}_T{T:03.0f}" for g, p, P, T in
                   zip(sel.gamma, sel.phi, sel.P_atm, sel.T_u)]
    cols = ["name", "wave", "gamma", "X_H2_F", "phi", "P_atm", "T_u", "sL", "Tb", "delta_T", "tau_F",
            "sigma", "Ze", "Le_eff.u.diff", "ZeLe", "ZeLe_min_of_row", "I0_schneider",
            "l_F_um", "Lx_mm", "Ly_mm", "dx_fine_um", "T_track", "t_end_ms", "cells_M", "steps", "core_h_est"]
    sel = sel.sort_values(["wave", "T_u", "phi", "P_atm", "gamma"])[cols]
    os.makedirs("cases", exist_ok=True)
    sel.to_csv("cases/phaseA_cases.csv", index=False)
    pd.DataFrame(VALIDATION).to_csv("cases/phaseA_validation.csv", index=False)
    drop[["gamma", "phi", "P_atm", "T_u", "status", "sL"]].sort_values(["T_u", "phi", "P_atm", "gamma"]) \
        .to_csv("cases/phaseA_excluded.csv", index=False)

    print(f"candidates {n_all}, excluded {len(drop)} (not converged or s_L < {SL_MIN*100:g} cm/s), kept {len(sel)}")
    print(sel.groupby("wave").agg(n=("name", "size"), core_h=("core_h_est", "sum"),
                                  steps_max=("steps", "max")).round(0).to_string())
    print("\nexcluded per (T_u, phi):")
    print(drop.groupby(["T_u", "phi"]).size().unstack(fill_value=0).to_string())
    print("\ncases per (T_u, phi, P) (gamma values kept):")
    print(sel.groupby(["T_u", "phi", "P_atm"]).gamma.apply(lambda g: " ".join(f"{x:g}" for x in g)).to_string())
    print("\nmost unstable gamma (A0, Ze(Le_eff-1)) per row:")
    m = sel[sel.ZeLe_min_of_row]
    print(m.pivot_table(index=["T_u", "phi"], columns="P_atm", values="gamma").to_string())
    print(f"\nl_F range {sel.l_F_um.min():.0f}-{sel.l_F_um.max():.0f} um, "
          f"t_end range {sel.t_end_ms.min():.2f}-{sel.t_end_ms.max():.1f} ms")


if __name__ == "__main__":
    main()
