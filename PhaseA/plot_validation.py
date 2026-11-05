"""Validation figure: dispersion relations of runs/val_L24_* against Lehmann et al. 2024, Fig. 3a
(X_H2,F = 0.5, phi = 0.5, T_u = 500 K, 1 and 10 bar), plus the pilot (Chabot conditions).

Lehmann's crosses were digitised from the PDF (rendered at 300 dpi, markers found by colour, see
cases/lehmann2024_fig3a_digitised.csv); accuracy ~ +/- 0.01 in k_bar and omega_bar.

  python plot_validation.py   -> figs/validation_lehmann_fig3a.png
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = {"1 bar": ("runs/val_L24_XH0.5_phi0.5_T500_P1", "#1f4fbf"),
        "10 bar": ("runs/val_L24_XH0.5_phi0.5_T500_P10", "#c0392b")}
PILOT = "runs/A_g0.43_phi0.50_P10.0_T750"


def kc_of(d):
    """Cut-off wavenumber: linear interpolation at the first sign change + -> - (all modes: near
    omega = 0 the amplitude barely changes, so r2 is low by construction)."""
    s = d.sort_values("k_bar")
    w, k = s.omega_bar.values, s.k_bar.values
    for i in range(len(w) - 1):
        if w[i] > 0 >= w[i + 1]:
            return k[i] + w[i] / (w[i] - w[i + 1]) * (k[i + 1] - k[i])
    return np.nan


def lehmann_scales(run, p_bar):
    """(k, omega) factors from our normalisation to Lehmann's (Table 1 of the supplement)."""
    m = json.load(open(os.path.join(HERE, run, "case.json")))
    L = pd.read_csv(os.path.join(HERE, "cases", "lehmann2024_table1.csv"))
    r = L[(L.Tu == 500) & (L.p == p_bar) & (L.phi == 0.5) & (L.XH2F == 0.5)].iloc[0]
    tauL = r.lF_um * 1e-6 / (r.sL_cms / 100)
    return r.lF_um * 1e-6 / m["l_F"], tauL / m["tau_F"], r


def load(run):
    f = os.path.join(HERE, run, "dispersion.csv")
    return pd.read_csv(f) if os.path.exists(f) else None


def main():
    ref = pd.read_csv(os.path.join(HERE, "cases", "lehmann2024_fig3a_digitised.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    for lab, (run, c) in RUNS.items():
        r = ref[ref.p_bar == float(lab.split()[0])]
        ax.plot(r.k_bar, r.omega_bar, "x", color=c, ms=7, mew=2, label=f"Lehmann 2024 (Zhang), {lab}")
        d = load(run)
        if d is None:
            continue
        ok = d.r2 > 0.9
        ax.plot(d.k_bar[ok], d.omega_bar[ok], "-o", color=c, lw=2, ms=4, label=f"This work (Mei), {lab}")
        ax.plot(d.k_bar[~ok], d.omega_bar[~ok], "o", color=c, mfc="none", ms=4)
        ax.plot(d.k_bar, d.omega_DL_bar, ":", color=c, lw=1.5)
        a, b, r = lehmann_scales(run, float(lab.split()[0]))
        ax.plot(d.k_bar[ok] * a, d.omega_bar[ok] * b, "--", color=c, lw=1.5,
                label=f"This work rescaled to Lehmann's $l_F$, $s_L$, {lab}")
        print(f"{lab}: ours omax {d.omega_bar[ok].max():.3f}, kc {kc_of(d):.2f} | in Lehmann scales "
              f"omax {d.omega_bar[ok].max()*b:.3f}, kc {kc_of(d)*a:.2f} | Lehmann Table 1 omax {r.omax}, kc {r.kc} "
              f"| digitised fig. omax {ref[ref.p_bar == float(lab.split()[0])].omega_bar.max():.3f}")
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlim(0, 2)
    ax.set_ylim(-1.2, 1.5)
    ax.set_xlabel(r"$\bar k = k\,l_F$")
    ax.set_ylabel(r"$\bar\omega = \omega\,\tau_F$")
    ax.set_title(r"$X_{H_2,F}$ = 0.5, φ = 0.5, $T_u$ = 500 K  (dotted: DL)", fontsize=10)
    ax.legend(fontsize=7, frameon=False, loc="lower left")
    ax.grid(alpha=0.25, lw=0.5)

    ax = axes[1]
    d = load(PILOT)
    if d is not None:
        ok = d.r2 > 0.9
        ax.plot(d.k_bar[ok], d.omega_bar[ok], "-o", color="#1f4fbf", lw=2, ms=4, label="This work (Mei)")
        ax.plot(d.k_bar[~ok], d.omega_bar[~ok], "o", color="#1f4fbf", mfc="none", ms=4, label="fit r² < 0.9")
        ax.plot(d.k_bar, d.omega_DL_bar, ":", color="#1f4fbf", lw=1.5, label="Darrieus–Landau")
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlim(0, 2)
    ax.set_ylim(-1.2, 1.5)
    ax.set_xlabel(r"$\bar k = k\,l_F$")
    ax.set_title("Pilot (Chabot): γ = 0.43, φ = 0.5, 10 atm, 750 K", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figs"), exist_ok=True)
    out = os.path.join(HERE, "figs", "validation_lehmann_fig3a.png")
    fig.savefig(out, dpi=140)
    print(out)


if __name__ == "__main__":
    main()
