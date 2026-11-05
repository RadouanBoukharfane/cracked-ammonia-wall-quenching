"""Mechanism test: dispersion relations with Mei and Zhang (Lehmann's mechanism) against Lehmann et al. 2024 Fig. 3a
(X_H2,F = 0.5, phi = 0.5, T_u = 500 K). Dashed: Zhang runs expressed in Lehmann's l_F, s_L (his Table S1).
  python plot_validation_zhang.py -> figs/validation_zhang.png"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from plot_validation import HERE, kc_of, lehmann_scales, load

ref = pd.read_csv(os.path.join(HERE, "cases", "lehmann2024_fig3a_digitised.csv"))
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, (P, tag) in zip(axes, ((1, "P1"), (10, "P10"))):
    r = ref[ref.p_bar == P]
    ax.plot(r.k_bar, r.omega_bar, "x", color="k", ms=7, mew=2, label="Lehmann 2024 (Zhang, FlameMaster scales)")
    for mech, suf, c in (("Mei", "", "#1f4fbf"), ("Zhang", "_zhang", "#c0392b")):
        run = f"runs/val_L24_XH0.5_phi0.5_T500_{tag}{suf}"
        d = load(run)
        ok = d.r2 > 0.9
        ax.plot(d.k_bar[ok], d.omega_bar[ok], "-o", color=c, lw=2, ms=3.5, label=f"This work, {mech}")
        if mech == "Zhang":
            a, b, row = lehmann_scales(run, float(P))
            ax.plot(d.k_bar[ok] * a, d.omega_bar[ok] * b, "--", color=c, lw=1.5,
                    label="This work, Zhang, in Lehmann's $l_F$, $s_L$")
            print(f"{P} bar Zhang in Lehmann scales: omax {d.omega_bar[ok].max()*b:.3f} (L: {r.omega_bar.max():.3f}), "
                  f"kc {kc_of(d)*a:.2f} (L table: {row.kc})")
        ax.plot(d.k_bar, d.omega_DL_bar, ":", color=c, lw=1)
    ax.axhline(0, color="0.5", lw=0.8)
    ax.set_xlim(0, 2.2)
    ax.set_ylim(-0.8, 1.4)
    ax.set_xlabel(r"$\bar k = k\,l_F$")
    ax.set_title(f"{P} bar  (dotted: Darrieus–Landau)", fontsize=10)
    ax.grid(alpha=0.25, lw=0.5)
axes[0].set_ylabel(r"$\bar\omega = \omega\,\tau_F$")
axes[0].legend(fontsize=7.5, frameon=False, loc="lower left")
fig.suptitle(r"$X_{H_2,F}$ = 0.5, φ = 0.5, $T_u$ = 500 K: mechanism and normalisation", fontsize=11)
fig.tight_layout()
out = os.path.join(HERE, "figs", "validation_zhang.png")
fig.savefig(out, dpi=140)
print(out)
