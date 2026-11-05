# Reproducibility package: unstable partially cracked ammonia flames at walls

Code, inputs and post-processed data for the paper

> L. Alzaben and R. Boukharfane, *Pulsed wall heat loads, delayed ammonia burnout and cold-wall nitrous oxide in unstable
> partially cracked ammonia flames*.

All simulations were run with PeleLMeX (low-Mach number, detailed kinetics, mixture-averaged transport with Soret
diffusion) on the Toubkal cluster (UM6P). This package lets you

1. **rebuild every figure and table of the paper that does not show a raw field** from the post-processed data
   included here (about 20 minutes on one core);
2. **re-run the post-processing and analyses** from the same data;
3. **re-run any simulation** from its exact input files, solver versions, problem set-ups and mechanisms.

Raw solver output (plotfiles and checkpoints, several TB) is not included; see *What is not included*.

## Layout

The package mirrors the working directory of the study, so that each script runs unchanged from its folder.

| folder | content |
|---|---|
| `A0/` | one-dimensional unstretched flames with Cantera (`a0lib.py`, `run_a0.py`): $s_L$, $l_F$, $l_D$, Ze, Le_eff, profiles; results in `A0/results/` |
| `PhaseA/` | linear stability (dispersion relations): case generator `make_case.py`, case lists `cases/`, front extraction `extract_front.py`, restart helper; one folder per case in `runs/` |
| `PhaseB0/` | saturated unstable free flames: `make_b0.py`, `post_b0.py` ($I_0$, $s_c$, $A$), `local_mixture.py` ($\phi_{loc}$, $\gamma_{loc}$, curvature) |
| `PhaseB/` | head-on quenching (unstable, planar, cold walls, crest mixtures, mechanism and grid checks): `make_hoq.py`, `make_crest_hoq.py`, `make_mech_hoq.py`, `post_hoq.py`, `emissions_tq.py`, `wall_statistics.py`, `crest_at_quench.py`, `coldwall_analysis.py`, `summary_tw.py`, `collapse_tw.py`, `no_budget.py`, `paper_figs_coldwall.py` |
| `analysis/` | analyses of the paper and of the review: stability database, Matalon theory, local events, local-mixture model, crest/trough and equal-local-age emissions, wall energy, $I_0$ marker, $\mathcal H$ sensitivity, sanity checks |
| `paper/` | `main.tex`, `refs.bib`, `make_figs.py` (all figures except the cold-wall ones) |
| `solver/` | solver commits (`VERSIONS.md`), local PelePhysics patch, CEPTR-generated mechanisms |
| `mechanisms/` | Mei et al. (2021) mechanism, its 31-species reduction (`reduce_for_pele.py`), Zhang et al. (2021) |
| `environment/` | Python packages (`requirements.txt`), cluster modules (`MODULES.md`) |
| `data_manifest.tsv` | list of the post-processed data files with size and SHA-256 |

Each run folder (`PhaseA/runs/<case>`, `PhaseB0/runs/<case>`, `PhaseB/runs/<case>`) contains its exact inputs
(`case.json`, `input.2d`, `run.sbatch`, `pmf.dat`) and its post-processed output. Case names encode the conditions,
e.g. `A_g0.43_phi0.50_P10.0_T750` = cracking ratio 0.43, $\phi = 0.5$, 10 atm, $T_u = 750$ K; HOQ names add the
initial state (`_t26.0` = saturated state at $26\,\tau_F$), `_planar`, the wall temperature (`_Tw500`), the initial
distance (`_d60`), `crest_..._q50` (planar flame at a crest mixture), `_mech_Zhang2021NH3` or `_L3` (grid check).

## Quick start: rebuild the figures

```bash
python3.12 -m venv .venv && . .venv/bin/activate
pip install -r environment/requirements.txt
export PELE_HOME=/path/to/PeleLMeX      # only needed for the Cantera steps that read the PelePhysics mechanisms
cd paper && python make_figs.py dispersion stability I0 mixture flux emissions regime local_quench crest_trough local_model
cd ../PhaseB && python paper_figs_coldwall.py
cd ../paper && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

`PELE_HOME` points to a PeleLMeX checkout in which `solver/PelePhysics_Mechanisms/*` have been copied to
`Submodules/PelePhysics/Mechanisms/` (the scripts load `mechanism.yaml` from there).

## Figures and tables

| paper item | produced by | from |
|---|---|---|
| dispersion relations | `paper/make_figs.py dispersion` | `PhaseA/runs/A_*_P10.0_T750/dispersion.csv` |
| stability map (database, 209 cases) | `analysis/stability_database.py`, `make_figs.py stability` | `PhaseA/runs/*/dispersion.csv`, `A0/results/mei2021/cases/*.json` |
| table of saturated flames, $I_0$ figure | `make_figs.py I0` | `PhaseB0/runs/*/b0_timeseries.csv`, `A0` cases |
| local mixture | `make_figs.py mixture` | `PhaseB0/runs/*/local_mixture.csv` |
| flux table, mean/maximum flux figure | `PhaseB/wall_statistics.py`, `make_figs.py flux` | `PhaseB/runs/*/hoq_wallflux.npz`, `hoq_timeseries.csv` |
| local-mixture reconstruction | `analysis/local_events.py 8 2.0`, `analysis/item7_local_model.py _lead2`, `make_figs.py local_model` | HOQ files above (+ plotfiles for the crest tips) |
| local quenching events, cumulative heat | `analysis/local_events.py`, `analysis/wall_energy.py`, `make_figs.py local_quench` | as above |
| near-wall species | `PhaseB/emissions_tq.py`, `make_figs.py emissions` | `PhaseB/emissions_tq.csv` |
| early/late quenching, NO budget | `analysis/conditioned_emissions.py`, `analysis/local_age_emissions.py`, `make_figs.py crest_trough` | `analysis/*.csv` |
| cold-wall collapse and NO budget | `PhaseB/summary_tw.py`, `collapse_tw.py`, `no_budget.py`, `paper_figs_coldwall.py` | `PhaseB/runs/*/coldwall.npz`, `PhaseB/*.csv` |
| regime map | `make_figs.py regime` | `PhaseB/tw_collapse.csv`, `emissions_tq.csv` |
| appendix values | `analysis/i0_marker.py`, `h_sensitivity.py`, `matalon.py`, `wave2_sanity.py` | `analysis/*.csv` |

Scripts marked as reading plotfiles or HOQ snapshots need the raw output; their results are included as CSV.

## Re-running the simulations

1. **Solver.** Check out the commits of `solver/VERSIONS.md`, apply `solver/pelephysics_transport.patch` to
   PelePhysics, copy `solver/PelePhysics_Mechanisms/*` into `Submodules/PelePhysics/Mechanisms/`, and build with the
   `GNUmakefile` and `pelelmex_prob.*` of `PhaseA/template`, `PhaseB0/build/10_000atm` and `PhaseB/build/*`
   (`make -j PDIR=10_000atm [MECH=Zhang2021NH3]`). Environment: `environment/MODULES.md` (toolchain foss/2024a: GCC 13.3, OpenMPI 5).
2. **1D flames.** `cd A0 && python run_a0.py` (Cantera, same kinetics and transport, with Soret).
3. **Linear stability.** `cd PhaseA && python make_case.py <case> [--cfl 0.3]` writes the run folder; submit
   `run.sbatch`; `extract_front.py <run>` gives `front.npz`, `dispersion.csv`, `quality.json`.
4. **Saturated flames.** `cd PhaseB0 && python make_b0.py ...`, then `post_b0.py`, `local_mixture.py`.
5. **Head-on quenching.** `cd PhaseB && python make_hoq.py unstable|planar --b0 ../PhaseB0/runs/<B0> [--tw 500] [--dist 60] [--levels 3]`;
   crest mixtures with `make_crest_hoq.py`, Zhang mechanism with `make_mech_hoq.py`; then `post_hoq.py`, `emissions_tq.py`,
   `wall_statistics.py`, `crest_at_quench.py`, `coldwall_analysis.py`.
6. **Analyses.** Scripts of `analysis/`, then figures as above.

The `run.sbatch` files are kept as run records: account, QoS and node exclusions are specific to Toubkal, and the
project path is written `${FQ_ROOT}`.

Approximate cost on Toubkal (56 cores per node): linear stability $2.2\times10^5$ core-hours (212 cases,
2 nodes, 6--11 h each), saturated flames $2.8\times10^4$, head-on quenching $2.4\times10^4$.
