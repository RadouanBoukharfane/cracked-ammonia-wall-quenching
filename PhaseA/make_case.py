"""Builds one Phase A case directory (Lehmann 2024 method) from cases/phaseA_cases.csv
or cases/phaseA_validation.csv.

  python make_case.py NAME [--nodes 2] [--hours 6] [--ntau 8] [--outdir runs] [--steps N]

Writes runs/NAME/: pmf.dat (Cantera 1D flame, MA + Soret, the PelePhysics YAML), input.2d,
case.json (scales), run.sbatch. --steps N: stop after N steps (scaling test).
"""
import argparse
import json
import os
import re
import sys

import cantera as ct
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "A0"))
import a0lib  # noqa: E402

MECH_ROOT = os.path.join(os.environ.get("PELE_HOME", "PeleLMeX"), "Submodules/PelePhysics/Mechanisms")
MECH = "Mei2021NH3"               # set by --mech (Zhang2021NH3 for the Lehmann comparison)
PELE_MECH = os.path.join(MECH_ROOT, MECH)
NMODES, A0_LF = 40, 1e-6          # Al Kassar perturbation (Lehmann eq. 17)
SEED = 20241015                   # same phases for every case
LX_LF, LY_LF, Y0_LF = 80.0, 20.0, 10.0
BAND_LF = (7.0, 13.0)             # 2 static levels (Lehmann)
NX, NY = 1024, 256                # base dx = l_F / 12.8
C_TRACK = 0.4
PDIRS = {(1.0, "atm"): "1_000atm", (5.0, "atm"): "5_000atm", (10.0, "atm"): "10_000atm",
         (1.0, "bar"): "0_987atm", (10.0, "bar"): "9_869atm"}


def find_case(name):
    for f in ("cases/phaseA_cases.csv", "cases/phaseA_validation.csv"):
        d = pd.read_csv(os.path.join(HERE, f))
        r = d[d.name == name]
        if len(r):
            return r.iloc[0].to_dict()
    raise SystemExit(f"case {name} not found")


def species_order():
    hdir = sorted(d for d in os.listdir(PELE_MECH) if d.endswith("atm"))[0]
    h = open(os.path.join(PELE_MECH, hdir, "mechanism.H")).read()
    return [s for _, s in sorted((int(i), s) for s, i in re.findall(r"#define (\w+)_ID (\d+)", h))]


def set_state(gas, case):
    if "gamma" in case and not pd.isna(case.get("gamma", np.nan)):
        a0lib.set_mixture(gas, case["gamma"], case["phi"], case["T_u"], case["P_atm"])
        return case["P_atm"] * ct.one_atm, (case["P_atm"], "atm")
    xh = case["X_H2_F"]           # Lehmann blend, no N2
    gas.set_equivalence_ratio(case["phi"], {"H2": xh, "NH3": 1 - xh}, a0lib.AIR)
    gas.TP = case["T_u"], case["P_bar"] * 1e5
    return case["P_bar"] * 1e5, (case["P_bar"], "bar")


def flame(case):
    gas = ct.Solution(os.path.join(PELE_MECH, "mechanism.yaml"))
    P, pkey = set_state(gas, case)
    g_eq = ct.Solution(os.path.join(PELE_MECH, "mechanism.yaml"))
    Tb, _ = a0lib.adiabatic_T(g_eq, gas)
    f = a0lib.solve_free_flame(gas, a0lib.initial_width(P / ct.one_atm, case["T_u"]), soret=True)
    dTdx = np.gradient(f.T, f.grid)
    lF = (Tb - f.T[0]) / np.max(np.abs(dTdx))
    return gas, f, P, pkey, Tb, lF, f.grid[np.argmax(np.abs(dTdx))]


def write_pmf(path, gas, f):
    order = species_order()
    assert sorted(order) == sorted(gas.species_names), "PelePhysics / Cantera species mismatch"
    X = np.array([f.X[gas.species_index(s)] for s in order]).T
    data = np.column_stack([f.grid * 100, f.T, f.velocity * 100, f.density * 1e-3, X])
    with open(path, "w") as h:
        h.write('VARIABLES = "X" "temp" "u" "rho" ' + " ".join(f'"{s}"' for s in order) + "\n")
        h.write(f" ZONE I={len(f.grid)} FORMAT=POINT\n")
        np.savetxt(h, data, fmt="%.16e")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--nodes", type=int, default=2)
    ap.add_argument("--hours", type=float, default=0.0, help="0 = estimated from sigma")
    ap.add_argument("--ntau", type=float, default=8.0)
    ap.add_argument("--steps", type=int, default=0)
    ap.add_argument("--outdir", default="runs")
    ap.add_argument("--tag", default="", help="suffix of the run directory")
    ap.add_argument("--mech", default="Mei2021NH3", help="PelePhysics mechanism directory")
    ap.add_argument("--a0", type=float, default=A0_LF, help="perturbation amplitude per mode, in l_F")
    ap.add_argument("--ode-rtol", type=float, default=1e-6)
    ap.add_argument("--cfl", type=float, default=0.5, help="0.3 removes the chemistry noise of some slow flames")
    ap.add_argument("--ode-atol", type=float, default=1e-5)
    a = ap.parse_args()

    global MECH, PELE_MECH
    MECH, PELE_MECH = a.mech, os.path.join(MECH_ROOT, a.mech)
    case = find_case(a.name)
    gas, f, P, pkey, Tb, lF, x_fl = flame(case)
    sL, Tu = f.velocity[0], f.T[0]
    tau = lF / sL
    T_track = Tu + C_TRACK * (Tb - Tu)
    pdir = PDIRS[pkey]
    fuel = "NH3" if case.get("gamma", 1.0) == 0.0 else "H2"
    phases = np.random.default_rng(SEED).uniform(0, 2 * np.pi, NMODES)

    run = os.path.join(HERE, a.outdir, a.name + a.tag)
    os.makedirs(run, exist_ok=True)
    write_pmf(os.path.join(run, "pmf.dat"), gas, f)
    meta = {"name": a.name, "case": case, "P_Pa": P, "pdir": pdir, "sL": sL, "Tu": Tu, "Tb_eq": Tb,
            "T_end_1d": float(f.T[-1]), "l_F": lF, "tau_F": tau, "T_track": T_track,
            "Lx": LX_LF * lF, "Ly": LY_LF * lF, "y0": Y0_LF * lF, "dx_fine": lF / 51.2,
            "nmodes": NMODES, "A0": a.a0 * lF, "ode_rtol": a.ode_rtol, "cfl": a.cfl, "seed": SEED, "ntau": a.ntau, "nodes": a.nodes, "mech": MECH}
    meta["sigma_1d"] = float(f.density[0] / f.density[-1])
    json.dump(meta, open(os.path.join(run, "case.json"), "w"), indent=1, default=float)

    stop = f"amr.max_step = {a.steps}\n" if a.steps else ""
    # exactly isotropic grid (PeleLMeX asserts dx == dy to 10 ulp): Ly = Lx * NY / NX (power of 2)
    Lx = float(LX_LF * lF)
    Ly = Lx * NY / NX
    # wall time: steps ~ ntau * 51.2 sigma / CFL (+ ramp), 8.5 s/step on 2 nodes measured
    # (scaling test 01/11/2026, 1.97 M cells, 10 atm), efficiency 0.86 per node doubling
    sigma = f.density[0] / f.density[-1]
    steps_est = a.ntau * 51.2 * sigma / 0.5 + 100
    if a.hours <= 0:
        a.hours = min(35.0, max(2.0, 1.3 * steps_est * 8.5 * (2.0 / a.nodes) ** 0.85 / 3600 + 0.5))
    inp = f"""# {a.name}: Lehmann 2024 dispersion-relation setup (see ../../README.md)
# l_F = {lF*1e6:.2f} um, s_L = {sL:.4f} m/s, tau_F = {tau*1e3:.4f} ms, T_track = {T_track:.1f} K
geometry.is_periodic = 1 0
geometry.coord_sys   = 0
geometry.prob_lo     = 0.0 0.0
geometry.prob_hi     = {Lx!r} {Ly!r}
peleLM.lo_bc = Interior Inflow
peleLM.hi_bc = Interior Outflow

amr.n_cell          = {NX} {NY}
amr.max_level       = 2
amr.ref_ratio       = 2 2
amr.regrid_int      = 1000000
amr.n_error_buf     = 2 2
amr.grid_eff        = 0.7
amr.blocking_factor = 16
amr.max_grid_size   = 64
amr.refinement_indicators = band
amr.band.in_box_lo = 0.0 {BAND_LF[0]*lF:.10e}
amr.band.in_box_hi = {Lx!r} {BAND_LF[1]*lF:.10e}

prob.P_mean   = {P:.6f}
prob.standoff = {Y0_LF*lF - x_fl:.10e}
prob.pertmag  = {a.a0*lF:.10e}
prob.pertlength = 0.0
prob.nmodes   = {NMODES}
prob.phases   = {" ".join(f"{p:.12f}" for p in phases)}
pmf.datafile  = "pmf.dat"

peleLM.v = 1
peleLM.use_wbar = 1
transport.use_soret = 1
peleLM.sdc_iterMax = 2
peleLM.num_init_iter = 3
peleLM.do_temporals = 1
peleLM.temporal_int = 10
peleLM.do_extremas = 1
peleLM.fuel_name = {fuel}
peleLM.divu_dt_rhoMin = 0.01
peleLM.plot_react = 0
peleLM.deltaT_iterMax = 40
peleLM.deltaT_tol = 1.0e-9
peleLM.deltaT_crashIfFailing = 0
peleLM.divu_dt_factor = 0.3

amr.stop_time = {a.ntau*tau:.10e}
{stop}amr.max_wall_time = {a.hours - 0.25:.2f}
amr.cfl = {a.cfl}
amr.dt_shrink = 0.01
amr.dt_change_max = 1.1
amr.plot_per = {0.1*tau:.10e}
amr.check_per = {1.0*tau:.10e}
amr.plot_file = plt
amr.check_file = chk

peleLM.chem_integrator = "ReactorCvode"
peleLM.use_typ_vals_chem = 1
ode.rtol = {a.ode_rtol:.1e}
ode.atol = {a.ode_atol:.1e}
cvode.solve_type = denseAJ_direct
cvode.max_order  = 4

nodal_proj.verbose = 0
nodal_proj.rtol = 1.0e-10
nodal_proj.atol = 1.0e-10
mac_proj.rtol = 1.0e-10
mac_proj.atol = 1.0e-12
"""
    open(os.path.join(run, "input.2d"), "w").write(inp)

    exe = os.path.join(HERE, "bin", f"PeleLMeX2d.gnu.MPI.{MECH}_{pdir}.ex")
    hh, mm = int(a.hours), int(round((a.hours - int(a.hours)) * 60))
    sb = f"""#!/bin/bash
#SBATCH --job-name={a.name[:40]}
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes={a.nodes}
#SBATCH --ntasks-per-node=56
#SBATCH --exclusive
#SBATCH --time={hh:02d}:{mm:02d}:00
#SBATCH --exclude=slurm-compute-h24b5-u30-svn3,slurm-compute-h21a5-u5-svn2,slurm-compute-h21c8-u5-svn4,slurm-compute-h23a5-u13-svn2,slurm-compute-h23a5-u13-svn4,slurm-compute-h23a5-u14-svn1,slurm-compute-h23a5-u14-svn3,slurm-compute-h21a5-u14-svn3
#SBATCH --output=slurm-%x-%j.out
module purge; module load foss/2024a; ulimit -s unlimited
cd "$SLURM_SUBMIT_DIR"; date
RESTART=""; last=$(ls -d chk* 2>/dev/null | sort | tail -1); [ -n "$last" ] && RESTART="amr.restart=$last"
mpirun -np $SLURM_NTASKS {exe} input.2d $RESTART > pele.log 2>&1; echo "pele exit=$?"
date
{os.path.join(HERE, "..", ".venv", "bin", "python")} {os.path.join(HERE, "extract_front.py")} . && echo "front extracted"
"""
    open(os.path.join(run, "run.sbatch"), "w").write(sb)
    print(f"{run}: sL={sL*100:.2f} cm/s l_F={lF*1e6:.1f} um tau_F={tau*1e3:.3f} ms "
          f"T_track={T_track:.0f} K pdir={pdir} fuel={fuel}")


if __name__ == "__main__":
    main()
