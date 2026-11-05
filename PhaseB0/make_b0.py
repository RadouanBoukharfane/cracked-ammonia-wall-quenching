"""Phase B0: freely propagating thermodiffusively unstable flame driven to the nonlinear (saturated)
regime, Schneider et al. 2025 Part I/II style, as the initial state of the 2D HOQ runs (Phase B).

  python make_b0.py NAME [--lx 100] [--ly 100] [--ntau 40] [--nodes 4] [--hours 24] [--a0 1e-2]

NAME is a Phase A case name (../PhaseA/cases/phaseA_cases.csv), e.g. A_g0.43_phi0.50_P10.0_T750.
Writes runs/B0_<...>/ : pmf.dat, input.2d, case.json, run.sbatch.

Setup: x periodic (lateral, Lx = 100 l_F), y inflow -> outflow (Ly = 100 l_F); base dx = l_F/5.12,
2 AMR levels on grad T -> dx_fine = l_F/20.48 (Schneider: delta_T/20). Multimode perturbation (50 modes,
a0 l_F each) to reach the nonlinear regime quickly. Active control keeps the leading edge of the
T_track isotherm at y = 30 l_F (inflow velocity follows the accelerating flame). Full plotfiles
(species, reaction rates, heat release) every tau_F, checkpoints every 2 tau_F (HOQ initial states).
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "PhaseA"))
import make_case as mc  # noqa: E402

NX = NY = 512                 # base grid -> l_F/5.12 for 100 l_F
NMODES = 50                   # lambda in [2, 100] l_F
Y_CTRL_LF = 30.0              # active-control target (leading edge), in l_F from the inlet


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--lx", type=float, default=100.0)
    ap.add_argument("--nx", type=int, default=None, help="base cells in x (default: keep dx = l_F/5.12); L_y stays NY cells")
    ap.add_argument("--ntau", type=float, default=40.0)
    ap.add_argument("--nodes", type=int, default=4)
    ap.add_argument("--hours", type=float, default=24.0)
    ap.add_argument("--a0", type=float, default=1e-2)
    ap.add_argument("--mech", default="Mei2021NH3")
    ap.add_argument("--outdir", default="runs")
    a = ap.parse_args()

    mc.MECH, mc.PELE_MECH = a.mech, os.path.join(mc.MECH_ROOT, a.mech)
    case = mc.find_case(a.name)
    gas, f, P, pkey, Tb, lF, x_fl = mc.flame(case)
    sL, Tu = f.velocity[0], f.T[0]
    tau = lF / sL
    T_track = Tu + mc.C_TRACK * (Tb - Tu)
    pdir = mc.PDIRS[pkey]
    fuel = "NH3" if case.get("gamma", 1.0) == 0.0 else "H2"
    phases = np.random.default_rng(mc.SEED + 1).uniform(0, 2 * np.pi, NMODES)
    nx = a.nx or int(round(NX * a.lx / 100.0))
    Lx = float(a.lx * lF)
    Ly = Lx * NY / nx                     # isotropic grid (PeleLMeX assert); L_y = 100 l_F for dx = l_F/5.12
    name = "B0_" + a.name[2:] + ("" if a.mech == "Mei2021NH3" else "_" + a.mech) + ("" if a.lx == 100 else f"_Lx{a.lx:.0f}")
    run = os.path.join(HERE, a.outdir, name)
    os.makedirs(run, exist_ok=True)
    mc.write_pmf(os.path.join(run, "pmf.dat"), gas, f)

    meta = {"name": name, "phaseA_case": a.name, "case": case, "P_Pa": P, "pdir": pdir, "sL": sL, "Tu": Tu,
            "Tb_eq": Tb, "l_F": lF, "tau_F": tau, "T_track": T_track, "Lx": Lx, "Ly": Ly,
            "dx_fine": Lx / nx / 4, "nmodes": NMODES, "A0": a.a0 * lF, "ntau": a.ntau, "nodes": a.nodes,
            "mech": a.mech, "sigma_1d": float(f.density[0] / f.density[-1]),
            "rho_u": float(f.density[0]), "Y_u": dict(zip(gas.species_names, map(float, f.Y[:, 0])))}
    json.dump(meta, open(os.path.join(run, "case.json"), "w"), indent=1, default=float)

    dT_tag = 0.02 * (Tb - Tu)          # adjacent-cell T jump: whole flame incl. preheat (triple-points: 15 K)
    inp = f"""# {name}: Phase B0, unstable free flame to saturation (Schneider-type), see ../../make_b0.py
# l_F = {lF*1e6:.2f} um, s_L = {sL:.4f} m/s, tau_F = {tau*1e3:.4f} ms, T_track = {T_track:.1f} K
geometry.is_periodic = 1 0
geometry.coord_sys   = 0
geometry.prob_lo     = 0.0 0.0
geometry.prob_hi     = {Lx!r} {Ly!r}
peleLM.lo_bc = Interior Inflow
peleLM.hi_bc = Interior Outflow

amr.n_cell          = {nx} {NY}
amr.max_level       = 2
amr.ref_ratio       = 2 2
amr.regrid_int      = 5
amr.n_error_buf     = 4 4
amr.grid_eff        = 0.7
amr.blocking_factor = 16
amr.max_grid_size   = 64
amr.refinement_indicators = gradT
amr.gradT.max_level = 2
amr.gradT.adjacent_difference_greater = {dT_tag:.1f}
amr.gradT.field_name = temp

prob.P_mean   = {P:.6f}
prob.standoff = {Y_CTRL_LF*lF - x_fl:.10e}
prob.pertmag  = {a.a0*lF:.10e}
prob.pertlength = 0.0
prob.nmodes   = {NMODES}
prob.phases   = {" ".join(f"{p:.12f}" for p in phases)}
pmf.datafile  = "pmf.dat"

active_control.on = 1
active_control.use_temp = 1
active_control.temperature = {T_track:.1f}
active_control.height = {Y_CTRL_LF*lF:.10e}
active_control.tau = {tau:.10e}
active_control.changeMax = {0.02*sL:.6f}
active_control.velMax = {10*sL:.6f}
active_control.flow_dir = 1
active_control.method = 1
active_control.pseudo_gravity = 1
active_control.v = 1

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
peleLM.deltaT_iterMax = 40
peleLM.deltaT_tol = 1.0e-9
peleLM.deltaT_crashIfFailing = 0
peleLM.divu_dt_factor = 0.3
amr.plot_speciesState = 1
amr.derive_plot_vars = mass_fractions
peleLM.plot_react = 1
peleLM.plot_heatRelease = 1

amr.stop_time = {a.ntau*tau:.10e}
amr.max_wall_time = {a.hours - 0.5:.2f}
amr.cfl = 0.5
amr.dt_shrink = 0.01
amr.dt_change_max = 1.1
amr.plot_per = {1.0*tau:.10e}
amr.check_per = {2.0*tau:.10e}
amr.plot_file = plt
amr.check_file = chk

peleLM.chem_integrator = "ReactorCvode"
peleLM.use_typ_vals_chem = 1
ode.rtol = 1.0e-6
ode.atol = 1.0e-5
cvode.solve_type = denseAJ_direct
cvode.max_order  = 4

nodal_proj.verbose = 0
nodal_proj.rtol = 1.0e-10
nodal_proj.atol = 1.0e-10
mac_proj.rtol = 1.0e-10
mac_proj.atol = 1.0e-12
"""
    open(os.path.join(run, "input.2d"), "w").write(inp)
    exe = os.path.join(HERE, "bin", f"PeleLMeX2d.gnu.MPI.{a.mech}_{pdir}.ex")
    hh, mm = int(a.hours), int(round((a.hours - int(a.hours)) * 60))
    sb = f"""#!/bin/bash
#SBATCH --job-name={name[:40]}
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
"""
    open(os.path.join(run, "run.sbatch"), "w").write(sb)
    print(f"{run}: sL={sL*100:.2f} cm/s l_F={lF*1e6:.1f} um tau_F={tau*1e3:.3f} ms Lx={Lx*1e3:.2f} mm "
          f"dx_fine={Lx/nx/4*1e6:.2f} um, gradT tag {dT_tag:.0f} K, pdir={pdir}")


if __name__ == "__main__":
    main()
