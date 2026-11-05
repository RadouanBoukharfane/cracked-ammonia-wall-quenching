"""Phase B: 2D head-on quenching (HOQ) against an isothermal inert wall (Schneider et al. 2025, Part I).

Unstable HOQ, initialised from a saturated Phase B0 state (plotfile with mass fractions):
  python make_hoq.py unstable --b0 ../PhaseB0/runs/B0_g0.43_phi0.50_P10.0_T750 \
                              --plt ../PhaseB0/runs/B0_.../species/from_chk04404/plt_sp04405 [--ntau 30]
Planar reference (narrow domain, TDI suppressed; equivalent to the 1D HOQ), initialised from the 1D flame:
  python make_hoq.py planar --b0 ../PhaseB0/runs/B0_g0.43_phi0.50_P10.0_T750 [--dist 30] [--ntau 45]

Wall at y = 0 (no-slip, T_w = T_u, zero species flux), outflow at y = L_y, x periodic. Wall frame: the
unburnt gas is at rest; the B0 inflow velocity V_in at the plotfile time is removed by patchFlowVariables
(prob.v_shift), the 1D flame velocity s_L is removed in the planar case. Grid as B0 (dx_fine = l_F/20.48),
plus the first 3 l_F above the wall always on the finest level. Plotfiles with species, reaction rates and heat
release every 0.5 tau_F (wall heat flux, quenching distance, emissions), checkpoints every 5 tau_F.
"""
import argparse
import json
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NXB = 512             # B0 base grid (L_x = 100 l_F)
NX_PLANAR = 16        # planar reference: 16 base cells = 3.125 l_F wide


def vin_at(b0, t):
    h = np.loadtxt(os.path.join(b0, "AC_History.dat"))
    return float(np.interp(t, h[:, 1], h[:, 2]))


def plt_time(plt):
    with open(os.path.join(plt, "Header")) as f:
        lines = f.read().split("\n")
    n = int(lines[1])
    return float(lines[2 + n + 1])     # after version, nvars, var names, dim -> time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["unstable", "planar"])
    ap.add_argument("--b0", required=True)
    ap.add_argument("--plt", help="B0 plotfile with Y(...) fields (unstable)")
    ap.add_argument("--dist", type=float, default=30.0, help="planar: initial flame distance from the wall, l_F")
    ap.add_argument("--ntau", type=float, default=None)
    ap.add_argument("--nodes", type=int, default=None)
    ap.add_argument("--hours", type=float, default=None)
    ap.add_argument("--tw", type=float, default=None, help="wall temperature [K] (default T_w = T_u)")
    ap.add_argument("--levels", type=int, default=2, help="AMR levels above the base grid (3: grid-convergence check, dx_fine halved)")
    a = ap.parse_args()

    b0 = os.path.abspath(a.b0)
    m = json.load(open(os.path.join(b0, "case.json")))
    lF, sL, tau, P, Tu = m["l_F"], m["sL"], m["tau_F"], m["P_Pa"], m["Tu"]
    inp_b0 = open(os.path.join(b0, "input.2d")).read()
    dT_tag = float(re.search(r"amr.gradT.adjacent_difference_greater = (\S+)", inp_b0).group(1))
    base = os.path.basename(b0).replace("B0_", "")
    Tw = a.tw if a.tw is not None else Tu
    tw_tag = "" if a.tw is None else f"_Tw{a.tw:.0f}"

    if a.mode == "unstable":
        plt = os.path.abspath(a.plt)
        t0 = plt_time(plt)
        vshift = vin_at(b0, t0)
        name = f"HOQ_{base}_t{t0 / tau:04.1f}{tw_tag}"
        nx, ny, Lx = NXB, NXB, m["Lx"]
        ntau = a.ntau or 30.0
        nodes, hours = a.nodes or 4, a.hours or 12.0
        init = f"""amr.initDataPlt = {plt}
peleLM.initDataPlt_reset_time = 1
peleLM.initDataPlt_patch_flow_variables = 1
prob.v_shift = {vshift:.10e}
prob.pertmag = 0.0
prob.standoff = 0.0"""
        meta_extra = {"init_plt": plt, "t0_b0_tauF": t0 / tau, "v_shift": vshift}
    else:
        name = f"HOQ_{base}_planar{tw_tag}" + ("" if a.dist == 30.0 else f"_d{a.dist:.0f}") + ("" if a.levels == 2 else f"_L{a.levels}")
        nx, ny = NX_PLANAR, NXB
        Lx = m["Lx"] * NX_PLANAR / NXB
        ntau = a.ntau or 45.0
        nodes, hours = a.nodes or 1, a.hours or 12.0
        x_fl = float(re.search(r"prob.standoff = (\S+)", inp_b0).group(1))     # standoff of B0 = 30 l_F - x_fl
        x_fl = 30.0 * lF - x_fl
        init = f"""prob.v_shift = {sL:.10e}
prob.pertmag = 0.0
prob.nmodes = 0
prob.standoff = {a.dist * lF - x_fl:.10e}"""
        meta_extra = {"dist_lF": a.dist, "v_shift": sL}
    Ly = Lx * ny / nx
    # after quenching the flow slows down and the convective dt grows until the divU check aborts:
    # cap dt at 1.5 x the pre-quench convective dt (CFL dx_fine / ((sigma - 1) s_L))
    max_dt = 1.5 * 0.5 * (Lx / nx / 2**a.levels) / ((m["sigma_1d"] - 1.0) * sL)
    run = os.path.join(HERE, "runs", name)
    os.makedirs(run, exist_ok=True)
    os.system(f"cp {os.path.join(b0, 'pmf.dat')} {run}/")
    meta = dict(m, name=name, mode=a.mode, Lx=Lx, Ly=Ly, ntau=ntau, nodes=nodes, T_wall=Tw, levels=a.levels,
                dx_fine=Lx / nx / 2**a.levels, **meta_extra)
    json.dump(meta, open(os.path.join(run, "case.json"), "w"), indent=1, default=float)

    inp = f"""# {name}: Phase B head-on quenching ({a.mode}), see ../../make_hoq.py
# l_F = {lF*1e6:.2f} um, s_L = {sL:.4f} m/s, tau_F = {tau*1e3:.4f} ms, T_wall = {Tw:.1f} K
geometry.is_periodic = 1 0
geometry.coord_sys   = 0
geometry.prob_lo     = 0.0 0.0
geometry.prob_hi     = {Lx!r} {Ly!r}
peleLM.lo_bc = Interior NoSlipWallIsotherm
peleLM.hi_bc = Interior Outflow

amr.n_cell          = {nx} {ny}
amr.max_level       = {a.levels}
amr.ref_ratio       = {" ".join(["2"] * a.levels)}
amr.regrid_int      = 5
amr.n_error_buf     = 4 4
amr.grid_eff        = 0.7
amr.blocking_factor = 16
amr.max_grid_size   = 64
amr.refinement_indicators = gradT wall
amr.gradT.max_level = {a.levels}
amr.gradT.adjacent_difference_greater = {dT_tag}
amr.gradT.field_name = temp
amr.wall.in_box_lo = 0.0 0.0
amr.wall.in_box_hi = {Lx!r} {3.0*lF:.10e}

prob.P_mean   = {P:.6f}
prob.T_wall   = {Tw:.6f}
prob.wall_lo  = 1
{init}
pmf.datafile  = "pmf.dat"

peleLM.v = 1
peleLM.use_wbar = 1
transport.use_soret = 1
peleLM.sdc_iterMax = 2
peleLM.num_init_iter = 3
peleLM.do_temporals = 1
peleLM.temporal_int = 10
peleLM.do_extremas = 1
peleLM.fuel_name = H2
peleLM.divu_dt_rhoMin = 0.01
peleLM.deltaT_iterMax = 40
peleLM.deltaT_tol = 1.0e-9
peleLM.deltaT_crashIfFailing = 0
peleLM.divu_dt_factor = 0.3
peleLM.plot_react = 1
peleLM.plot_heatRelease = 1
amr.plot_speciesState = 1
amr.derive_plot_vars = mass_fractions

amr.stop_time = {ntau*tau:.10e}
amr.max_wall_time = {hours - 0.5:.2f}
amr.cfl = 0.5
amr.max_dt = {max_dt:.6e}
amr.dt_shrink = 0.01
amr.dt_change_max = 1.1
amr.plot_per = {0.5*tau:.10e}
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
    exe = os.path.join(HERE, "bin", f"PeleLMeX2d.gnu.MPI.{m['mech']}_{m['pdir']}.ex")
    hh, mm = int(hours), int(round((hours - int(hours)) * 60))
    open(os.path.join(run, "run.sbatch"), "w").write(f"""#!/bin/bash
#SBATCH --job-name={name[:40]}
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes={nodes}
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
""")
    print(f"{run}: Lx = {Lx/lF:.2f} l_F, v_shift = {meta['v_shift']:.4f} m/s ({meta['v_shift']/sL:.2f} s_L), "
          f"{ntau} tau_F, {nodes} nodes")


if __name__ == "__main__":
    main()
