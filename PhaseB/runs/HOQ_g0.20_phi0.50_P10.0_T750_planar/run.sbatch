#!/bin/bash
#SBATCH --job-name=HOQ_g0.20_phi0.50_P10.0_T750_planar
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --exclusive
#SBATCH --time=12:00:00
#SBATCH --exclude=slurm-compute-h24b5-u30-svn3,slurm-compute-h21a5-u5-svn2,slurm-compute-h21c8-u5-svn4,slurm-compute-h23a5-u13-svn2,slurm-compute-h23a5-u13-svn4,slurm-compute-h23a5-u14-svn1,slurm-compute-h23a5-u14-svn3
#SBATCH --output=slurm-%x-%j.out
module purge; module load foss/2024a; ulimit -s unlimited
cd "$SLURM_SUBMIT_DIR"; date
RESTART=""; last=$(ls -d chk* 2>/dev/null | sort | tail -1); [ -n "$last" ] && RESTART="amr.restart=$last"
mpirun -np $SLURM_NTASKS ${FQ_ROOT}/PhaseB/bin/PeleLMeX2d.gnu.MPI.Mei2021NH3_10_000atm.ex input.2d $RESTART > pele.log 2>&1; echo "pele exit=$?"
date
