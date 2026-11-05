#!/bin/bash
#SBATCH --job-name=A_g0.60_phi0.40_P01.0_T750
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=56
#SBATCH --exclusive
#SBATCH --time=07:44:00
#SBATCH --exclude=slurm-compute-h24b5-u30-svn3,slurm-compute-h21a5-u5-svn2,slurm-compute-h21c8-u5-svn4
#SBATCH --output=slurm-%x-%j.out
module purge; module load foss/2024a; ulimit -s unlimited
cd "$SLURM_SUBMIT_DIR"; date
RESTART=""; last=$(ls -d chk* 2>/dev/null | sort | tail -1); [ -n "$last" ] && RESTART="amr.restart=$last"
mpirun -np $SLURM_NTASKS ${FQ_ROOT}/PhaseA/bin/PeleLMeX2d.gnu.MPI.Mei2021NH3_1_000atm.ex input.2d $RESTART > pele.log 2>&1; echo "pele exit=$?"
date
${FQ_ROOT}/PhaseA/../.venv/bin/python ${FQ_ROOT}/PhaseA/extract_front.py . && echo "front extracted"
