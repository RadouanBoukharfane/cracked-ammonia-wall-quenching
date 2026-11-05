#!/bin/bash
#SBATCH --job-name=B0_build
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --exclude=slurm-compute-h24b5-u30-svn3,slurm-compute-h21a5-u5-svn2,slurm-compute-h21c8-u5-svn4
#SBATCH --output=slurm-%x-%j.out
module purge; module load foss/2024a
cd "$SLURM_SUBMIT_DIR"
make COMP=gnu MECH=Mei2021NH3 PDIR=10_000atm -j 56 > build.log 2>&1 && cp PeleLMeX2d.gnu.MPI.Mei2021NH3_10_000atm.ex ../../bin/ && echo "build ok"
