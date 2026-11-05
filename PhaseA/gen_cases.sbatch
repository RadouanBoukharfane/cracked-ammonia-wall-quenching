#!/bin/bash
#SBATCH --job-name=A_gen
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --exclude=slurm-compute-h24b5-u30-svn3,slurm-compute-h21a5-u5-svn2,slurm-compute-h21c8-u5-svn4,slurm-compute-h23a5-u13-svn2,slurm-compute-h23a5-u13-svn4,slurm-compute-h23a5-u14-svn1,slurm-compute-h23a5-u14-svn3,slurm-compute-h21a5-u14-svn3
#SBATCH --output=slurm-%x-%j.out
# usage: sbatch gen_cases.sbatch LISTFILE [make_case options, e.g. --cfl 0.3]  -> runs/<name>/, 15 min max per case
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
xargs -a "$1" -P 56 -I{} timeout 900 ../.venv/bin/python make_case.py {} "${@:2}" 2>&1 | grep -E "sL=|Error|rror"
echo done
