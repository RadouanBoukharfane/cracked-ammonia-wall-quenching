#!/bin/bash
#SBATCH --job-name=B0_post
#SBATCH --account=msdaphyik-frgj4qkkjni-premium-cpu
#SBATCH --qos=premium-cpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --exclude=slurm-compute-h24b5-u30-svn3,slurm-compute-h21a5-u5-svn2,slurm-compute-h21c8-u5-svn4
#SBATCH --output=slurm-%x-%j.out
# usage: sbatch post_b0.sbatch RUN_DIR
cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1
../.venv/bin/python post_b0.py "$1" 28 > "$1/post_b0.log" 2>&1; echo "exit=$?"
